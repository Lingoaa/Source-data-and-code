rm(list=ls())
setwd("E:/AOM_tumor/filter_data/cellchat")


##加载R包
{
  library(CellChat)
  library(Seurat)
  library(patchwork)
  library(ggplot2)
  library(openxlsx)
}

# Part1 数据输入和处理以及CellChat对象的初始化
## 导入数据
imm <- readRDS("colon_imm_seurat.rds")

bc <- readRDS("5bc_seurat.rds")

tc <- readRDS("TC5_seurat.rds")


rm(NT_BT)

#合并数据
all <- merge(imm, y = c(bc, tc), add.cell.ids = c("4K", "8K"), project = "PBMC12K",merge.data = TRUE)

table(NT_BT$orig.ident)
table(NT_BT$group)
table(NT_BT$hash.ID)
table(NT_BT$cell_type)

saveRDS(NT_BT, file = "NT_BT_seurat.rds")

#提取normal
N_BT <- subset(x = NT_BT, subset = group == "normal")

saveRDS(N_BT, file = "N_BT_seurat.rds")

#提取tumor
T_BT <- subset(x = NT_BT, subset = group == "tumor")

saveRDS(T_BT, file = "T_BT_seurat.rds")


####重新导入数据分析
rm(list=ls())
setwd("E:/spleen_alldata/cellchat/tumor_BT")

T_BT <-  readRDS("T_BT_seurat.rds")
T_BT
table(T_BT$orig.ident)
table(T_BT$cell_type)

#分析wt
tb_ko <- T_BT[,T_BT@meta.data[["orig.ident"]]=='T-ko']
exp <- tb_ko@assays$RNA@data
meta <- tb_ko@meta.data

index1 <- which(colnames(exp) %in% rownames(meta))
data.input <- exp[,index1]
dim(data.input)

## 创建CellChat对象
cellchat <- createCellChat(object = data.input,
                           meta = meta,
                           group.by = "cell_type")

## 导入配体受体相互作用数据库
CellChatDB <- CellChatDB.mouse  #CellChatDB.human,CellChatDB.mouse
showDatabaseCategory(CellChatDB)
dplyr::glimpse(CellChatDB$interaction)# Show the structure of the database
# use a subset of CellChatDB for cell-cell communication analysis
#CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
cellchat@DB <- CellChatDB

## 数据预处理
cellchat <- subsetData(cellchat)# subset the expression data of signaling genes for saving computation cost
#future::plan("multisession", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
cellchat <- projectData(cellchat, PPI.human)# project gene expression data onto PPI network (optional)
#推测的每个配体-受体对的细胞间通信网络和每个信号通路分别存储在“net”和“netP”槽中。

# Part2 细胞间通讯网络推断
## 计算通讯概率并推断网络
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)

## 提取推断得到的细胞通讯网络
df.net <- subsetCommunication(cellchat)
#df.net <- subsetCommunication(all_cellchat, signaling = 'WNT')#netp返回通路级别数据/sources.use返回信号起点终点/signaling 返回特定信号通路

## 在信号通路水平上推断细胞间通讯
cellchat <- computeCommunProbPathway(cellchat)

## 计算整合网络
cellchat <- aggregateNet(cellchat)
groupSize <- as.numeric(table(cellchat@idents))

mat <- cellchat@net$weight
par(mfrow = c(3,4), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

# 保存CellChat对象
cellchat@netP$pathways

saveRDS(cellchat, file = "cellchat_T_BT_wt.rds")

rm(list=ls())


## 导入数据
cellchat_wt <- readRDS(file = "cellchat_T_BT_wt.rds")
cellchat_ko <- readRDS(file = "cellchat_T_BT_ko.rds")


#合并cellchat对象
cco.list <- list(T_wt=cellchat_wt, T_ko=cellchat_ko)

cellchat <- mergeCellChat(cco.list, add.names = names(cco.list), cell.prefix = TRUE)


#所有细胞群总体观：通讯数量与强度对比
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "count")
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
p <- gg1 + gg2
p

ggsave("Overview_number_strength.png", p, width = 6, height = 4)


#数量与强度差异网络图
par(mfrow = c(1,2))
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

#数量与强度差异热图
par(mfrow = c(1,1))
h1 <- netVisual_heatmap(cellchat)
h2 <- netVisual_heatmap(cellchat, measure = "weight")
h1+h2

#细胞互作数量对比网络图
par(mfrow = c(1,2))
weight.max <- getMaxWeight(cco.list, attribute = c("idents","count"))
for (i in 1:length(cco.list)) {
  netVisual_circle(cco.list[[i]]@net$count, weight.scale = T, label.edge= F, 
                   edge.weight.max = weight.max[2], edge.width.max = 12, 
                   title.name = paste0("Number of interactions - ", names(cco.list)[i]))
}


#保守和特异性信号通路的识别与可视化
## 通路信号强度对比分析
gg1 <- rankNet(cellchat, mode = "comparison", stacked = T, do.stat = TRUE)
gg2 <- rankNet(cellchat, mode = "comparison", stacked = F, do.stat = TRUE)
p <- gg1 + gg2
p
ggsave("Compare_pathway_strengh.png", p, width = 10, height = 6)


#气泡图展示所有配体受体对的差异
levels(cellchat@idents$joint)
p <- netVisual_bubble(cellchat, sources.use = c(3), targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14),  comparison = c(1, 2), angle.x = 45)
p
ggsave("Compare_bubble.png", p, width = 12, height = 8)

#气泡图展示上调或下调的配体受体对
p1 <- netVisual_bubble(cellchat, sources.use = c(3), targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14), comparison = c(1, 2), 
                       max.dataset = 2, title.name = "Increased signaling in T_ko", angle.x = 45, remove.isolate = T)
p2 <- netVisual_bubble(cellchat, sources.use = c(4,5), targets.use = c(1,2,3,6), comparison = c(1, 2), 
                       max.dataset = 1, title.name = "Decreased signaling in TIL", angle.x = 45, remove.isolate = T)
pc <- p1 + p2
ggsave("Increased signaling in T_ko.png", p1, width = 12, height = 7)


#####单通路展示
## 查看重要通信的信号通路
cellchat_ko@netP$pathways
cellchat_wt@netP$pathways

color <- c("#c04327","#dd6f56","#e79a88","#c491a5","#8c8ec1","#2a4986",
           "#8189b0")

color <- c('#800020','#9A3671','#E1C392','#25BCCD','#95A238','#FADA5E',
           '#B05923',"#c04327","#dd6f56","#e79a88","#c491a5","#8c8ec1","#2a4986",
           "#8189b0") 
source <- c('Cd8+ NKT',
            'Effector Memory Cd4+ TC',
            'Follicular BC',
            'GC BC',
            'IFI+ BC',
            'MZ BC',
            'Naive Cd4+ TC',
            'Naive Cd8+ TC',
            'Plasma cell',
            'Plasmblast',
            'Pro/Pre BC',
            'Tfh',
            'Transitional BC',
            'Tregs')
target <- c('Cd8+ NKT',
            'Effector Memory Cd4+ TC',
            'Follicular BC',
            'GC BC',
            'IFI+ BC',
            'MZ BC',
            'Naive Cd4+ TC',
            'Naive Cd8+ TC',
            'Plasma cell',
            'Plasmblast',
            'Pro/Pre BC',
            'Tfh',
            'Transitional BC',
            'Tregs')


## 挑选用于可视化的信号通路
CD40 <- c("CD40")
MHCI <- c("MHC-I")

CD22 <- c("CD22")
CD45 <- c("CD45")
GALECTIN <- c("GALECTIN")
PECAM1 <- c("PECAM1")
SELPLG <- c("SELPLG")
LAMININ <- c("LAMININ")
SEMA4 <- c("SEMA4")
TGFb <- c("TGFb")
CD80 <- c("CD80")
SEMA7 <- c("SEMA7")

## 弦图
#CD52
netVisual_chord_gene(cellchat_ko, big.gap = 20,small.gap = 5,
                     targets.use = target,transparency = 0.5,
                     signaling = MHCI,show.legend = T,color.use = color,
                     title.name = 'MHCI signaling pathway interaction in ko')

netVisual_chord_gene(cellchat_wt, big.gap = 20,small.gap = 5,
                     targets.use = target,transparency = 0.5,
                     signaling = MHCI,show.legend = T,color.use = color,
                     title.name = 'MHCI signaling pathway interaction in wt')



## 热图
#CD52
netVisual_heatmap(cellchat_ko, color.heatmap = 'Reds',color.use = color,
                  signaling = MHCI,targets.use = target,
                  sources.use = source,
                  title.name = 'MHC-I signaling pathway interaction in ko')

p2 = netVisual_heatmap(cellchat_wt, color.heatmap = 'Reds',color.use = color,
                  signaling = MHCI,targets.use = target,
                  sources.use = source,
                  title.name = 'MHC-I signaling pathway interaction in wt')
p2

ggsave("MHCI_wt_heatmap.pdf",p2,width = 12, height = 7)

## 点图

netVisual_bubble(cellchat_ko,
                 signaling = "MHC-I",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'MHC-I signaling pathway in ko')
netVisual_bubble(cellchat_wt,
                 signaling = "MHC-I",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'MHC-I signaling pathway in wt')

#MIF
netVisual_bubble(cellchat_ko,
                 signaling = "MIF",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'MIF signaling pathway in ko')
netVisual_bubble(cellchat_wt,
                 signaling = "MIF",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'MIF signaling pathway in ko')

#CD22
netVisual_bubble(cellchat_ko,
                 signaling = "CD22",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'CD22 signaling pathway in ko')
netVisual_bubble(cellchat_wt,
                 signaling = "CD22",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'wt-CD22 signaling pathway in wt')
#CD45
netVisual_bubble(cellchat_ko,
                 signaling = "CD45",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'CD45 signaling pathway in ko')
netVisual_bubble(cellchat_wt,
                 signaling = "CD45",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'wt-CD45 signaling pathway in wt')
#GALECTIN
netVisual_bubble(cellchat_ko,
                 signaling = "GALECTIN",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'GALECTIN signaling pathway in ko')
netVisual_bubble(cellchat_wt,
                 signaling = "GALECTIN",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'GALECTIN signaling pathway in wt')
#PECAM1
netVisual_bubble(cellchat_ko,
                 signaling = "PECAM1",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'PECAM1 signaling pathway in ko')
netVisual_bubble(cellchat_wt,
                 signaling = "PECAM1",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'PECAM1 signaling pathway in wt')

##SELPLG
netVisual_bubble(cellchat_ko,
                 signaling = "SELPLG",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'SELPLG signaling pathway in ko')
netVisual_bubble(cellchat_wt,
                 signaling = "SELPLG",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'SELPLG signaling pathway in wt')

#LAMININ
netVisual_bubble(cellchat_ko,
                 signaling = "LAMININ",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'LAMININ signaling pathway in ko')
netVisual_bubble(cellchat_wt,
                 signaling = "LAMININ",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'LAMININ signaling pathway in wt')
##SEMA4
netVisual_bubble(cellchat_ko,
                 signaling = "SEMA4",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'SEMA4 signaling pathway in ko')
netVisual_bubble(cellchat_wt,
                 signaling = "SEMA4",targets.use =  target, 
                 angle.x = 45, remove.isolate =T,
                 title.name = 'SEMA4 signaling pathway in wt')



#导出prob数据用于标注在热图上
KO_MHCI <- cellchat_ko@netP[["prob"]][,,'MHC-I'] 
WT_MHCI <- cellchat_wt@netP[["prob"]][,,'MHC-I'] 

KO_MIF <- cellchat_ko@netP[["prob"]][,,'MIF'] 
WT_MIF <- cellchat_wt@netP[["prob"]][,,'MIF'] 

KO_CD22 <- cellchat_ko@netP[["prob"]][,,'CD22'] 
WT_CD22 <- cellchat_wt@netP[["prob"]][,,'CD22'] 

KO_CD45 <- cellchat_ko@netP[["prob"]][,,'CD45'] 
WT_CD45 <- cellchat_wt@netP[["prob"]][,,'CD45'] 

KO_GALECTIN <- cellchat_ko@netP[["prob"]][,,'GALECTIN'] 
WT_GALECTIN <- cellchat_wt@netP[["prob"]][,,'GALECTIN']

KO_PECAM1 <- cellchat_ko@netP[["prob"]][,,'PECAM1'] 
WT_PECAM1 <- cellchat_wt@netP[["prob"]][,,'PECAM1']

KO_SELPLG <- cellchat_ko@netP[["prob"]][,,'SELPLG'] 
WT_SELPLG <- cellchat_wt@netP[["prob"]][,,'SELPLG']

list_of_datasets <- list('ko_MHCI' = KO_MHCI, 
                         'wt_MHCI' = WT_MHCI)

list_of_datasets <- list('ko_MHCI' = KO_CD52, 
                         'wt_MHCI' = WT_CD52,
                         'ko_MIF' = KO_MIF,
                         'wt_MIF' = WT_MIF,
                         'ko_CD22' = KO_CD22,
                         'wt_CD22' = WT_CD22,
                         'ko_CD45' = KO_CD45,
                         'wt_CD45' = WT_CD45,
                         'ko_GALECTIN' = KO_GALECTIN,
                         'wt_GALECTIN' = WT_GALECTIN,
                         'ko_PECAM1' = KO_PECAM1,
                         'wt_PECAM1' = WT_PECAM1,
                         'ko_SELPLG' = KO_SELPLG,
                         'wt_SELPLG' = WT_SELPLG)

write.xlsx(list_of_datasets, "MHCI_cellchat_pathways_prob.xlsx",rowNames=TRUE)
