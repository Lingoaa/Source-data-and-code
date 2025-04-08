setwd("E:/spleen_alldata/NT_spleen_TC")

rm(list=ls())

##加载R包
{
  library(densityClust)
  library(scran)
  library(Seurat)
  library(tidyverse)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
  library(DoubletFinder)
  library(BiocSingular)
  library(sctransform)
  library(glmGamPoi)
  library(scater)
  library(scDblFinder)
  library(gridExtra)
  library(SingleR)
  library(cowplot)
  library(harmony)
  library(devtools)
  library(tidyverse)
  library(scCATCH)
  library(mindr)
  
}

set.seed(123) 

###导入数据
#######normal
##导入数据
exp <- readRDS("NT_spleen_TC_exp.rds")

meta <- read.csv("NT_spleen_TC_meta.csv", row.names = 1)


##初始化Seurat对象
TC <- CreateSeuratObject(counts =exp , project = "T_cell", min.cells = 3) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = TC@var.genes, npcs = 20, verbose = FALSE)


##添加分组信息
#分组信息
orig.ident <- meta$orig.ident
names(orig.ident) <- colnames(x = TC)
TC <- AddMetaData(
  object = TC,
  metadata = orig.ident,
  col.name = "orig.ident"
)
table(TC@meta.data$orig.ident)

#样本信息
hash.ID <- meta$hash.ID
names(hash.ID) <- colnames(x = TC)
TC <- AddMetaData(
  object = TC,
  metadata = hash.ID,
  col.name = "hash.ID"
)
table(TC@meta.data$hash.ID)

#group信息
group <- meta$group
names(group) <- colnames(x = TC)
TC <- AddMetaData(
  object = TC,
  metadata = group,
  col.name = "group"
)
table(TC@meta.data$group)

save(TC, file = "SBC.RData")
load("SBC.RData")

table(TC@meta.data$orig.ident)
table(TC@meta.data$group)


p1 <- DimPlot(object = TC, reduction = "pca", pt.size = .1,group.by = "hash.ID")
p2 <- VlnPlot(object = TC, features = "PC_1", pt.size = .1,group.by = "hash.ID")
plot_grid(p1,p2)
p = plot_grid(p1,p2)
ggsave("sbc_plot.png", p ,width=9 ,height=6)


##Run Harmony
TC <- TC %>%
  RunHarmony("hash.ID", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(TC, 'harmony')
harmony_embeddings[1:5, 1:5]

p1 <- DimPlot(object = TC, reduction = "harmony", pt.size = .1,group.by = "hash.ID")
p2 <- VlnPlot(object = TC, features = "harmony_1", pt.size = .1,group.by = "hash.ID")
plot_grid(p1,p2)

p = plot_grid(p1,p2)

ggsave("TC_plot_harmony.png", p ,width=9 ,height=6)

##Downstream analysis
TC <- FindNeighbors(TC, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.2)

TC <- RunUMAP(TC, reduction = "harmony", dims = 1:20)

DimPlot(TC, reduction = "umap",label = T,pt.size = 0.6) 

##找高变基因
markers <- FindAllMarkers(object = TC, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)   

all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)

top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

write.csv(top20, "NT_TC_markers_.csv")

write.csv(all.markers, "TC_markers.csv")


##DotPlot
genes_to_check = c("Cd3d","Cd4","Ly6c1","Igfbp4","Ccr7","Sugct","Selenop",
                   "Cd8b1","Cd8a","Dapl1","Rgcc","Tagln2",
                   "Ly6c2","Ccl5","Cd7","Nkg7","Ctsw",
                   "Lgals1", "S100a4", "Itgb1",
                   "Foxp3","Ctla4","Izumo1r","Ikzf2",
                   "Cd40lg","Icos","Il4", "Gata3", "Bcl6")
p <- DotPlot(TC, features = genes_to_check,assay='RNA' ,cols = c(low="white",high="darkred"))
p


new.cluster.ids <- c("0" = "Naive Cd4+ TC",
                     "1" = "Naive Cd8+ TC",
                     "2" = "Cd8+ NKT",
                     "3" = "Effector Memory Cd4+ TC",
                     "4" = "Tregs",
                     "5" = "Tfh",
                     "6" = "Cd8+ NKT",
                     "7" = "other",
                     "8" = "other")

names(new.cluster.ids) <- levels(TC)
TC <- RenameIdents(TC, new.cluster.ids)

#查看细胞数
table(TC@meta.data$orig.ident,TC@active.ident)

ident <- data.frame(TC@active.ident)

head(ident)

TC@meta.data$cell_type <- ident$TC.active.ident

table(TC@meta.data$cell_type)
table(TC@meta.data$hash.ID, TC@meta.data$cell_type)


##去除非T细胞
TC_scRNA <- subset(x=TC, idents =c("Naive Cd4+ TC","Naive Cd8+ TC","Tregs","Cd8+ NKT","Effector Memory Cd4+ TC","Tfh"))

TC_scRNA

table(TC_scRNA@meta.data$orig.ident)

table(TC_scRNA@meta.data$orig.ident,TC_scRNA@active.ident)

rm(scRNA)

#FeaturePlot
library(viridis)
pal <- viridis(n = 15, option = "C", direction = -1)

p1 = FeaturePlot(TC_scRNA,features = c("Ly6c1","Igfbp4","Ccr7","Selenop",
                                       "Dapl1","Rgcc",
                                       "Ly6c2","Ccl5","Nkg7","Ctsw",
                                       "Lgals1", "S100a4", "Itgb1",
                                       "Foxp3","Ctla4",
                                       "Cd40lg","Icos","Il4"),cols = pal, order = T, ncol =6)

ggsave("marker_FeaturePlot.png", p1, width=18,height=7)

##ggplot绘图
##提取前两个主成分数据
# extact PC ranges
pc12 <- Embeddings(object = TC_scRNA,reduction = 'umap') %>%
  data.frame()

# check
head(pc12,3)

##构造坐标轴需要的标签和位置信息
# get botomn-left coord
lower <- floor(min(min(pc12$UMAP_1),min(pc12$UMAP_2))) - 2

# get relative line length
linelen <- abs(0.3*lower) + lower

# mid point
mid <- abs(0.3*lower)/2 + lower

# axies data
axes <- data.frame(x = c(lower,lower,lower,linelen),y = c(lower,linelen,lower,lower),
                   group = c(1,1,2,2),
                   label = rep(c('UMAP_2','UMAP_1'),each = 2))

# axies label
label <- data.frame(lab = c('UMAP_2','UMAP_1'),angle = c(90,0),
                    x = c(lower - 3,mid),y = c(mid,lower - 2.5))

##可视化
# plot
color <- c("#eb998b", "#4DBBD5B2","#A9D179","#84CAC0","#6a51a3","#a14462")

plot1 <- DimPlot(TC_scRNA, reduction = 'umap', label = T,repel = TRUE,
                 pt.size = .6, cols = color,label.size = 5) +
  NoAxes() + 
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab))

plot1

color <- c("darkred","#a14462","#F5AE6B","#6a51a3","#4DBBD5B2","#b0d45d","#7E6148B2","#0072B2","#356d67","#9e9ac8")
color <- c("#35212e","#562e3c","#a14462","#eb998b",
"#fddbc8","#42465c","#356d67","#4c9568",
"#7fb961","#b0d45d")

plot2 <- DimPlot(TC_scRNA, reduction = 'umap', label = F,repel = TRUE,
                 pt.size = .5, cols = color,label.size = 4,group.by = 'hash.ID') +
  NoAxes() +
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab))

plot2

ggsave('imm_group_v2.png',plot2,width = 10,height = 8)


plot3 <- DimPlot(TC_scRNA, reduction = 'umap', label = F,repel = TRUE,
                 pt.size = .8, cols = color,label.size = 4,split.by  = "orig.ident",ncol = 2) +
  NoAxes() +
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab))

plot3


####提取normal/tumor
table(TC_scRNA@meta.data$group)

#提取normal
N_TC <- subset(x = TC_scRNA, subset = group == "normal")
#N_TC_exp <- as.matrix(N_TC@assays$RNA@counts)
#saveRDS(N_TC_exp, file = "N_TC_exp.rds")
#write.csv(N_TC@meta.data, "N_TC_meta.csv")

saveRDS(N_TC, file = "N_TC_seurat.rds")

rm(N_TC_exp)

#提取tumor
T_TC <- subset(x = TC_scRNA, subset = group == "tumor")

#T_TC_exp <- as.matrix(T_TC@assays$RNA@counts)
#saveRDS(T_TC_exp, file = "T_TC_exp.rds")
#write.csv(T_TC@meta.data, "T_TC_meta.csv")

saveRDS(T_TC, file = "T_TC_seurat.rds")

rm(T_TC_exp)


####不同分组各细胞比例统计图

##加载R包
options(stringsAsFactors = F)
{
  library(densityClust)
  library(scran)
  library(Seurat)
  library(tidyverse)
  library(dplyr)
  library(patchwork)
  library(ggplot2)
  library(DoubletFinder)
  library(BiocSingular)
  library(sctransform)
  library(glmGamPoi)
  library(scater)
  library(scDblFinder)
  library(gridExtra)
  library(SingleR)
  library(cowplot)
  library(harmony)
  library(devtools)
  library(tidyverse)
  library(scCATCH)
  library(mindr)
  library(Seurat)
  library(ggplot2)
  library(clustree)
  library(cowplot)
  library(dplyr)
  library(tidyr)# 使用的gather & spread
  library(reshape2) # 使用的函数 melt & dcast 
  library (gplots) 
  
}


#加载数据
load("")
TC <- readRDS("TC_seurat_harmony.rds")
TC

phe=T_TC@meta.data
colnames(phe)

tb=table(phe$cell_type,
         phe$hash.ID)

head(tb)

#统计细胞数量
balloonplot(tb, main ="Tumor-spleen T cells", xlab ="celltype", ylab="sample",
            label = T, show.margins = T)


bar_data <- as.data.frame(tb)

bar_per <- bar_data %>% 
  group_by(Var2) %>%
  mutate(sum(Freq)) %>%
  mutate(percent = Freq / `sum(Freq)`)
head(bar_per) 
write.csv(bar_per,file = "celltype_by_group_percent.csv")
col =c("#BD6263","#8EA325","#3FA116","#CE2820","#9265C1","#885649","#2EBEBE")

col <- c("#a14462","#eb998b","#9265C1","#2EBEBE","#b0d45d",
           "#356d67","#42465c",
          "#b0d45d", "#7fb961","#562e3c")

ggplot(bar_per, aes(x = percent, y = Var2)) +
  geom_bar(aes(fill = Var1) , stat = "identity") + coord_flip() +
  theme(axis.ticks = element_line(linetype = "blank"),
        legend.position = "right",
        panel.grid.minor = element_line(colour = NA,linetype = "blank"), 
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(colour = NA)) +
  labs(y = "% Relative cell source", fill = NULL)+labs(x = NULL)+
  scale_fill_manual(values=col)

ggsave("N_tc_spleen.png",width = 8,height = 6)





