setwd("E:/AOM_tumor/filter_data/B_cell")
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
  library(viridis)
  
}

set.seed(123)  

##导入数据
exp <- readRDS("BC5_exp.rds")

meta <- read.csv("BC5_meta.csv", row.names = 1)


##初始化Seurat对象
bc <- CreateSeuratObject(counts =exp , project = "b_cell", min.cells = 3) %>%
  Seurat::NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(selection.method = "vst", nfeatures = 2000) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(pc.genes = bc@var.genes, npcs = 20, verbose = FALSE)


##添加分组信息
#分组信息
orig.ident <- meta$orig.ident
names(orig.ident) <- colnames(x = bc)
bc <- AddMetaData(
  object = bc,
  metadata = orig.ident,
  col.name = "orig.ident"
)
table(bc@meta.data$orig.ident)

#样本信息
group <- meta$group
names(group) <- colnames(x = bc)
bc <- AddMetaData(
  object = bc,
  metadata = group,
  col.name = "group"
)
table(bc@meta.data$group)

save(bc, file = "5BC.RData")

load("5BC.RData")

table(bc@meta.data$orig.ident)


p1 <- DimPlot(object = bc, reduction = "pca", pt.size = .1,group.by = "orig.ident")
p2 <- VlnPlot(object = bc, features = "PC_1", pt.size = .1,group.by = "orig.ident")
plot_grid(p1,p2)
p = plot_grid(p1,p2)
ggsave("bc_plot.png", p ,width=9 ,height=6)


##Run Harmony
bc <- bc %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(bc, 'harmony')
harmony_embeddings[1:5, 1:5]

p1 <- DimPlot(object = bc, reduction = "harmony", pt.size = .1,group.by = "orig.ident")
p2 <- VlnPlot(object = bc, features = "harmony_1", pt.size = .1,group.by = "orig.ident")
plot_grid(p1,p2)


##Downstream analysis
bc <- FindNeighbors(bc, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.2)

bc <- RunUMAP(bc, reduction = "harmony", dims = 1:20)

DimPlot(bc, reduction = "umap",label = T,pt.size = 0.5) 

##找高变基因
markers <- FindAllMarkers(object = bc, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)   

all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)

top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

write.csv(top20, "5bc_markers_0.2.csv")

table(bc@meta.data$orig.ident,bc@meta.data$RNA_snn_res.0.3)

##DotPlot
genes_to_check = c("Ccr7","Cd79a","Ebf1","Bank1",
                   "Jchain","Igha","Igkc","Iglv1",
                   "Sox5","Top2a","Pdzd2","Cecr2","Nt5e","Basp1")

p <- DotPlot(bc, features = genes_to_check,assay='RNA',cols = c(low="white",high="darkred"),
             col.min = 0, dot.min = 0)  + RotatedAxis() + coord_flip()
p
ggsave("marker_dotplot.png", p, width=7,height=7)


new.cluster.ids <- c("0" = "Follicular BC",
                     "1" = "Plasm cells",
                     "2" = "Follicular BC",
                     "3" = "Follicular BC",
                     "4" = "GC BC",
                     "5" = "Plasm cells",
                     "6" = "Follicular BC")

names(new.cluster.ids) <- levels(bc)
bc <- RenameIdents(bc, new.cluster.ids)

#查看细胞数
table(bc@meta.data$orig.ident,bc@active.ident)

ident <- data.frame(bc@active.ident)

head(ident)

bc@meta.data$cell_type <- ident$bc.active.ident

table(bc@meta.data$cell_type)
table(bc@meta.data$orig.ident, bc@meta.data$cell_type)

write.csv(bc@meta.data, "bc_meta.csv")

saveRDS(bc, file = "5bc_seurat.rds")


##导入数据
bc <- readRDS("bc_seurat_harmony.rds")

##ggplot绘图
##提取前两个主成分数据
# extact PC ranges
pc12 <- Embeddings(object = bc,reduction = 'umap') %>%
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
color <- c("#f39c90", "#4DBBD5B2","#A9D179","#00A087B2","#CC79A7")
color <- c("#f39c90","#A9D179","#CC79A7")

plot1 <- DimPlot(bc, reduction = 'umap', label = T,repel = TRUE,
                 pt.size = .8, cols = color, label.size = 4) +
  NoAxes() + NoLegend() +
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab))
plot1

plot2 <- DimPlot(bc, reduction = 'umap', label = F,repel = TRUE,
                 pt.size = .8, cols = color,split.by  = "group", label.size = 4,ncol = 2) +
  NoAxes() + 
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab))
plot2
ggsave('bc_celltype_group.png',plot1,width = 10,height = 10)


##热图marker可视化
##DotPlot
genes_to_check = c("Ccr7","Cd79a","Ebf1","Bank1",
                   "Jchain","Igha","Igkc","Iglv1",
                   "Sox5","Top2a","Pdzd2","Cecr2","Nt5e","Basp1")

p <- DotPlot(bc, features = genes_to_check,assay='RNA',cols = c(low="white",high="darkred"),
             col.min = 0, dot.min = 0)  + RotatedAxis() + coord_flip()
p
ggsave("marker_dotplot.png", p, width=7,height=7)


##小提琴图

genes_to_check = c("Ccr7","Cd79a","Ebf1","Bank1",
                   "Jchain","Igha","Igkc","Iglv1",
                   "Sox5","Top2a","Pdzd2","Cecr2","Nt5e","Basp1")
markers <- CaseMatch(genes_to_check, rownames(bc))
markers <- as.character(markers)
VlnPlot(bc, features = markers, pt.size = 0, group.by = 'cell_type', stack = T)+NoLegend()


pal <- viridis(n = 15, option = "C", direction = -1)

p1 = FeaturePlot(bc,features = c("Ccr7","Cd79a","Ebf1","Bank1",
                             "Jchain","Igha","Igkc","Iglv1",
                             "Sox5","Top2a","Pdzd2","Cecr2","Nt5e","Basp1"),cols = pal, order = T, ncol =4)

ggsave("marker_FeaturePlot.png", p1, width=12,height=8)



##提取Fob
table(bc@meta.data$cell_type)

Fob <- subset(x=bc, idents =c("Follicular BC"))
Fob_exp <- as.matrix(Fob@assays$RNA@counts)
saveRDS(Fob_exp, file = "Fob_exp.rds")
write.csv(Fob@meta.data, "Fob_meta.csv")
saveRDS(Fob, file = "Fob_seurat.rds")

