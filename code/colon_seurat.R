setwd("E:/AOMDSS/filter_data")
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

load("5AOMDSS.RData")

scRNA

## QC
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^mt-")

scRNA[['percent.ribo']] <- PercentageFeatureSet(scRNA, pattern = "^Rp[sl]")
head(scRNA$percent.ribo)

Idents(scRNA) = scRNA$orig.ident

VlnPlot(scRNA, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.ribo'), pt.size = 0, ncol = 2)

table(scRNA$orig.ident,scRNA$percent.ribo < 25)


plot1 <- FeatureScatter(scRNA,
                        feature1 = "nCount_RNA",
                        feature2 = "percent.mt")

plot2 <- FeatureScatter(scRNA,
                        feature1 = "nCount_RNA",
                        feature2 = "nFeature_RNA")

plot3 <- FeatureScatter(scRNA, 
                        feature1 = "nFeature_RNA", 
                        feature2 = "percent.ribo")

p = plot1 + plot2 + plot3

p

##QC
scRNA <- subset(scRNA, subset = nFeature_RNA > 200 & nFeature_RNA < 4000 & percent.mt < 7 & percent.ribo < 25)

scRNA


#Normalizing the data
scRNA <- NormalizeData(scRNA, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of highly variable features (feature selection)
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 3000)

#Scaling the data
all.genes <- rownames(scRNA)
scRNA <- ScaleData(scRNA, features = all.genes)

#Perform linear dimensional reduction
scRNA <- RunPCA(scRNA, features = VariableFeatures(object = scRNA))

p1 <- DimPlot(object = scRNA, reduction = "pca", pt.size = .1, group.by = "orig.ident")
p2 <- VlnPlot(object = scRNA, features = "PC_1", group.by = "orig.ident", pt.size = .1)
plot_grid(p1,p2)


##Run Harmony
scRNA <- scRNA %>%
  RunHarmony("orig.ident", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(scRNA, 'harmony')
harmony_embeddings[1:5, 1:5]

p1 <- DimPlot(object = scRNA, reduction = "harmony", pt.size = .3, group.by = "orig.ident")
p2 <- VlnPlot(object = scRNA, features = "harmony_1", group.by = "orig.ident", pt.size = .3)
p = plot_grid(p1,p2)
p

ggsave('降维.png' , p , width = 12,height = 9)

scRNA

scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.3)

scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:20)

DimPlot(scRNA, reduction = "umap",label = T,pt.size = 0.5) 

#scRNA <- RunTSNE(scRNA, reduction = "harmony", dims = 1:20)


##去双胞方法
vars.genes <- scRNA@assays[["RNA"]]@var.features

spleen.sce <- as.SingleCellExperiment(scRNA)

dbl.dens <- computeDoubletDensity(spleen.sce, subset.row=vars.genes, d=ncol(reducedDim(spleen.sce)))

spleen.sce$DoubletScore <- dbl.dens

dbl.calls <- doubletThresholding(data.frame(score=dbl.dens),method="griffiths", returnType="call", )

scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.3)
scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:20)
DimPlot(scRNA, reduction = "umap",label = T,pt.size = 0.5)

spleen.sce.dbl <- scDblFinder(spleen.sce, clusters=colData(spleen.sce)$seurat_clusters) 

##结果可视化
p1 <- plotUMAP(spleen.sce, colour_by = "DoubletScore")

p2 <- plotUMAP(spleen.sce.dbl, colour_by="scDblFinder.score")

p3 <- plotColData(spleen.sce, x="seurat_clusters", y="DoubletScore", colour_by=I(dbl.calls))

p4 <- plotColData(spleen.sce.dbl, x="seurat_clusters", y="scDblFinder.score", colour_by=I(spleen.sce.dbl$scDblFinder.class))

layout_matrix <- rbind(c(1, 2),c(3, 4))

p = grid.arrange(p1, p2, p3, p4, layout_matrix = layout_matrix)

ggsave("去双胞_v2.png", p, width=8, height=7)

summary(spleen.sce.dbl$scDblFinder.class == "doublet" & dbl.calls == "doublet")

scRNA <- scRNA[, !(spleen.sce.dbl$scDblFinder.class == "doublet" & dbl.calls == "doublet")]

scRNA

table(scRNA@meta.data$orig.ident)

saveRDS(scRNA, file = "colon_KOdoublet.rds")

load("colon_KOdoublet.rds")

rm(spleen.sce,spleen.sce.dbl)

##Downstream analysis
scRNA <- FindNeighbors(scRNA, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.2)

scRNA <- RunUMAP(scRNA, reduction = "harmony", dims = 1:20)

##找高变基因
markers <- FindAllMarkers(object = scRNA, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)   

all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)

top20 = all.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)

write.csv(top20, "5colon_markers_0.3.csv")

genes_to_check = c("S100a9","S100a8","Retnlg","G0s2","Hdc",
                   "Plac8","Inhba","Plcb1","Fn1",
                   "Frmd4b","Mrc1","C1qc","C1qb","C1qa","Adgre1",
                   "Cd3g","Cd3e","Cd3d","Camk4","Skap1","Rora",
                   "Cst3","Crip1","Tbc1d8","Flt3",
                   "Cd79a","Cd79b",
                   "Epcam","Krt8","Krt19","Lgals4",
                   "Dcn","Col3a1","Mgp",
                   "Cpa3","Cd200r3", "Gata2")

p <- DotPlot(scRNA, features = genes_to_check,assay='RNA' )  + RotatedAxis() + coord_flip()
p
ggsave("5colon_dotplot.png", p, width=8,height=8)

#细胞类型注释
new.cluster.ids <- c("0" = "Neutrophils",
                     "1" = "Monocytes",
                     "2" = "Neutrophils",
                     "3" = "Macrophages",
                     "4" = "T cells",
                     "5" = "DC",
                     "6" = "B cells",
                     "7" = "Fibroblasts",
                     "8" = "Epithelial cells",
                     "9" = "B cells",
                     "10" = "Mast cells")


names(new.cluster.ids) <- levels(scRNA)
scRNA <- RenameIdents(scRNA, new.cluster.ids)

table(scRNA@meta.data$orig.ident,scRNA@active.ident)


##去除非免疫细胞
imm_scRNA <- subset(x=scRNA, idents =c("Neutrophils","Monocytes","Macrophages","T cells","DC","B cells","Mast cells"))

imm_scRNA

table(imm_scRNA@meta.data$orig.ident)

table(imm_scRNA@meta.data$orig.ident,imm_scRNA@active.ident)

rm(scRNA)

#查看细胞数
table(imm_scRNA@meta.data$orig.ident,imm_scRNA@active.ident)

ident <- data.frame(imm_scRNA@active.ident)

head(ident)

imm_scRNA@meta.data$cell_type <- ident$imm_scRNA.active.ident

table(imm_scRNA@meta.data$cell_type)

write.csv(imm_scRNA@meta.data, "5colon_imm_meta.csv")



##ggplot绘图
##提取前两个主成分数据
# extact PC ranges
pc12 <- Embeddings(object = imm_scRNA,reduction = 'umap') %>%
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
color <- c("#eb998b", "#4DBBD5B2","#A9D179","#84CAC0","#F5AE6B","#6a51a3","#a14462")

plot1 <- DimPlot(imm_scRNA, reduction = 'umap', label = T,repel = TRUE,
                 pt.size = .2, cols = color,label.size = 4) +
  NoAxes() + NoLegend() +
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab))

plot1

plot2 <- DimPlot(imm_scRNA, reduction = 'umap', label = T,repel = TRUE,
                 pt.size = .2, cols = color,label.size = 4,group.by = 'group') +
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


plot3 <- DimPlot(imm_scRNA, reduction = 'umap', label = F,repel = TRUE,
                 pt.size = .2, cols = color,label.size = 4,split.by  = "orig.ident") +
  NoAxes() +
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab))

plot3

ggsave('imm_ident.png',plot3,width = 12,height = 8)


##所有细胞marker
genes_to_check = c("S100a9","S100a8","Retnlg","G0s2","Hdc",
                   "Plac8","Inhba","Plcb1","Fn1",
                   "Frmd4b","Mrc1","C1qc","C1qb","C1qa","Adgre1",
                   "Cd3g","Cd3e","Cd3d","Camk4","Skap1","Rora",
                   "Cst3","Crip1","Tbc1d8","Flt3",
                   "Cd79a","Cd79b",
                   "Cpa3","Cd200r3","Gata2")

p <- DotPlot(imm_scRNA, features = genes_to_check,assay='RNA',cols = c(low="white",high="darkred"),
             col.min = 0, dot.min = 0 )  + RotatedAxis() + coord_flip()
p
ggsave("5imm_marker_dotplot.png", p, width=10,height=10)

##小提琴图
genes_to_check = c("G0s2",
                   "Inhba","Plcb1",
                   "Cd3g","Cd3e","Camk4",
                   "Mrc1","Adgre1",
                   "Flt3",
                   "Cd79a","Cd79b",
                   "Cpa3","Gata2")

markers <- CaseMatch(genes_to_check, rownames(imm_scRNA))
markers <- as.character(markers)
VlnPlot(imm_scRNA, features = markers, pt.size = 0, group.by = 'cell_type',stack = T)+NoLegend()

p = VlnPlot(imm_scRNA, features = markers, pt.size = 0, group.by = 'cell_type', stack = T)+NoLegend()


ggsave("5imm_marker_VlnPlot.png", p, width=15,height=7)


#FeaturePlot
library(viridis)
pal <- viridis(n = 15, option = "C", direction = -1)

p1 = FeaturePlot(imm_scRNA,features = c("S100a8","G0s2","Hdc",
                               "Plac8","Plcb1","Fn1",
                               "Cd3g","Cd3e","Cd3d","Camk4",
                               "Frmd4b","Mrc1","Adgre1",
                               "Flt3",
                               "Cd79a","Cd79b",
                               "Cpa3","Gata2"),cols = pal, order = T, ncol =6)




ggsave("5imm_marker_FeaturePlot.png", p1, width=16,height=8)


saveRDS(imm_scRNA, file = "colon_imm_seurat.rds")

scRNA <- readRDS("colon_imm_seurat.rds")


##提取B/T cell
table(imm_scRNA@meta.data$cell_type)

#提取B cell
BC5 <- subset(x=imm_scRNA, idents =c("B cells"))
BC5_exp <- as.matrix(BC5@assays$RNA@counts)
saveRDS(BC5_exp, file = "BC5_exp.rds")
write.csv(BC5@meta.data, "BC5_meta.csv")
saveRDS(BC5, file = "BC5_seurat.rds")

rm(ident,label,all.genes)
rm(BC5,BC5_exp)

#提取T cell
TC5 <- subset(x=imm_scRNA, cell_type == "T cells")
TC5_exp <- as.matrix(TC5@assays$RNA@counts)
saveRDS(TC5_exp, file = "TC5_exp.rds")
write.csv(TC5@meta.data, "TC5_meta.csv")
saveRDS(TC5, file = "TC5_seurat.rds")

rm(ident,label,all.genes)
rm(TC5,TC5_exp)

saveRDS(imm_scRNA, file = "5colon_imm_seurat.rds")

