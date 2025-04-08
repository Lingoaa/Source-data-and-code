rm(list = ls())
setwd('/home/rongge/QQ')

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

##加载数据
load("NTimmu_spleen.RData")

DefaultAssay(scRNA_spleen) <- "RNA"


## QC
scRNA_spleen[["percent.mt"]] <- PercentageFeatureSet(scRNA_spleen, pattern = "^mt-")

scRNA_spleen[['percent.ribo']] <- PercentageFeatureSet(scRNA_spleen, pattern = "^Rp[sl]")

Idents(scRNA_spleen) = scRNA_spleen$hash.ID

VlnPlot(scRNA_spleen, features = c("nFeature_RNA", "nCount_RNA", "percent.mt",'percent.ribo'),pt.size = 0,ncol = 2)

table(scRNA_spleen$hash.ID,scRNA_spleen$percent.ribo < 35)

table(scRNA_spleen$hash.ID,scRNA_spleen$percent.mt < 5)

plot1 <- FeatureScatter(scRNA_spleen,
                        feature1 = "nCount_RNA",
                        feature2 = "percent.mt")

plot2 <- FeatureScatter(scRNA_spleen,
                        feature1 = "nCount_RNA",
                        feature2 = "nFeature_RNA")

plot3 <- FeatureScatter(scRNA_spleen, 
                        feature1 = "nFeature_RNA", 
                        feature2 = "percent.ribo")

p = plot1 + plot2 + plot3

p

ggsave('plot_aom_spleen.png', p , width = 12,height = 9)

##QC
#rm(scRNA)
scRNA_spleen <- subset(scRNA_spleen,subset =  nCount_RNA < 30000 & nFeature_RNA>200 & nFeature_RNA < 6000 & percent.mt < 10 & percent.ribo < 40)

scRNA_spleen

#Normalizing the data
scRNA_spleen <- NormalizeData(scRNA_spleen, normalization.method = "LogNormalize", scale.factor = 10000)

#Identification of highly variable features (feature selection)
scRNA_spleen <- FindVariableFeatures(scRNA_spleen, selection.method = "vst", nfeatures = 2000)

#Scaling the data
all.genes <- rownames(scRNA_spleen)

scRNA_spleen <- ScaleData(scRNA_spleen, features = all.genes)

#Perform linear dimensional reduction
scRNA_spleen <- RunPCA(scRNA_spleen, features = VariableFeatures(object = scRNA_spleen))

p1 <- DimPlot(object = scRNA_spleen, reduction = "pca", pt.size = .5, group.by = "hash.ID")
p2 <- VlnPlot(object = scRNA_spleen, features = "PC_1", group.by = "hash.ID", pt.size = .3)
plot_grid(p1,p2)

##Run Harmony
scRNA_spleen <- scRNA_spleen %>%
  RunHarmony("hash.ID", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(scRNA_spleen, 'harmony')
harmony_embeddings[1:5, 1:5]

p1 <- DimPlot(object = scRNA_spleen, reduction = "harmony", pt.size = .3, group.by = "hash.ID")
p2 <- VlnPlot(object = scRNA_spleen, features = "harmony_1", group.by = "hash.ID", pt.size = .3)
plot_grid(p1,p2)

ggsave('降维_spleen.png' , p , width = 12,height = 9)

scRNA_spleen


##去双胞方法
vars.genes <- scRNA_spleen@assays[["RNA"]]@var.features

spleen.sce <- as.SingleCellExperiment(scRNA_spleen)


# 方法一：基于模拟方法
dbl.dens <- computeDoubletDensity(spleen.sce, subset.row=vars.genes, d=ncol(reducedDim(spleen.sce)))

spleen.sce$DoubletScore <- dbl.dens

dbl.calls <- doubletThresholding(data.frame(score=dbl.dens),method="griffiths", returnType="call", )

# 方法二：基于聚类结果
scRNA_spleen <- FindNeighbors(scRNA_spleen, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.1)
scRNA_spleen <- RunUMAP(scRNA_spleen, reduction = "harmony", dims = 1:20)
DimPlot(scRNA_spleen, reduction = "umap",label = T,pt.size = 0.5)

spleen.sce.dbl <- scDblFinder(spleen.sce, clusters=colData(spleen.sce)$seurat_clusters) 

layout_matrix <- rbind(c(1, 2),c(3, 4))

p = grid.arrange(p1, p2, p3, p4, layout_matrix = layout_matrix)

ggsave("去双胞_spleen.png", p, width=8, height=7)

summary(spleen.sce.dbl$scDblFinder.class == "doublet" & dbl.calls == "doublet")

scRNA_spleen <- scRNA_spleen[, !(spleen.sce.dbl$scDblFinder.class == "doublet" & dbl.calls == "doublet")]

scRNA_spleen 

save(scRNA_spleen, file = "alldata_spleen_KOdoublet.RData")

##Downstream analysis
scRNA_spleen <- FindNeighbors(scRNA_spleen, reduction = "harmony", dims = 1:20) %>% FindClusters(resolution = 0.1)

scRNA_spleen <- RunUMAP(scRNA_spleen, reduction = "harmony", dims = 1:20)

DimPlot(scRNA_spleen, reduction = "umap",label = T,pt.size = 0.1) 

library(viridis)


#细胞类型注释
new.cluster.ids <- c("0" = "B cells",
                                "1" = "T cells",
                                "2" = "other",
                                "3" = "Neutrophils",
                                "4" = "Monocytes",
                                "5" = "Macrophages",
                                "6" = "NK",
                                "7" = "DC",
                                "8" = "Mast cells",
                                "9" = "B cells")

names(new.cluster.ids) <- levels(scRNA_spleen)
scRNA_spleen <- RenameIdents(scRNA_spleen, new.cluster.ids)

table(scRNA_spleen@meta.data$orig.ident,scRNA_spleen@active.ident)

##去除非免疫细胞
imm_scRNA <- subset(x=scRNA_spleen, idents =c("Neutrophils","Monocytes","Macrophages","T cells","DC","B cells","Mast cells","NK"))

imm_scRNA

table(imm_scRNA@meta.data$hash.ID,imm_scRNA@active.ident)

rm(scRNA)

#查看细胞数
table(imm_scRNA@meta.data$hash.ID,imm_scRNA@active.ident)

ident <- data.frame(imm_scRNA@active.ident)

head(ident)

imm_scRNA@meta.data$cell_type <- ident$imm_scRNA.active.ident

table(imm_scRNA@meta.data$cell_type)

write.csv(imm_scRNA@meta.data, "allimmu_spleen_meta.csv")

saveRDS(imm_scRNA, file = "all_spleen_imm_seurat.rds")

imm_exp <- as.matrix(imm_scRNA@assays$RNA@counts)

saveRDS(imm_exp, file = "allimm_spleen_exp.rds")

rm(scRNA_spleen)

##ggplot画图
##绘图
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

# plot

color <- c("#a14462","#F5AE6B","#6a51a3","#4DBBD5B2","#b0d45d","#7E6148B2","#0072B2","#6a51a3","#9e9ac8")

plot1 <- DimPlot(imm_scRNA, reduction = 'umap', label = T,
                 pt.size = .1, cols = color) +
  NoAxes() + 
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab))
plot1

ggsave('spleen_Celltype_umap.png',plot1,width = 8,height = 8)

color <- c("#a14462","#F5AE6B","#b0d45d","#6a51a3","#9e9ac8")

plot2 <- DimPlot(imm_scRNA, reduction = 'umap', label = F,
                 pt.size = .1, cols = color, split.by  = "orig.ident",ncol = 2) +
  NoAxes() + 
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab))
plot2


##绘制 DotPlot
genes_to_check = c("Cd79b","Cd79a","Ebf1","Ly6d","Bank1",
                   "Cd3d","Cd3e","Cd3g",
                   "S100a8","S100a9",
                   "S100a4","Ms4a6c","Fn1","F13a1",
                   "C1qc","C1qa","Vcam1","C1qb",
                   "Ccl5","Klrd1","Nkg7",
                   "Itgax","Xcr1","Ccnd1","Naaa","Flt3",
                   "Mcpt8", "Cd200r3", "Gata2", "Cpa3")

p <- DotPlot(imm_scRNA, features = genes_to_check,assay='RNA',cols = c(low="white",high="darkred"),
             col.min = 0, dot.min = 0 )  + RotatedAxis() + coord_flip()
p
ggsave("marker_dotplot.png", p, width=8,height=8)


#FeaturePlot
library(viridis)
pal <- viridis(n = 15, option = "C", direction = -1)

#FeaturePlot
p1 = FeaturePlot(imm_scRNA,features = c("Cd79b","Cd79a",
                                        "Cd3d","Cd3e","Cd3g",
                                        "S100a8","S100a9",
                                        "S100a4","Ms4a6c",
                                        "C1qc","C1qa","Vcam1",
                                        "Xcr1","Naaa",
                                        "Ccl5","Klrd1","Nkg7",
                                        "Gata2"),cols = pal, order = T, ncol =6)
p1

ggsave("imm_marker_FeaturePlot.png", p1, width=16,height=7)

rm(scRNA_spleen)

##提取normal/tumor
table(imm_scRNA@meta.data$group)

#提取normal
N_spleen_imm <- subset(x = imm_scRNA, subset = group == "normal")
N_spleen_imm_exp <- as.matrix(N_spleen_imm@assays$RNA@counts)
saveRDS(N_spleen_imm_exp, file = "N_spleen_imm_exp.rds")
write.csv(N_spleen_imm@meta.data, "N_spleen_imm_meta.csv")
#saveRDS(N_spleen_imm, file = "N_spleen_imm_seurat.rds")

rm(N_spleen_imm_exp)

plot1 <- DimPlot(T_spleen_imm, reduction = 'umap', label = T,
                 pt.size = .1, cols = color) +
  NoAxes() + 
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab))
plot1

#提取tumor
T_spleen_imm <- subset(x = imm_scRNA, subset = group == "tumor")
T_spleen_imm_exp <- as.matrix(T_spleen_imm@assays$RNA@counts)
saveRDS(T_spleen_imm_exp, file = "T_spleen_imm_exp.rds")
write.csv(T_spleen_imm@meta.data, "T_spleen_imm_meta.csv")
#saveRDS(T_spleen_imm, file = "T_spleen_imm_seurat.rds")

rm(T_spleen_imm_exp)

##提取B/T cell
BC <- subset(x=imm_scRNA, cell_type == "B cells")
TC <- subset(x=imm_scRNA, cell_type == "T cells")

BC_exp <- as.matrix(BC@assays$RNA@counts)
TC_exp <- as.matrix(TC@assays$RNA@counts)

saveRDS(BC_exp, file = "NT_spleen_BC_exp.rds")
saveRDS(TC_exp, file = "NT_spleen_TC_exp.rds")

write.csv(BC@meta.data, "NT_spleen_BC_meta.csv")
write.csv(TC@meta.data, "NT_spleen_TC_meta.csv")