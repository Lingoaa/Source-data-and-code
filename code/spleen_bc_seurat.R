setwd("E:/spleen_alldata/NT_spleen_BC")

rm(list=ls())
set.seed(123)  

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

##导入数据
exp <- readRDS("NT_spleen_BC_exp.rds")

meta <- read.csv("NT_spleen_BC_meta.csv", row.names = 1)


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
hash.ID <- meta$hash.ID
names(hash.ID) <- colnames(x = bc)
bc <- AddMetaData(
  object = bc,
  metadata = hash.ID,
  col.name = "hash.ID"
)
table(bc@meta.data$hash.ID)

#group信息
group <- meta$group
names(group) <- colnames(x = bc)
bc <- AddMetaData(
  object = bc,
  metadata = group,
  col.name = "group"
)
table(bc@meta.data$group)

save(bc, file = "SBC.RData")
load("SBC.RData")

table(bc@meta.data$orig.ident)


p1 <- DimPlot(object = bc, reduction = "pca", pt.size = .1,group.by = "hash.ID")
p2 <- VlnPlot(object = bc, features = "PC_1", pt.size = .1,group.by = "hash.ID")
plot_grid(p1,p2)
p = plot_grid(p1,p2)
ggsave("sbc_plot.png", p ,width=9 ,height=6)


##Run Harmony
bc <- bc %>%
  RunHarmony("hash.ID", plot_convergence = TRUE)

harmony_embeddings <- Embeddings(bc, 'harmony')
harmony_embeddings[1:5, 1:5]

p1 <- DimPlot(object = bc, reduction = "harmony", pt.size = .1,group.by = "hash.ID")
p2 <- VlnPlot(object = bc, features = "harmony_1", pt.size = .1,group.by = "hash.ID")
plot_grid(p1,p2)



##Downstream analysis
bc <- FindNeighbors(bc, reduction = "harmony", dims = 1:15) %>% FindClusters(resolution = 0.3)

bc <- RunUMAP(bc, reduction = "harmony", dims = 1:15)

DimPlot(bc, reduction = "umap",label = T,pt.size = 0.5) 

##找高变基因
markers <- FindAllMarkers(object = bc, test.use="wilcox" ,
                          only.pos = TRUE,
                          logfc.threshold = 0.25)   

all.markers =markers %>% dplyr::select(gene, everything()) %>% subset(p_val<0.05)

top30 = all.markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)

write.csv(top30, "bc_markers.csv")

new.cluster.ids <- c( "0" = "Follicular BC",
                      "1" = "MZ BC",
                      "2" = "GC BC",
                      "3" = "Pro/Pre BC",
                      "4" = "Plasmblast",
                      "5" = "IFI+ BC",
                      "6" = "Transitional BC",
                      "7" = "Plasma cell",
                      "8" = "Follicular BC")

names(new.cluster.ids) <- levels(bc)
bc <- RenameIdents(bc, new.cluster.ids)

#查看细胞数
table(bc@meta.data$orig.ident,bc@active.ident)
table(bc@meta.data$hash.ID,bc@active.ident)

ident <- data.frame(bc@active.ident)

head(ident)

bc@meta.data$cell_type <- ident$bc.active.ident

table(bc@meta.data$cell_type)
table(bc@meta.data$hash.ID, bc@meta.data$cell_type)

write.csv(bc@meta.data, "NT_sbc_meta.csv")

saveRDS(bc, file = "NT_spleen_bc_seurat.rds")

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
color <- c("#f39c90", "#4DBBD5B2","#A9D179","#00A087B2","#F5AE6B","#CC79A7","#4387B5","#7E6148B2")
color <- c("#f39c90","#4DBBD5B2", "#A9D179", '#4f9568', "#F5AE6B", '#c57aa9', '#4387B5', '#7E6148B2')

plot1 <- DimPlot(bc, reduction = 'umap', label = F,repel = TRUE,
                 pt.size = .5, cols = color,label.size = 3) +
  NoAxes() +
  theme(aspect.ratio = 1) +
  geom_line(data = axes,
            aes(x = x,y = y,group = group),
            arrow = arrow(length = unit(0.1, "inches"),
                          ends="last", type="closed")) +
  geom_text(data = label,
            aes(x = x,y = y,angle = angle,label = lab))
plot1



##热图marker可视化
##DotPlot
genes_to_check = c("Ighd","Ccr7","Cd79a",
                   "Cr2","Dtx1","Myc","Pik3r4","Tm6sf1",
                   "Mif","Nme1","Nme2","Eif5a","Ncl",
                   "Vpreb3","Cd24a","Sox4",
                   "Ctla4","Cd44","Id2","Sox5","Actn1","Pstpip2",
                   "Ifi27l2a","Ifit3","Ifi213","Isg15",
                   "Fcer2a","Zfp318","Pgap1",
                   "Jchain","Iglc1","Ighm","Igkc","Mzb1")

p <- DotPlot(bc, features = genes_to_check,assay='RNA',cols = c(low="white",high="darkred"),
             col.min = 0, dot.min = 0)  + RotatedAxis() + coord_flip()
p
ggsave("marker_dotplot.png", p, width=9,height=7)

#FeaturePlot
library(viridis)
pal <- viridis(n = 15, option = "C", direction = -1)

p2 = FeaturePlot(bc,features = c("Ighd","Ccr7","Cd79a",
                                       "Cr2","Dtx1",
                                       "Mif","Nme1",
                                       "Vpreb3","Sox4",
                                       "Ctla4","Cd44",
                                       "Ifi27l2a","Ifit3",
                                       "Fcer2a","Zfp318",
                                       "Jchain","Iglc1"),cols = pal, order = T, ncol =5)


ggsave("marker_FeaturePlot.png", p2, width=12,height=10)


##提取normal/tumor
table(bc@meta.data$group)

#提取normal
N_bc <- subset(x = bc, subset = group == "normal")
saveRDS(N_bc, file = "N_bc_seurat.rds")


#提取tumor
T_bc <- subset(x = bc, subset = group == "tumor")
saveRDS(T_bc, file = "T_bc_seurat.rds")

rm(T_bc_exp)

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
rm(phe)
phe=T_bc@meta.data
colnames(phe)
tb=table(phe$cell_type,
         phe$hash.ID)

head(tb)

#统计细胞数量
balloonplot(tb, main ="Tumor spleen B cells", xlab ="celltype", ylab="sample",
            label = T, show.margins = T)


bar_data <- as.data.frame(tb)

bar_per <- bar_data %>% 
  group_by(Var2) %>%
  mutate(sum(Freq)) %>%
  mutate(percent = Freq / `sum(Freq)`)
head(bar_per) 
write.csv(bar_per,file = "celltype_by_group_percent.csv")
col =c("#2873B3","#2EBEBE","#8264CC","#74B346","#167153","#F1CC2F","#e8743c","#7D4444","#A14462")

col = c("#E64B35B2", "#2EBEBE","#8264CC","#74B346","#167153","#F1CC2F","#7D4444","#7E6148B2",
        "#8491B4B2", "#91D1C2B2", "#7E6148B2")
ggplot(bar_per, aes(x = percent, y = Var2)) +
  geom_bar(aes(fill = Var1) , stat = "identity") + coord_flip() +
  theme(axis.ticks = element_line(linetype = "blank"),
        legend.position = "right",
        panel.grid.minor = element_line(colour = NA,linetype = "blank"), 
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(colour = NA)) +
  labs(y = "% Relative cell source", fill = NULL)+labs(x = NULL)+
  scale_fill_manual(values=col)

ggsave("NTbc_spleen.png",width = 8,height = 4)


###不同分组中细胞亚群的比例
##分为ko和wt组

table(phe$cell_type)
table(phe$hash.ID)
phe$outcone =  paste(phe$hash.ID,phe$orig.ident,sep = '_')


phe_immune = phe[phe$cell_type %in% c('Follicular BC','MZ BC','GC BC','Plasmblast',
                                      'Pro/Pre BC','IFI+ BC','Transitional BC','Plasma cell'),]

tb=table(phe_immune$cell_type,
         phe_immune$outcone)

head(tb)
balloonplot(tb)
bar_data <- as.data.frame(tb)


bar_per <- bar_data %>% 
  group_by(Var2) %>%
  mutate(sum(Freq)) %>%
  mutate(percent = Freq / `sum(Freq)`)
head(bar_per) 
bar_per$Var2
group = ifelse(str_detect(bar_per$Var2,"T-ko"),"T-ko","T-wt")

bar_per$group = group

bar_per$group = factor(bar_per$group,levels =c("T-ko","T-wt"))
table(bar_per$Var1)
bar_per$Var1 = factor(bar_per$Var1,levels =c('Follicular BC','MZ BC','GC BC','Plasmblast',
                                             'Pro/Pre BC','IFI+ BC','Transitional BC','Plasma cell'))

###拼接7个箱线图
library(gridExtra)
library(ggpubr)

p <- list()

col = c("#BD6263","#8EA325")
for(i in 1:8){
  # 选择当前细胞类型的数据
  cell <- bar_per[bar_per$Var1 == unique(bar_per$Var1)[i],]
  
  p[[i]] <- ggplot(cell, aes(x = group, y = percent, fill = group)) +
    geom_boxplot(size = 0.5, fill = "white", outlier.fill = "white", outlier.color = "white",color = col) +
    geom_jitter(aes(fill = group), width = 0.2, shape = 21, size = 2) +
    scale_fill_manual(values = col) +
    ggtitle(paste0(unique(bar_per$Var1)[i], " ")) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(colour = "black", family = "Times", size = 14),
          axis.text.y = element_text(family = "Times", size = 14, face = "plain"),
          axis.title.y = element_text(family = "Times", size = 14, face = "plain"),
          axis.title.x = element_text(family = "Times", size = 14, face = "plain"),
          plot.title = element_text(family = "Times", size = 15, face = "bold", hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    ylab("Fraction (%) of immune cells") + xlab(" ")
  # 移除y轴标题和刻度标签
  compaired = list(c("T-ko", "T-wt"))
  p[[i]] <- p[[i]] + theme(axis.title.y = element_blank())+
    geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = T,test = t.test,tip_length = 0)
}

combined_plot <- grid.arrange(grobs = p, ncol = 4, left = "Fraction (%) of cells") # 设置共同的纵轴标题和拼图列数

print(combined_plot)

ggsave("Tumor_boxplot_celltype.png", combined_plot, width = 13,height = 9)







