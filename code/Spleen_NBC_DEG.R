setwd("E:/spleen_alldata/NT_spleen_BC/normal_SBC")

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

##导入数据
bc <- readRDS("N_bc_seurat.rds")

table(bc@meta.data$cell_type)

Fob <- subset(x=bc, cell_type == "Follicular BC")

##DEGs
#设置min.pct = 0.5参数过滤掉那些在50%以下细胞中检测到的基因 
DefaultAssay(Fob) <- "RNA"

deg_bc <- FindMarkers(Fob, min.pct = 0.1, 
                      logfc.threshold = 0.25,
                      group.by = "orig.ident",
                      ident.1 ="N-ko",
                      ident.2= "N-wt",
                      only.pos = T) 

write.csv(deg_bc,"Normal_Fob_DEG.csv")


##计算差异基因
dge.sample <- FindMarkers(Fob, ident.1 = 'N-ko',ident.2 = 'N-wt',
                          group.by = 'orig.ident',logfc.threshold = 0,min.pct = 0)

sig_dge.all <- subset(dge.sample, p_val_adj<0.05&abs(avg_log2FC)>0.25) #所有差异基因
View(sig_dge.all)

write.csv(sig_dge.all,'Normal_Fob_allDEG.csv')



#火山图可视化
library(EnhancedVolcano)

EnhancedVolcano(sig_dge.all,
                lab = rownames(sig_dge.all),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.25,
                pointSize = 3.0,
                labSize = 6.0,
                title = 'Fob_DEG KO/WT')
##ggplot
sig_dge.all$threshold<-as.factor(ifelse(sig_dge.all$log2FoldChange <= -0.25,'Down',ifelse(res$padj<0.05 & res$log2FoldChange >= 0.25,'Up','Not')))
sig_dge.all<-data.frame(sig_dge.all)
p3<-ggplot(data=sig_dge.all, aes(x=log2FoldChange, y=-log10(padj), colour=threshold, fill=threshold)) + 
  scale_color_manual(values=c("red", "grey","blue"))+
  geom_point(alpha=0.6,size=2) +
  #xlim(c(-6, 6)) +
  #ylim(c(0, 300)) +
  theme_bw(base_size = 12, base_family = "Times") +
  geom_vline(xintercept=c(-2,2),lty=1,col=c("blue","red"),lwd=0.6)+
  geom_hline(yintercept = -log10(0.05),lty=2,col="black",lwd=0.6)+
  theme(legend.position="right",
        panel.grid=element_blank(),
        legend.title = element_blank(),
        legend.text= element_text(face="bold", color="black",family = "Times", size=8),
        plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(color="black", size=12),
        axis.text.y = element_text(color="black", size=12),
        axis.title.x = element_text(face="bold", color="black", size=12),
        axis.title.y = element_text(face="bold",color="black", size=12))+
  labs(x="log2 (Fold Change)",y="-log10 (p-value)",title="bc_deg WT vs KO")