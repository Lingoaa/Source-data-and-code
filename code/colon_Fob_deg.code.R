setwd("E:/NT_spleen_bc/N_slpeen_bc")

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
  library(ggpubr)
  library(Seurat)
  library(stringr)
  
}

##导入数据
N_bc <- readRDS("N_bc_seurat.rds")

N_bc

table(N_bc$cell_type)


############分析B细胞大类的基因差异
##提取Pro-B cell
PC <- subset(x=N_bc, cell_type == "Plasma cell")

table(PC@meta.data$orig.ident)

PC$patientoutcone =  paste(PC$hash.ID, PC$orig.ident,sep = '_')


# 通过Seurat的AverageExpression函数计算每个样本基因的平均表达量
avg <- AverageExpression(object = PC, assay = "RNA", group.by = "patientoutcone")
a = avg$RNA

#想要呈现的基因

C1 = c('Dhx9',"Ighg1","Ighg3","Igha","Jchain")

#可视化
gene_boxplot <- function(gene){
  b = as.data.frame(a[gene,])
  colnames(b) = gene
  phe = b
  group = ifelse(str_detect(rownames(phe),"N-ko"),"N-ko","N-wt")
  phe$group = group
  
  phe$group = factor(phe$group,levels =c("N-ko","N-wt"))
  table(phe$group)
  colnames(phe)
  
  col = c("#DF5E38","#7B7DC4")
  P_T <- ggplot(phe,aes(x=group,y=!!rlang::sym(gene),fill=group)) +
    geom_boxplot(size=0.5,fill="white",outlier.fill="white",outlier.color="white",color=col) +
    geom_jitter(aes(fill=group),width =0.2,shape = 21,size=2) +
    scale_fill_manual(values = col)+
    ggtitle(gene) +
    theme_bw() +
    theme(panel.grid = element_blank(), panel.background = element_blank(),
          axis.line = element_line(color = 'black')) +
    ylab("Average signature score") +xlab(" ")+
    theme(legend.position = "none",
          axis.text.x = element_text(colour = "black", family = "Times", size = 14),
          axis.text.y = element_text(family = "Times", size = 14, face = "plain"),
          axis.title.y = element_text(family = "Times", size = 14, face = "plain"),
          axis.title.x = element_text(family = "Times", size = 14, face = "plain"),
          plot.title = element_text(family = "Times", size = 17, face = "bold", hjust = 0.5),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank())
  
  P_T
  compaired = list(c("N-ko", "N-wt"))
  P_T = P_T+geom_signif(comparisons = compaired,step_increase = 0.1,map_signif_level = T,test = t.test,tip_length = 0)
  P_T
  
}

P_C1 = gene_boxplot(C1[1])+gene_boxplot(C1[2])+gene_boxplot(C1[3])+gene_boxplot(C1[4])+
  gene_boxplot(C1[5])+gene_boxplot(C1[6])+gene_boxplot(C1[7])+gene_boxplot(C1[8])+gene_boxplot(C1[9])

P_C2 = gene_boxplot(C2[1])+gene_boxplot(C2[2])+gene_boxplot(C2[3])+gene_boxplot(C2[4])+
  gene_boxplot(C2[5])+gene_boxplot(C2[6])+gene_boxplot(C2[7])+gene_boxplot(C2[8])

P_C1/P_C2
p = P_C1/P_C2

ggsave('pro_gene_outcome.png', p,width = 10,height = 14)
ggsave('pro_gene_boxplot_C1.png', P_C1, width = 10,height = 12)
ggsave('pro_gene_boxplot_C2.png', P_C2, width = 10,height = 12)



##############差异基因分析
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

##先分析大类B细胞
table(N_bc@meta.data$cell_type)

Mem <- subset(x=bc, cell_type == "Memory B cells")

##DEGs
#设置min.pct = 0.5参数过滤掉那些在50%以下细胞中检测到的基因 
DefaultAssay(N_bc) <- "RNA"

#方法一
deg_bc <- FindMarkers(N_bc, min.pct = 0.1, 
                      logfc.threshold = 0.25,
                      group.by = "orig.ident",
                      ident.1 ="N-wt",
                      ident.2= "N-wt",
                      only.pos = T) 

write.csv(deg_bc,"bc_DEG.csv")


##方法二：计算差异基因
dge.sample <- FindMarkers(N_bc, ident.1 = 'N-wt',ident.2 = 'N-wt',
                          group.by = 'orig.ident',logfc.threshold = 0,min.pct = 0)

sig_dge.all <- subset(dge.sample, p_val_adj<0.05&abs(avg_log2FC)>0.25) #所有差异基因
View(sig_dge.all)

write.csv(sig_dge.all,'bc_dge_0.25.csv')



#火山图可视化
library(EnhancedVolcano)

EnhancedVolcano(deg_bc,
                lab = rownames(deg_bc),
                x = 'avg_log2FC',
                y = 'p_val_adj',
                pCutoff = 0.05,
                FCcutoff = 0.25,
                pointSize = 3.0,
                labSize = 6.0,
                title = 'aom_deg_bc WT vs KO')
##ggplot
sig_dge.all$threshold<-as.factor(ifelse(deg_bc$log2FoldChange <= -0.25,'Up',ifelse(res$padj<0.05 & res$log2FoldChange >= 0.25,'Down','Not')))
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
