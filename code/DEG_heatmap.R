setwd("E:/spleen_alldata/NT_spleen_BC/normal_SBC/GC")

rm(list=ls())

##加载R包

{
  library(Seurat)
  library(tidyverse)
  library(reshape2)
  library(psych)
  library(ggpubr)
  library(ggthemes)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(VennDiagram)
  library(corrplot)
  library(ggdendro)
  library(survival)
  library(survminer)
  library(glmnet)
  library(pheatmap)
  
}

##导入数据
bc <- readRDS("N_bc_seurat.rds")

table(bc@meta.data$cell_type)

GC <- subset(x=bc, cell_type == "GC BC")
table(GC$orig.ident)


#单细胞差异分析
markers_bc <- FindMarkers(GC, min.pct = 0.1, 
                      logfc.threshold = 0,
                      slot = 'data',
                      test.use = 'wilcox',
                      group.by = "orig.ident",
                      ident.1 ="N-ko",
                      random.seed = 1) 


##绘制火山图
DESeq2_result <- markers_bc
DESeq2_result$logP <- -log(DESeq2_result$p_val_adj,base = 20)
DESeq2_result$Group <- 'non-significant'
DESeq2_result$Group[which((DESeq2_result$p_val_adj<0.05) & (DESeq2_result$avg_log2FC>0.25))] <- 'highly expressed in CKO'
DESeq2_result$Group[which((DESeq2_result$p_val_adj<0.05) & (DESeq2_result$avg_log2FC< (-0.25)))] <- 'highly expressed in WT'
table(DESeq2_result$Group)

DESeq2_result$Label <- ''
DESeq2_result <- DESeq2_result[order(DESeq2_result$avg_log2FC,decreasing = TRUE),]
up_genes <- head(rownames(DESeq2_result)[which(DESeq2_result$Group=='highly expressed in CKO')],10)
down_genes <- tail(rownames(DESeq2_result)[which(DESeq2_result$Group=='highly expressed in WT')],10)
DEGs_top10 <- c(up_genes,down_genes)
DESeq2_result$Label[match(DEGs_top10,rownames(DESeq2_result))] <- DEGs_top10


library(ggpubr)

plot1 <- ggscatter(DESeq2_result,
                   x='avg_log2FC',font.family = 'sans',
                   y='logP',
                   color = 'Group',
                   palette = c('red','blue','#BBBBBB'),
                   size = 2,
                   label = 'Label',
                   repel = T,
                   font.label = 8)+theme_base()+
  geom_hline(yintercept = 1,linetype='dashed')+
  geom_vline(xintercept = c(-0.25,0.25),linetype='dashed')+
  theme(legend.position = 'top')+
  geom_point(size = 4, alpha = 0)+
  xlab('Log fold change')+ylab('-Log20(adj P-value)')
plot1

ggsave('GC_deg_volcano.pdf',plot1,width = 9,height = 6)
colon_DEG_Fob <-  filter(DESeq2_result,Group!='non-significant')
write.csv(colon_DEG_Fob,'GC_sig_0.05.csv')
DEGs <- rownames(colon_DEG_Fob)#差异基因


#所有差异基因go富集分析
go_BP <- enrichGO( DEGs,
                   OrgDb = 'org.Mm.eg.db',
                   keyType = 'SYMBOL',
                   ont = "BP",
                   readable = T,
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   minGSSize = 5,
                   maxGSSize = 500)

go_fil <- go_BP@result
go_fil <- go_fil %>% subset(.,qvalue<0.05)
go_fil <- go_fil[order(go_fil$Count,decreasing = T),]
write.csv(go_fil,'GC_BP_alldeg.csv')

##go结果可视化
go_fil2 <- go_fil[1:20,]
datbar <- go_fil2[,c(2,6,9)]
datbar$Term <- factor(datbar$Description,levels = rev(datbar$Description))

datbar <- datbar %>% rename(FDR = p.adjust)

plot2 <- ggplot(datbar,aes(x=Count,y=Term,fill=FDR))+
  geom_bar(stat = 'identity',mapping = aes(fill=FDR))+
  theme_test()+
  scale_colour_gradient(low="#abddff",high="#5066a1",aesthetics = "fill")+
  theme(axis.text = element_text(face='bold',size = 8,angle = 0,family = 'sans',colour = 'black'),
        legend.position = c(0.9,0.5))+
  xlab('Gene Counts')+ylab('')
plot2      
ggsave('GC_alldeg_BP.pdf',plot2,width = 9,height = 6)



#####################################################
##分出KO与WT差异基因
sig_dge.up <- subset(markers_bc, p_val_adj<0.05&avg_log2FC>0.25)
sig_dge.up <- sig_dge.up[order(sig_dge.up$avg_log2FC,decreasing = T),]
write.csv(sig_dge.up,'GC_dge_ko.csv')

sig_dge.down <- subset(markers_bc, p_val_adj<0.05&avg_log2FC< -0.25)
sig_dge.down <- sig_dge.down[order(sig_dge.down$avg_log2FC,decreasing = T),]
write.csv(sig_dge.down,'GC_dge_wt.csv')


##GO
#ko
egoALL <- enrichGO(gene = row.names(sig_dge.up),
                   OrgDb = 'org.Mm.eg.db',
                   keyType = 'SYMBOL',
                   ont = "ALL", #设置为ALL时BP, CC, MF都计算
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05,
                   qvalueCutoff = 0.05,
                   readable = T,
                   minGSSize = 5,
                   maxGSSize = 500)
ego_all <- data.frame(egoALL)

write.csv(egoALL, 'GC_enrichGO_all_ko.csv')
View(ego_all)

#wt_up_deg
egoALL_wt <- enrichGO(gene = row.names(sig_dge.down),
                      OrgDb = 'org.Mm.eg.db',
                      keyType = 'SYMBOL',
                      ont = "ALL", #设置为ALL时BP, CC, MF都计算
                      pAdjustMethod = "BH",
                      pvalueCutoff = 0.05,
                      qvalueCutoff = 0.05,
                      readable = T,
                      minGSSize = 5,
                      maxGSSize = 500)
ego_all_wt <- data.frame(egoALL_wt)
write.csv(egoALL_wt, 'GC_enrichGO_all_wt.csv')

