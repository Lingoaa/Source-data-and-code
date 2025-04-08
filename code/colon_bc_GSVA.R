setwd("D:/Desktop/RData/mice_scRNA/GEX+ADT/spleen_identify/B_cell")

##清空环境
rm(list=ls())

##加载R包
library(Seurat)
library(GSVA)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(dplyr)
library(readxl)

##加载数据
bc <- readRDS("bc_seurat_harmony.rds")

table(bc@meta.data$orig.ident)

meta <- bc@meta.data[,c("orig.ident","group")]#分组信息，为了后续作图
bc <- as.matrix(bc@assays$RNA@counts)#提取count矩阵

library(msigdbr)

msigdbr_species() 

mouse <- msigdbr(species = "Mus musculus")

mouse[1:5,1:5]

table(mouse$gs_cat) 

#提取GO文件
table(mouse$gs_subcat)

mouse_GO_bp = msigdbr(species = "Mus musculus",
                      category = "C5", #GO在C5
                      subcategory = "GO:BP") %>% 
  dplyr::select(gs_name,gene_symbol)#这里可以选择gene symbol，也可以选择ID，根据自己数据需求来，主要为了方便
mouse_GO_bp_Set = mouse_GO_bp %>% split(x = .$gene_symbol, f = .$gs_name)#后续gsva要求是list，所以将他转化为list

#表达矩阵数据有了，通路信息有了，就可以进行GSVA分析了

T_gsva <- gsva(expr = bc, 
               gset.idx.list = mouse_GO_bp_Set,
               kcdf="Poisson", #查看帮助函数选择合适的kcdf方法 
               parallel.sz = 5)

write.table(T_gsva, '5bc_gsva.xlsx', row.names=T, col.names=NA, sep="\t")

T_gsva <- read_excel("D:/Desktop/RData/mice_scRNA/GEX+ADT/spleen_identify/B_cell/Fob/Fob_gsva.xlsx")

#用limma包进行差异分析可以
library(edgeR)
library(limma)

table(bc@meta.data$group)

#设置分组
group <- c(rep("WT", 2167), rep("CKO", 189)) %>% as.factor()#设置分组，对照在前
desigN <- model.matrix(~ 0 + group) #构建比较矩阵
colnames(desigN) <- levels(group)

#设置阈值
#logFCcutoff <- log2(1.5)
logFCcutoff <- log2(1.2)
PvalueCutoff <- 0.05      #这里使用未调整p值

#limma通路差异分析:mcao vs sham
fit = lmFit(T_gsva, desigN)
fit2 <- eBayes(fit)
diff=topTable(fit2,adjust='fdr',coef=2,number=Inf)#校准按照需求修改 ？topTable

DEgeneSets <- diff[(diff$P.Value < PvalueCutoff & (diff$logFC>logFCcutoff | diff$logFC < (-logFCcutoff))),]

write.csv(DEgeneSets, file = "5bc_GSVA_pathway_1.2.xlsx")

DEgeneSets <- read_excel("Naive_Cd4_GSVA_DEG.xlsx")

#热图展示
pheatmap::pheatmap(T_gsva[rownames(DEgeneSets,),],
                   show_rownames = T,
                   show_colnames = T)
p = pheatmap::pheatmap(T_gsva[rownames(DEgeneSets,),],
                       show_rownames = T,
                       show_colnames = T)
ggsave('pro_gsva.png', p,width = 25,height = 15)


#最后对差异的感兴趣的通路进行可视化
ko_up <- c("GOBP_PHAGOCYTOSIS_RECOGNITION",
           "GOBP_COMPLEMENT_ACTIVATION",
           "GOBP_HUMORAL_IMMUNE_RESPONSE_MEDIATED_BY_CIRCULATING_IMMUNOGLOBULIN")
wt_up <- c("GOBP_REGULATION_OF_GAP_JUNCTION_ASSEMBLY",
           "GOBP_MESENCHYMAL_STEM_CELL_PROLIFERATION",
           "GOBP_POSITIVE_REGULATION_OF_MESENCHYMAL_STEM_CELL_PROLIFERATION",
           "GOBP_REGULATION_OF_LYMPHANGIOGENESIS",
           "GOBP_MAST_CELL_PROLIFERATION",
           "GOBP_INTERCELLULAR_TRANSPORT",
           "GOBP_NEGATIVE_REGULATION_OF_CELL_FATE_SPECIFICATION",
           "GOBP_POSITIVE_REGULATION_OF_CELL_FATE_COMMITMENT",
           "GOBP_POSITIVE_REGULATION_OF_BICELLULAR_TIGHT_JUNCTION_ASSEMBLY")
TEST <- c(ko_up, wt_up)
diff$ID <- rownames(diff) 
Q <- diff[TEST,]
group1 <- c(rep("ko_up", 3), rep("wt_up", 9)) 
df <- data.frame(ID = Q$ID, score = Q$t,group=group1 )
# 按照t score排序
sortdf <- df[order(df$score),]
sortdf$ID <- factor(sortdf$ID, levels = sortdf$ID)#增加通路ID那一列

ggplot(sortdf, aes(ID, score,fill=group)) + geom_bar(stat = 'identity',alpha = 0.7) + 
  coord_flip() + 
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank())+
  theme(panel.border = element_rect(size = 0.6))+
  labs(x = "",
       y="t value of GSVA score")+
  scale_fill_manual(values = c("red","blue"))#设置颜色

p = ggplot(sortdf, aes(ID, score,fill=group)) + geom_bar(stat = 'identity',alpha = 0.7) + 
  coord_flip() + 
  theme_bw() + #去除背景色
  theme(panel.grid =element_blank())+
  theme(panel.border = element_rect(size = 0.6))+
  labs(x = "",
       y="t value of GSVA score")+
  scale_fill_manual(values = c("#008020","#08519C"))#设置颜色

ggsave('PC_gsva.png', p, width = 20,height = 8)
