setwd("E:/spleen_alldata/NT_spleen_BC")
rm(list=ls())


##加载R包
{
  library(Seurat)
  library(cowplot)
  library(harmony)
  library(dplyr)
  library(tidyverse)
  library(patchwork)
  library(BiocManager)
  library(monocle)
  library(patchwork)
}

###导入数据
bc <- readRDS("bc_seurat_harmony.rds")

expr_matrix <- as(as.matrix(bc@assays$RNA@counts), 'sparseMatrix')

p_data <- bc@meta.data

p_data$celltype <- bc@active.ident 

f_data <- data.frame(gene_short_name = row.names(bc),row.names = row.names(bc))

#构建CDS对象
pd <- new('AnnotatedDataFrame', data = p_data)
fd <- new('AnnotatedDataFrame', data = f_data)

#将p_data和f_data从data.frame转换Annotated Data Frame对象。
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())

cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

cds <- detectGenes(cds, min_expr = 3)
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds), num_cells_expressed >= 10))
length(expressed_genes)


###第五步 轨迹定义基因选择及可视化和构建轨迹
#1.使用seurat选择的高变基因
express_genes <- VariableFeatures(bc)
cds <- setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)

#2.使用clusters差异表达基因
deg.cluster <- FindAllMarkers(bc)
express_genes <- subset(deg.cluster,p_val_adj<0.05)$gene
cds <- setOrderingFilter(cds, express_genes)
plot_ordering_genes(cds)

#3.使用monocle选择的高变基因
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table, mean_expression >= 0.1 & dispersion_empirical >= 1 * dispersion_fit)$gene_id
cds <- setOrderingFilter(cds, disp.genes)
plot_ordering_genes(cds)


#4.”dpFeature“法
diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr="~cell_type",cores=1)
head(diff)

deg <- subset(diff, qval < 0.01) #选出2829个基因
deg <- deg[order(deg$qval,decreasing=F),]
head(deg)

ordergene <- rownames(deg)
cds <- setOrderingFilter(cds, ordergene)
plot_ordering_genes(cds)

#选择的用于排序的基因数目一般在2000左右比较合适，gene数太多的话也可以选择top基因
ordergene <- row.names(deg)[order(deg$qval)][1:400]



##第六步降维（首先选择用于细胞排序的基因，然后使用反向图嵌入(DDRTree)算法对数据进行降维）
cds <- reduceDimension(cds, max_components = 2,
                       method = 'DDRTree')

cds <- orderCells(cds)

#使用root_state参数可以设置拟时间轴的根，如下面的拟时间着色图中可以看出，左边是根。根据state图可以看出，根是State1，若要想把另一端设为根，可以按如下操作
cds <- orderCells(cds, root_state = 5) #把State5设成拟时间轴的起始点


##保存meta信息和拟时序对象
saveRDS(cds,file = "bc_monocle_WT_DEG.rds")
write.csv(pData(cds),"bc_monocleMeta_WT_DEG.csv")

##绘图
pData(cds)$State <- as.factor(pData(cds)$State)
pData(cds)$pseudotime <- as.factor(pData(cds)$pseudotime)

p1 <- plot_cell_trajectory(cds, color_by="Pseudotime")
p1
ggsave("Pseudotime.pdf",p1,width = 8,height = 6)


plot_cell_trajectory(cds, color_by="State", show_tree=F, 
                     show_branch_points = F,show_cell_names =F,
                     show_state_number = F,show_backbone = T)

p2 <- plot_cell_trajectory(cds, color_by="cell_type", show_tree=F, 
                           show_branch_points = F,show_cell_names =F,
                           show_state_number = F,show_backbone = T)+
  scale_colour_manual(values = c("#f39c90","#4d62ae","#b37eb7"))
p2
ggsave("Pseudotime_celltype.pdf",p2,width = 8,height = 6)


p3 <- plot_cell_trajectory(cds, color_by="pseudotime", show_tree=F, 
                           show_branch_points = F,show_cell_names =F,
                           show_state_number = F,show_backbone = T)+
  scale_colour_manual(values = c("#c04327","#dd6f56","#e79a88",
                                 "#c491a5","#8c8ec1","#2a4986"))
p3
ggsave("Pseudotime_after.pdf",p3,width = 8,height = 6)

