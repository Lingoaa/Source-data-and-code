setwd("E:/all_scRNA")
rm(list=ls())


##加载R包
BiocManager::install("ggsci")

{
  library(Seurat)
  library(dplyr)
  library(reticulate)
  library(sctransform)
  library(cowplot)
  library(ggplot2)
  library(viridis)
  library(tidyr)
  library(magrittr)
  library(reshape2)
  library(readxl)
  library(stringr)
  library(cowplot)
  library(scales)
  library(readr)
  library(gplots)
  library(tibble)
  library(grid)
  library(rlang)
  library(plotrix)
  library(ggsci)
  
}

load()

phe=imm_scRNA@meta.data
colnames(phe)
tb=table(phe$cell_type,
         phe$tissue)

head(tb)

library(gplots)

#统计细胞数量
balloonplot(tb, main ="Immune_cells", xlab ="celltype", ylab="sample",
            label = T, show.margins = T)


bar_data <- as.data.frame(tb)

bar_per <- bar_data %>% 
  group_by(Var2) %>%
  mutate(sum(Freq)) %>%
  mutate(percent = Freq / `sum(Freq)`)
head(bar_per) 
write.csv(bar_per,file = "celltype_by_group_percent.csv")
col =c("#BD0026", "#8EA325","#20B2AA","#8B008B", "#a14462","#E69F00",
       "#0072B2", "#167153", "#223D6C", "#8c510a")
ggplot(bar_per, aes(x = percent, y = Var2)) +
  geom_bar(aes(fill = Var1) , stat = "identity") + coord_flip() +
  theme(axis.ticks = element_line(linetype = "blank"),
        legend.position = "right",
        panel.grid.minor = element_line(colour = NA,linetype = "blank"), 
        panel.background = element_rect(fill = NA),
        plot.background = element_rect(colour = NA)) +
  labs(y = "% Relative cell source", fill = NULL)+labs(x = NULL)+
  scale_fill_manual(values=col)

ggsave("celltype_percent_tissue.png",width = 8,height = 8)


###不同分组中细胞亚群的比例
##分为ko和wt组

table(phe$cell_type)
table(phe$hash.ID)
phe$outcone =  paste(phe$hash.ID,phe$orig.ident,sep = '_')


phe_immune = phe[phe$cell_type %in% c('B cells','T cells','Monocytes','Neutrophils','Macrophages','DC','NK','Mast cells'),]

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
bar_per$Var1 = factor(bar_per$Var1,levels =c('B cells','T cells','Monocytes','Neutrophils',
                                             'Macrophages','DC','NK','Mast cells'))

###拼接7个箱线图
library(gridExtra)
library(ggpubr)

p <- list()

col = c("#F5AE6B","#b0d45d")
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

ggsave("tumor_spleen_boxplot_celltype.png", combined_plot, width = 13,height = 9)