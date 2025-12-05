##免疫细胞浸润
#####Cibersort
rm(list=ls())
load("I:/免疫亚型分型/MOVICS/STAD.multiomics.rda")
setwd("I:/免疫亚型分型/result3/Cibersort/")
surv.info<-STAD_multiomics$surv.info
#STAD.FPKM<-STAD_multiomics$STAD.FPKM
source('Cibersort.R')
###tpm
result <- CIBERSORT("LM22.txt","STAD.FPKM.txt", perm = 1000, QN = T)
save(result,file = "CIBERSORT_result.rda")
load("CIBERSORT_result.rda")
####整理成矩阵形式
TCGA_cibersort_result<-result[,-(23:25)]
TCGA_cibersort_result<-data.frame(TCGA_cibersort_result)
TCGA_cibersort_result$Cluster<-surv.info$Cluster
library(tidyr)
results <- gather(TCGA_cibersort_result, key = "ImmCells", value = "value",
                  -"Cluster")

##箱线图可视化##
library(dplyr)
library(tibble)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(ggsci)
colnames(results)
p<-ggplot(results,aes(ImmCells,value,fill = Cluster)) +
  geom_boxplot(outlier.shape = 21,color = "black") +
  theme_bw() +
  labs(x = "Cell Type", y = "Estimated Proportion") +
  theme(legend.position = "top") +
  theme(axis.text.x = element_text(angle=80,vjust = 0.5,size = 6,face = "italic",colour = 'black'),
        axis.text.y = element_text(face = "italic",size = 9,colour = 'black'))+
  scale_fill_nejm()+
  stat_compare_means(label = "p.signif",size=3,method = "wilcox.test")##热图展示
ggsave("Cibersort_boxplot.pdf",height = 4.5,width = 6)
##热图展示

heat_map_res<-result[,1:22]
heat_map_res<-t(heat_map_res)

surv.info_cluster<-surv.info[order(surv.info$Cluster),]
heat_map_res<-heat_map_res[,rownames(surv.info_cluster)]

mycol<-colorRampPalette(c("#0000FF", "white", "#FF3030"))(100)
#构建分组信息
annotation_col<-data.frame(clust=surv.info_cluster$Cluster)
rownames(annotation_col) <- rownames(surv.info_cluster)
library(pheatmap)
p<-pheatmap(heat_map_res
         ,color = mycol,
         cluster_rows = T,cluster_cols = F,
         scale = 'row',show_colnames = F,show_rownames = T,
         annotation_col = annotation_col )
ggsave("Cibersort_heatmap.pdf",p,width = 6,height = 2)
##免疫细胞堆叠直方图(意义不大，不画了)
library(IOBR)
p1 <- res %>%
  ggplot(aes(sample,value))+
  geom_bar(stat = "identity",position = "stack",aes(fill=cell.type))+
  labs(x=NULL)+
  scale_y_continuous(expand = c(0,0))+
  scale_fill_manual(values = palette4,name=NULL)+ # iobr还给大家准备了几个色盘，贴心！
  theme_bw()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "bottom"
  )
p1
