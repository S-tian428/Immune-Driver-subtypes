options(stringsAsFactors = FALSE)
###清空环境
rm(list=ls())
# ####GEO外部验证集
setwd("G:/免疫亚型分型/result6/GSE84437/Oncopredict")
resultPtype <- read.csv('./calcPhenotype_Output/DrugPredictions.csv', 
                        header = T ,stringsAsFactors = F ,check.names = F)
colnames(resultPtype)<-gsub("_.*","",colnames(resultPtype))
rownames(resultPtype)<-resultPtype[,1]
resultPtype<-resultPtype[,-1]
load("G:/免疫亚型分型/result6/GSE84437/yau.ntp.pred.Rdata")
####匹配样本亚型信息
clust.res<-yau.ntp.pred$clust.res
clust.res$clust<-paste0("CS",clust.res$clust)
index<-match(rownames(resultPtype),rownames(clust.res))
resultPtype$Cluster<-clust.res$clust
# ####免疫相关药物柱形图
Key_drug<- c("Osimertinib","Gefitinib","Crizotinib",
             "Niraparib", "Olaparib", "Talazoparib","Venetoclax") 
Key_drug_boxplot<-resultPtype[,Key_drug]
Key_drug_boxplot$Cluster<-resultPtype$Cluster
#####
library(ggplot2)
library(ggsci)
library(ggpubr)
setwd("G:/免疫亚型分型/result6/Oncopredict/GSE84437")
j=1
for (j in 1:(ncol(Key_drug_boxplot)-1)) {
  p = ggplot(Key_drug_boxplot, aes(x = Cluster, y = Key_drug_boxplot[,j], fill=Cluster)) + 
    labs(x="", y="IC50",title = colnames(Key_drug_boxplot)[j]) +
    theme_bw(base_size = 9) + ###去除背景颜色
    #geom_text(aes(label=Hugo_Symbol),check_overlap = TRUE,size=6) +
    geom_boxplot(width=0.5, outlier.colour = NA) +
    # 不显示离群点
    geom_jitter(width = 0.1) +
    theme(axis.text = element_text(size = 20))+ 
    #theme(axis.title.x=element_text(vjust=2, size=20,face = "bold"))+
    # 添加散点
    # geom_line(aes(group=Hugo_Symbol),
    #           color="black", alpha=1,linetype=2,
    #           linewidth=0.8) +
    theme(axis.title.y=element_text(vjust=2, size=20))+
    scale_fill_manual(values = c('#fdd692','#ec7357'))+
    stat_compare_means( method = "wilcox.test")
  p
  ggsave(paste0(colnames(Key_drug_boxplot)[j],".pdf"),p,height = 6,width = 6)
}

options(stringsAsFactors = FALSE)
###清空环境
rm(list=ls())
# ####TCGA
setwd("G:/免疫亚型分型/result5/Oncopredict")
resultPtype <- read.csv('./calcPhenotype_Output/DrugPredictions.csv', 
                        header = T ,stringsAsFactors = F ,check.names = F)
colnames(resultPtype)<-gsub("_.*","",colnames(resultPtype))
rownames(resultPtype)<-resultPtype[,1]
resultPtype<-resultPtype[,-1]
load("G:/免疫亚型分型/MOVICS/cmoic.brca.rda")
####匹配样本亚型信息
clust.res<-cmoic.brca$clust.res
clust.res$clust<-paste0("CS",clust.res$clust)
index<-match(rownames(resultPtype),rownames(clust.res))
resultPtype$Cluster<-clust.res$clust
# ####免疫相关药物柱形图
Key_drug<- c("Camptothecin","Osimertinib","Gefitinib","Crizotinib",
             "Niraparib", "Olaparib", "Talazoparib","Venetoclax") 
Key_drug_boxplot<-resultPtype[,Key_drug]
Key_drug_boxplot$Cluster<-resultPtype$Cluster
#####
library(ggplot2)
library(ggsci)
library(ggpubr)
setwd("G:/免疫亚型分型/result6/Oncopredict/TCGA")
j=1
for (j in 1:(ncol(Key_drug_boxplot)-1)) {
  p = ggplot(Key_drug_boxplot, aes(x = Cluster, y = Key_drug_boxplot[,j], fill=Cluster)) + 
    labs(x="", y="IC50",title = colnames(Key_drug_boxplot)[j]) +
    theme_bw(base_size = 9) + ###去除背景颜色
    #geom_text(aes(label=Hugo_Symbol),check_overlap = TRUE,size=6) +
    geom_boxplot(width=0.5, outlier.colour = NA) +
    # 不显示离群点
    geom_jitter(width = 0.1) +
    theme(axis.text = element_text(size = 20))+ 
    #theme(axis.title.x=element_text(vjust=2, size=20,face = "bold"))+
    # 添加散点
    # geom_line(aes(group=Hugo_Symbol),
    #           color="black", alpha=1,linetype=2,
    #           linewidth=0.8) +
    theme(axis.title.y=element_text(vjust=2, size=20))+
    scale_fill_manual(values = c('#fdd692','#ec7357'))+
    stat_compare_means( method = "wilcox.test")
  p
  ggsave(paste0(colnames(Key_drug_boxplot)[j],".pdf"),p,height = 6,width = 6)
}

