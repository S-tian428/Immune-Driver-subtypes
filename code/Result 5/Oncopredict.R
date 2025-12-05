rm(list=ls())
####亚型分类结果
load("H:/免疫亚型分型/result6/GSE34942_GSE15459/yau.ntp.pred.Rdata")
####数据集表达结果
load("H:/免疫亚型分型/result6/GSE34942/GSE34942_NormExp.rda")
load("H:/免疫亚型分型/result6/GSE15459/GSE15459_NormExp.rda")
GEO_exp<-cbind(GSE34942_expression,GSE15459_expression)
GPL <- getGEO('GPL570', destdir=".")
GPL1<-Table(GPL)
probe2symbol_ALL<-GPL1[,c(1,11)]
dat<-as.data.frame(GEO_exp)
#View(head(dat))
dat$ID<-rownames(dat)
gene_exp<-merge(dat,probe2symbol_ALL,by="ID")
#gene_exp$`Gene Symbol`
####将一行拆分成多行
library(tidyr)
gene_exp <-gene_exp %>% as_tibble() %>% 
  separate_rows(`Gene Symbol`, sep = " /// ")
#####去除重复的基因名
gene_exp <- gene_exp[!duplicated(gene_exp$`Gene Symbol`), ]
gene_exp<-gene_exp[,-1]
library(tidyverse) 
rownames(gene_exp)<-NULL
gene_exp<-column_to_rownames(gene_exp, "Gene Symbol") 
subtype_sample<-yau.ntp.pred$clust.res
#BiocManager::install("GEOquery")
#library(GEOquery)
#####oncoPredict
library(oncoPredict)
trainingExprData=readRDS(file="H:/免疫亚型分型/result5/Oncopredict/DataFiles/Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
trainingExprData[1:4,1:4]
trainingPtype=readRDS(file="H:/免疫亚型分型/result5/Oncopredict/DataFiles/Training Data/GDSC2_Res.rds")
trainingPtype<-exp(trainingPtype)
setwd("H:/免疫亚型分型/result6/GSE34942_GSE15459/Oncopredict/")
calcPhenotype(trainingExprData=trainingExprData,
              trainingPtype=trainingPtype,
              testExprData=as.matrix(gene_exp),
              batchCorrect="eb",
              powerTransformPhenotype=TRUE,
              removeLowVaryingGenes=0.2,
              minNumSamples=10,
              selection=1,
              printOutput=TRUE,
              pcr=FALSE,
              removeLowVaringGenesFrom="homogenizeData",
              report_pc=FALSE,
              cc=FALSE,
              percent=80,
              rsq=FALSE)
#####读取结果
resultPtype <- read.csv('./calcPhenotype_Output/DrugPredictions.csv', 
                        header = T ,stringsAsFactors = F ,check.names = F)
colnames(resultPtype)<-gsub("_.*","",colnames(resultPtype))
rownames(resultPtype)<-resultPtype[,1]
resultPtype<-resultPtype[,-1]
index<-match(rownames(resultPtype),subtype_sample$samID)
resultPtype$Cluster<-subtype_sample$clust[index]
resultPtype<-resultPtype[!is.na(resultPtype$Cluster),]
i=1
wilcox_test_result<-NULL
for (i in 1:(ncol(resultPtype)-1)) {
  wt<- wilcox.test(resultPtype[,i] ~ Cluster, data = resultPtype)
  wt_p<-data.frame(Drug=colnames(resultPtype)[i],P=wt$p.value)
  wilcox_test_result<-rbind(wilcox_test_result,wt_p)
}
wilcox_test_result<-wilcox_test_result[wilcox_test_result$P<0.05,]
a<-aggregate(resultPtype[,1:198],by=list(resultPtype$Cluster),mean)
rownames(a)<-paste0("CS",  a[,1])
a<-a[,-1]
IC50_mean<-t(a)
IC50_mean<-data.frame(IC50_mean,stringsAsFactors = F)
index1<-match(wilcox_test_result$Drug,rownames(IC50_mean))
result<-cbind(wilcox_test_result,IC50_mean[index1,])
result$Difference=result$CS2-result$CS1
result$response_cluster<-"CS1"
index_CS2<-which(result$Difference<0)
result$response_cluster[index_CS2]<-"CS2"
write.table(result,"Wilcox_test.txt",sep="\t",quote=FALSE,row.names = F)
####免疫相关药物柱形图
Immune_therapy_drugs <- c( "Nutlin-3a (-)", 
                         "Talazoparib","Venetoclax",
                          "JQ1")
index<-which(result$Drug%in%Immune_therapy_drugs)
Imm_result<-result[index,]
boxplot_drug<-Imm_result$Drug[Imm_result$response_cluster%in%"CS2"]
resultPtype_box<-resultPtype[,boxplot_drug]
resultPtype_box$Cluster<-resultPtype$Cluster
#####
library(ggplot2)
#install.packages("ggsci")
library(ggsci)
#install.packages("ggpubr")
library(ggpubr)
j=1
for (j in 1:(ncol(resultPtype_box)-1)) {
  p = ggplot(resultPtype_box, aes(x = Cluster, y = resultPtype_box[,j], fill=Cluster)) + 
    labs(x="", y="IC50",title = colnames(resultPtype_box)[j]) +
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
    scale_fill_manual(values = c('#fdd692','#0f7a3f'))+
    stat_compare_means( method = "wilcox.test",size=7)
  p
  ggsave(paste0(colnames(resultPtype_box)[j],".pdf"),p,height = 6,width = 6)
}
####经典传统放射治疗
Radio_therapy_drugs<-c("Cytarabine","Docetaxel",
                       "5-Fluorouracil","Paclitaxel")
resultPtype_box_Radio<-resultPtype[,Radio_therapy_drugs]
resultPtype_box_Radio$Cluster<-resultPtype$Cluster

#####
library(ggplot2)
library(ggsci)
library(ggpubr)
j=1
for (j in 1:(ncol(resultPtype_box_Radio)-1)) {
  p = ggplot(resultPtype_box_Radio, aes(x = Cluster, y = resultPtype_box_Radio[,j], fill=Cluster)) + 
    labs(x="", y="IC50",title = colnames(resultPtype_box_Radio)[j]) +
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
    scale_fill_manual(values = c('#fdd692','#0f7a3f'))+
    stat_compare_means( method = "wilcox.test",size=7)
  p
  ggsave(paste0(colnames(resultPtype_box_Radio)[j],".pdf"),p,height = 6,width = 6)
}
