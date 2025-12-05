options(stringsAsFactors = FALSE)
###清空环境
rm(list=ls())
library(GEOquery)
library(limma)
#library(affy)
#读取验证集数据GSE62254
setwd("G:/免疫亚型分型/result6/GSE62254")
GSE62254 = getGEO('GSE62254', destdir=".", AnnotGPL = F, getGPL = F) 
save(GSE62254,file = 'GSE62254.gset.Rdata')
load("GSE62254.gset.Rdata")
######exp1(GPL8432)######
###查看平台信息
GSE62254[["GSE62254_series_matrix.txt.gz"]]@annotation
# options('download.file.method.GEOquery'='libcurl')
# options('download.file.method.GEOquery'='auto') 
# unlink("./GPL8432.soft")
GPL <- getGEO('GPL570', destdir=".")
GPL1<-Table(GPL)
write.table(GPL1,"GPL570.txt",sep = "\t",quote=FALSE,row.names = F)
library(tidyr)
GPL2 <-GPL1 %>% as_tibble() %>% 
  separate_rows(`Gene Symbol`, sep = " /// ")
probe2symbol<-GPL2[,c(1,11)]####gene symbol名字有重复，需要保留第一个基因名

###临床数据(没用)
clinical<-read.table("Survival clinical.txt",
                     header = T,sep="\t",
                     stringsAsFactors=FALSE,quote = "",
                     fill = T)
clinical$OS.status[clinical$OS.status==0|clinical$OS.status==1]<-0
clinical$OS.status[clinical$OS.status!=0]<-1
clinical$Tumor.ID<-paste0("T",clinical$Tumor.ID)
####匹配GEO ID
GEO_clinical<-pData(GSE62254[[1]])
index<-match(clinical$Tumor.ID,GEO_clinical$title)
clinical$geo_accession<-GEO_clinical$geo_accession[index]
#表达谱数据
dat<-exprs(GSE62254[[1]])
dat<-as.data.frame(dat)
#View(head(dat))
dat$ID<-rownames(dat)
###所有基因表达谱
GSE62254_exp<-merge(dat,probe2symbol,by="ID")
#####去除重复的基因名
GSE62254_exp <- GSE62254_exp[!duplicated(GSE62254_exp$`Gene Symbol`), ]
GSE62254_exp<-GSE62254_exp[,-1]
library(tidyverse) 
rownames(GSE62254_exp)<-NULL
GSE62254_exp<-column_to_rownames(GSE62254_exp, "Gene Symbol")

Clinical_outcome<-clinical[clinical$geo_accession%in%colnames(GSE62254_exp),]
GSE62254_exp<-GSE62254_exp[,Clinical_outcome$geo_accession]
#####标准化后的表达谱
save(GSE62254_exp,file = "GSE62254_NormExp.rda")
colnames(Clinical_outcome)[3]<-"fustat"
colnames(Clinical_outcome)[5]<-"futime"
Clinical_outcome$futime<-Clinical_outcome$futime*30
rownames(Clinical_outcome)<-Clinical_outcome$geo_accession
save(Clinical_outcome,file = "GSE62254_clin.info.rda")

load("G:/免疫亚型分型/MOVICS/Specific_markers.rda")
#NTP预测
library(MOVICS)
getwd()
yau.ntp.pred <- runNTP(expr       = GSE62254_exp,
                       templates  = templates, # the template has been already prepared in runMarker()
                       scaleFlag  =T, # scale input data (by default)
                       centerFlag = T, # center input data (by default)
                       distance = "spearman",
                       ####用"cosine", "pearson", "spearman", "kendall"
                       ###方法分别试一下
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR YAU1",
                       height = 5,
                       width=5)
str(Clinical_outcome)
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = Clinical_outcome,
                     convt.time       = "d", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR YAU") 

save(yau.ntp.pred, file ="yau.ntp.pred.Rdata")
load("yau.ntp.pred.Rdata")
table(yau.ntp.pred$clust.res)
######Oncopredict
trainingExprData=readRDS(file="G:/免疫亚型分型/result5/Oncopredict/DataFiles/Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
trainingExprData[1:4,1:4]
trainingPtype=readRDS(file="G:/免疫亚型分型/result5/Oncopredict/DataFiles/Training Data/GDSC2_Res.rds")
trainingPtype<-exp(trainingPtype)
trainingPtype[1:4,1:4]
setwd("G:/免疫亚型分型/result6/GSE62254/Oncopredict")
####加载Oncopredict前需要先加载下面这些包
library(fgsea)
library(SummarizedExperiment)
library(DelayedArray)
library(dbplyr)
library(oncoPredict)
calcPhenotype(trainingExprData=trainingExprData,
              trainingPtype=trainingPtype,
              testExprData=as.matrix(GSE62254_exp),
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
resultPtype <- read.csv('./calcPhenotype_Output/DrugPredictions.csv', 
                        header = T ,stringsAsFactors = F ,check.names = F)
colnames(resultPtype)<-gsub("_.*","",colnames(resultPtype))
rownames(resultPtype)<-resultPtype[,1]
resultPtype<-resultPtype[,-1]
####匹配样本亚型信息
clust.res<-yau.ntp.pred$clust.res
clust.res$clust<-paste0("CS",clust.res$clust)
index<-match(rownames(resultPtype),rownames(clust.res))
resultPtype$Cluster<-clust.res$clust
i=1
wilcox_test_result<-NULL
for (i in 1:(ncol(resultPtype)-1)) {
  wt<- wilcox.test(resultPtype[,i] ~ Cluster, data = resultPtype)
  wt_p<-data.frame(Drug=colnames(resultPtype)[i],P=wt$p.value)
  wilcox_test_result<-rbind(wilcox_test_result,wt_p)
}

wilcox_test_result<-wilcox_test_result[wilcox_test_result$P<0.05,]
a<-aggregate(resultPtype[,1:198],by=list(resultPtype$Cluster),mean)
rownames(a)<-a[,1]
a<-a[,-1]
IC50_mean<-t(a)
IC50_mean<-data.frame(IC50_mean,stringsAsFactors = F)
index1<-match(wilcox_test_result$Drug,rownames(IC50_mean))
result<-cbind(wilcox_test_result,IC50_mean[index1,])
str(result)
result$Difference=result$CS2-result$CS1
result$response_cluster<-"CS1"
index_CS2<-which(result$Difference<0)
result$response_cluster[index_CS2]<-"CS2"
write.table(result,"Wilcox_test.txt",sep="\t",quote=FALSE,row.names = F)

resultPtype_box<-resultPtype[,result$Drug]
resultPtype_box<-cbind(resultPtype_box,Cluster=resultPtype$Cluster)

#####
library(ggplot2)
library(ggsci)
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
    stat_compare_means( method = "wilcox.test")
  p
  ggsave(paste0(colnames(resultPtype_box)[j],".pdf"),p,height = 6,width = 6)
}
