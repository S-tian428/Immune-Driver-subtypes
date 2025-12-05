options(stringsAsFactors = FALSE)
###清空环境
rm(list=ls())
library(GEOquery)
library(limma)
#library(affy)
#读取验证集数据GSE26253
setwd("G:/免疫亚型分型/result6/GSE26253")
GSE26253 = getGEO('GSE26253', destdir=".", AnnotGPL = F, getGPL = F) 
save(GSE26253,file = 'GSE26253.gset.Rdata')
load("GSE26253.gset.Rdata")
######exp1(GPL8432)######
###查看平台信息
GSE26253[["GSE26253_series_matrix.txt.gz"]]@annotation
# options('download.file.method.GEOquery'='libcurl')
# options('download.file.method.GEOquery'='auto') 
# unlink("./GPL8432.soft")
GPL <- getGEO('GPL8432', destdir=".")
GPL1<-Table(GPL)
write.table(GPL1,"GPL8432.txt",sep = "\t",quote=FALSE,row.names = F)
probe2symbol_allgene<-GPL1[,c(1,12)]
###临床数据(没用)
clinical<-pData(GSE26253[[1]])
#表达谱数据
dat<-exprs(GSE26253[[1]])
# #limma背景矫正
# dat<- backgroundCorrect(dat, method="normexp")
# #limma归一化
# dat<- normalizeBetweenArrays(dat, method="quantile")#329
GSE26253_expression<-log2(dat+1)
GSE26253_expression<-as.data.frame(GSE26253_expression)
GSE26253_expression$ID<-rownames(GSE26253_expression)
###所有基因表达谱
GSE26253_expression<-merge(GSE26253_expression,probe2symbol_allgene,by="ID")
#####去除重复的基因名
GSE26253_expression <- GSE26253_expression[!duplicated(GSE26253_expression$Symbol), ]
rownames(GSE26253_expression)<-GSE26253_expression$Symbol
####共有36个maker基因的表达
GSE26253_expression<-GSE26253_expression[,-c(1,ncol(GSE26253_expression))]
#####标准化后的表达谱
Clinical_outcome<-clinical[clinical$geo_accession%in%colnames(GSE26253_expression),]
GSE26253_expression<-GSE26253_expression[,Clinical_outcome$geo_accession]
save(GSE26253_expression,file = "GSE26253_NormExp.rda")
#load("H:/免疫亚型分型/MOVICS/STAD.multiomics.rda")
#STAD.surv<-STAD_multiomics$surv.info
GSE26253_clin.info<-data.frame(Stage=Clinical_outcome$`pathological stage:ch1`,
                               futime=Clinical_outcome$`recurrence free survival time (month):ch1`,
                               fustat=Clinical_outcome$`status (0=non-recurrence, 1=recurrence):ch1`
)
rownames(GSE26253_clin.info)<-Clinical_outcome$geo_accession
save(GSE26253_clin.info,file = "GSE26253_clin.info.rda")
#NTP预测
load("G:/免疫亚型分型/MOVICS/Specific_markers.rda")

library(MOVICS)
getwd()
yau.ntp.pred <- runNTP(expr       = GSE26253_expression,
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
str(GSE26253_clin.info)
GSE26253_clin.info$futime<-as.numeric(GSE26253_clin.info$futime)
GSE26253_clin.info$futime<-GSE26253_clin.info$futime*30
GSE26253_clin.info$fustat<-as.numeric(GSE26253_clin.info$fustat)

surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = GSE26253_clin.info,
                     convt.time       = "d", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR YAU") 

save(yau.ntp.pred, file ="yau.ntp.pred.Rdata")
table(yau.ntp.pred$clust.res)
######Oncopredict
trainingExprData=readRDS(file="G:/免疫亚型分型/result5/Oncopredict/DataFiles/Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
trainingExprData[1:4,1:4]
trainingPtype=readRDS(file="G:/免疫亚型分型/result5/Oncopredict/DataFiles/Training Data/GDSC2_Res.rds")
trainingPtype<-exp(trainingPtype)
trainingPtype[1:4,1:4]
setwd("G:/免疫亚型分型/result6/GSE26253/Oncopredict")
####加载Oncopredict前需要先加载下面这些包
library(fgsea)
library(SummarizedExperiment)
library(DelayedArray)
library(dbplyr)
library(oncoPredict)
calcPhenotype(trainingExprData=trainingExprData,
              trainingPtype=trainingPtype,
              testExprData=as.matrix(GSE26253_expression),
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
