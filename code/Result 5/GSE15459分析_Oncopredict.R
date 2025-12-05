options(stringsAsFactors = FALSE)
###清空环境
rm(list=ls())
#BiocManager::install("GEOquery")
library(GEOquery)
library(limma)
#library(affy)
#读取验证集数据GSE15459
setwd("G:/免疫亚型分型/result6/GSE15459")
GSE15459 = getGEO('GSE15459', destdir=".", AnnotGPL = F, getGPL = F) 
save(GSE15459,file = 'GSE15459.gset.Rdata')
load("GSE15459.gset.Rdata")
######exp1(GPL570)######
###查看平台信息
GSE15459[["GSE15459_series_matrix.txt.gz"]]@annotation
GPL <- getGEO('GPL570', destdir=".")
GPL1<-Table(GPL)
library(tidyr)
GPL2 <-GPL1 %>% as_tibble() %>% 
  separate_rows(`Gene Symbol`, sep = " /// ")
probe2symbol<-GPL2[,c(1,11)]####gene symbol名字有重复，需要保留第一个基因名
###表达
dat<-exprs(GSE15459[[1]])
#limma背景矫正
dat<- backgroundCorrect(dat, method="normexp")
#limma归一化
dat<- normalizeBetweenArrays(dat, method="quantile")#329
dat<-log2(dat+1)
dat<-as.data.frame(dat)
dat$ID<-rownames(dat)
dat<-merge(dat,probe2symbol,by="ID")
#####去除重复的基因名
dat <- dat[!duplicated(dat$`Gene Symbol`), ]
rownames(dat)<-dat$`Gene Symbol`
dat<-dat[,-c(1,ncol(dat))]
####和marker基因匹配上的基因名保留
clinical<-pData(GSE15459[[1]])
#筛选出有生存的样本
Clinical_outcome<-read.table("Clinical_outcome.txt", 
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE)
Clinical_outcome<-Clinical_outcome[Clinical_outcome$GSM.ID%in%colnames(dat),]
GSE15459_exp<-dat[,Clinical_outcome$GSM.ID]
save(GSE15459_exp,file = "GSE15459_NormExp.rda")

#load("H:/免疫亚型分型/MOVICS/STAD.multiomics.rda")
#STAD.surv<-STAD_multiomics$surv.info
GSE15459_clin.info<-data.frame(fustat=Clinical_outcome$Outcome..1.dead.,
                               futime=Clinical_outcome$Overall.Survival..Months...,
                               Subtype=Clinical_outcome$Subtype,
                               Age=Clinical_outcome$Age_at_surgery,
                               Gender=Clinical_outcome$Gender,
                               Stage=Clinical_outcome$Stage)
rownames(GSE15459_clin.info)<-Clinical_outcome$GSM.ID
GSE15459_clin.info$futime<-GSE15459_clin.info$futime*30
save(GSE15459_clin.info,file = "GSE15459_clin.info.rda")
load("G:/免疫亚型分型/MOVICS/Specific_markers.rda")

#load("G:/免疫亚型分型/MOVICS/Specific_markers.rda")
#NTP预测
library(MOVICS)
getwd()
yau.ntp.pred <- runNTP(expr       = GSE15459_exp,
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
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = GSE15459_clin.info,
                     convt.time       = "d", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR YAU") 

save(yau.ntp.pred, file ="yau.ntp.pred.Rdata")
#load("yau.ntp.pred.Rdata")
######Oncopredict
trainingExprData=readRDS(file="G:/免疫亚型分型/result5/Oncopredict/DataFiles/Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
trainingExprData[1:4,1:4]
trainingPtype=readRDS(file="G:/免疫亚型分型/result5/Oncopredict/DataFiles/Training Data/GDSC2_Res.rds")
trainingPtype<-exp(trainingPtype)
trainingPtype[1:4,1:4]
setwd("G:/免疫亚型分型/result6/GSE15459/Oncopredict")
####加载Oncopredict前需要先加载下面这些包
library(fgsea)
library(SummarizedExperiment)
library(DelayedArray)
library(dbplyr)
library(oncoPredict)
calcPhenotype(trainingExprData=trainingExprData,
              trainingPtype=trainingPtype,
              testExprData=as.matrix(GSE15459_exp),
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
getwd()
