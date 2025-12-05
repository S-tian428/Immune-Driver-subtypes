options(stringsAsFactors = FALSE)
###清空环境
rm(list=ls())
library(GEOquery)
library(limma)
#library(affy)
#读取验证集数据GSE84437
setwd("G:/免疫亚型分型/result6/GSE84437")
GSE84437 = getGEO('GSE84437', destdir=".", AnnotGPL = F, getGPL = F) 
save(GSE84437,file = 'GSE84437.gset.Rdata')
load("GSE84437.gset.Rdata")
######exp1(GPL8432)######
###查看平台信息
GSE84437[["GSE84437_series_matrix.txt.gz"]]@annotation
# options('download.file.method.GEOquery'='libcurl')
# options('download.file.method.GEOquery'='auto') 
# unlink("./GPL8432.soft")
GPL <- getGEO('GPL6947', destdir=".")
GPL1<-Table(GPL)
write.table(GPL1,"GPL6947.txt",sep = "\t",quote=FALSE,row.names = F)
probe2symbol_allgene<-GPL1[,c(1,14)]
###临床数据
clinical<-pData(GSE84437[[1]])
#表达谱数据
dat<-exprs(GSE84437[[1]])
dat<-as.data.frame(dat)
#View(head(dat))
dat$ID<-rownames(dat)
###所有基因表达谱
dat<-merge(dat,probe2symbol_allgene,by="ID")
#####去除重复的基因名
dat <- dat[!duplicated(dat$Symbol), ]
dat<-dat[-1,]
rownames(dat)<-dat$Symbol
####共有36个maker基因的表达
dat<-dat[,-c(1,ncol(dat))]
#limma背景矫正
dat<- backgroundCorrect(dat, method="normexp")
# #limma归一化
dat<- normalizeBetweenArrays(dat, method="quantile")#329
#####将NA设置成0
GSE84437_expression<-log2(dat+1)
#####标准化后的表达谱
save(GSE84437_expression,file = "GSE84437_NormExp.rda")


#筛选出有生存的样本
Clinical_outcome<-clinical[clinical$geo_accession%in%colnames(GSE84437_expression),]
GSE84437_expression<-GSE84437_expression[,Clinical_outcome$geo_accession]
#GSE84437.marker.exp<-GSE84437.marker.exp[,colSums(GSE84437.marker.exp)>0]

save(GSE84437_expression,file="GSE84437.surv.exp.rda")#192个样本

GSE84437_clin.info<-data.frame(futime=Clinical_outcome$`duration overall survival:ch1`,
                               fustat=Clinical_outcome$`death:ch1`
)
rownames(GSE84437_clin.info)<-Clinical_outcome$geo_accession
#GSE84437_clin.info<-GSE84437_clin.info[!is.na(GSE84437_clin.info$futime),]
save(GSE84437_clin.info,file = "GSE84437_clin.info.rda")
#NTP预测
library(MOVICS)
getwd()
yau.ntp.pred <- runNTP(expr       = GSE84437_expression,
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

GSE84437_clin.info$futime<-as.numeric(GSE84437_clin.info$futime)
GSE84437_clin.info$futime<-GSE84437_clin.info$futime*30
GSE84437_clin.info$fustat<-as.numeric(GSE84437_clin.info$fustat)

surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = GSE84437_clin.info,
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
setwd("G:/免疫亚型分型/result6/GSE84437/Oncopredict")
####加载Oncopredict前需要先加载下面这些包
library(fgsea)
library(SummarizedExperiment)
library(DelayedArray)
library(dbplyr)
library(oncoPredict)
calcPhenotype(trainingExprData=trainingExprData,
              trainingPtype=trainingPtype,
              testExprData=as.matrix(GSE84437_expression),
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
