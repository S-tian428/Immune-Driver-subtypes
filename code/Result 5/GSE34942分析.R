options(stringsAsFactors = FALSE)
###清空环境
rm(list=ls())
#BiocManager::install("GEOquery")
library(GEOquery)
library(limma)
#library(affy)
#读取验证集数据GSE34942
setwd("G:/免疫亚型分型/result6/GSE34942")
GSE34942 = getGEO('GSE34942', destdir=".", AnnotGPL = F, getGPL = F) 
save(GSE34942,file = 'GSE34942.gset.Rdata')
load("GSE34942.gset.Rdata")
######exp1(GPL570)######
###查看平台信息
GSE34942[["GSE34942_series_matrix.txt.gz"]]@annotation
GPL <- getGEO('GPL570', destdir=".")
GPL1<-Table(GPL)
library(tidyr)
GPL2 <-GPL1 %>% as_tibble() %>% 
  separate_rows(`Gene Symbol`, sep = " /// ")
probe2symbol<-GPL2[,c(1,11)]####gene symbol名字有重复，需要保留第一个基因名
###临床数据(没用)
clinical<-pData(GSE34942[[1]])
#表达谱数据
dat<-exprs(GSE34942[[1]])
#limma背景矫正
GSE34942_expression<-log2(dat+1)
GSE34942_expression<-as.data.frame(GSE34942_expression)
#View(head(dat))
GSE34942_expression$ID<-rownames(GSE34942_expression)
GSE34942_expression<-merge(GSE34942_expression,probe2symbol,by="ID")
#####去除重复的基因名
GSE34942_expression <- GSE34942_expression[!duplicated(GSE34942_expression$`Gene Symbol`), ]
rownames(GSE34942_expression)<-GSE34942_expression$`Gene Symbol`
GSE34942_expression<-GSE34942_expression[,-c(1,ncol(GSE34942_expression))]
#筛选出有生存的样本
Clinical_outcome<-read.table("Clinical_outcome.txt", 
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE)
Clinical_outcome<-Clinical_outcome[Clinical_outcome$GSM.ID%in%colnames(GSE34942_expression),]

GSE34942_expression<-GSE34942_expression[,Clinical_outcome$GSM.ID]
save(GSE34942_expression,file = "GSE34942_NormExp.rda")

#load("H:/免疫亚型分型/MOVICS/STAD.multiomics.rda")
#STAD.surv<-STAD_multiomics$surv.info
GSE34942_clin.info<-data.frame(fustat=Clinical_outcome$Outcome..1.dead.,
                               futime=Clinical_outcome$Overall.Survival..Months...,
                               Subtype=Clinical_outcome$Subtype,
                               Age=Clinical_outcome$Age_at_surgery,
                               Gender=Clinical_outcome$Gender,
                               Stage=Clinical_outcome$Stage)
rownames(GSE34942_clin.info)<-Clinical_outcome$GSM.ID
GSE34942_clin.info$futime<-GSE34942_clin.info$futime*30
save(GSE34942_clin.info,file = "GSE34942_clin.info.rda")

load("G:/免疫亚型分型/MOVICS/Specific_markers.rda")
#NTP预测
library(MOVICS)
getwd()
yau.ntp.pred <- runNTP(expr       = GSE34942_expression,
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
                     surv.info        = GSE34942_clin.info,
                     convt.time       = "d", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR YAU") 

save(yau.ntp.pred, file ="yau.ntp.pred.Rdata")
