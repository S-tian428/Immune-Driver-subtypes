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
write.table(GPL1,"GPL570.txt",sep = "\t",quote=FALSE,row.names = F)

####gene symbol名字有重复，需要保留第一个基因名
####和marker基因匹配上的基因名保留
load("G:/免疫亚型分型/MOVICS/Specific_markers.rda")

intersectgene<-unique(GPL1$`Gene Symbol`[which(GPL1$`Gene Symbol`%in%templates$probe)])
#####直接可以匹配的基因及探针
GPL_inter<-GPL1[GPL1$`Gene Symbol`%in%templates$probe,]

diff_gene<-setdiff(templates$probe,intersectgene)
GPL2<-GPL1[grep("///",GPL1$`Gene Symbol`),]
rownames(GPL2)<-NULL
GPL_diff_result<-NULL
i=998
library(stringr)
for (i in 1:nrow(GPL2)) {
  GPL_gene<-unlist(str_split(GPL2$`Gene Symbol`[i],' /// '))
  a<-intersect(GPL_gene,diff_gene)
  if(length(a)>0){
    GPL_diff<-GPL2[i,]
    GPL_diff$`Gene Symbol`<-a[1]
    GPL_diff_result<-rbind(GPL_diff_result,GPL_diff) 
  }
}
GPL_result<-rbind(GPL_inter,GPL_diff_result)
####提取平台探针和基因名文件
probe2symbol_df<-GPL_result[,c(1,11)]
###临床数据(没用)
clinical<-pData(GSE15459[[1]])
#表达谱数据
dat<-exprs(GSE15459[[1]])
#limma背景矫正
dat<- backgroundCorrect(dat, method="normexp")
#limma归一化
dat<- normalizeBetweenArrays(dat, method="quantile")#329
GSE15459_expression<-log2(dat+1)
save(GSE15459_expression,file = "GSE15459_NormExp.rda")
dat<-as.data.frame(GSE15459_expression)
#View(head(dat))
dat$ID<-rownames(dat)
marker_exp<-merge(dat,probe2symbol_df,by="ID")
#####去除重复的基因名
marker_exp <- marker_exp[!duplicated(marker_exp$`Gene Symbol`), ]
rownames(marker_exp)<-marker_exp$`Gene Symbol`
marker_exp<-marker_exp[,-c(1,ncol(marker_exp))]


#筛选出有生存的样本
Clinical_outcome<-read.table("Clinical_outcome.txt", 
                             header = TRUE, sep = "\t", stringsAsFactors = FALSE)
GSE15459.marker.exp<-marker_exp[,colnames(marker_exp)%in%Clinical_outcome$GSM.ID]
save(GSE15459.marker.exp,file="GSE15459.marker.exp.rda")#192个样本
Clinical_outcome<-Clinical_outcome[Clinical_outcome$GSM.ID%in%colnames(GSE15459.marker.exp),]

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
#NTP预测
library(MOVICS)
getwd()
yau.ntp.pred <- runNTP(expr       = GSE15459.marker.exp,
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
######药物反应预测