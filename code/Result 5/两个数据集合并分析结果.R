options(stringsAsFactors = FALSE)
###清空环境
rm(list=ls())
library(GEOquery)
setwd("G:/免疫亚型分型/result6/GSE34942_GSE15459")
#########NTP########
load("G:/免疫亚型分型/result6/GSE34942/GSE34942_NormExp.rda")
GSE34942_expression<-dat
load("G:/免疫亚型分型/result6/GSE15459/GSE15459_NormExp.rda")
####合并表达谱
index1<-match(rownames(GSE34942_expression),rownames(GSE15459_expression))
ALL_expression<-cbind(GSE34942_expression,GSE15459_expression[index1,])
#library(GEOquery)
GPL <- getGEO('GPL570', destdir=".")
GPL1<-Table(GPL)
#GPL1$`Gene Symbol`
library(tidyr)
GPL2 <-GPL1 %>% as_tibble() %>% 
  separate_rows(`Gene Symbol`, sep = " /// ")
probe2symbol<-GPL2[,c(1,11)]
ALL_expression<-as.data.frame(ALL_expression)
ALL_expression$ID<-rownames(ALL_expression)
ALL_expression<-merge(ALL_expression,probe2symbol,by="ID")
####提取平台探针和基因名文件
###临床数据(没用)
#####去除重复的基因名
ALL_expression <- ALL_expression[!duplicated(ALL_expression$`Gene Symbol`), ]
rownames(ALL_expression)<-ALL_expression$`Gene Symbol`
ALL_expression<-ALL_expression[,-c(1,ncol(ALL_expression))]
dim(ALL_expression)
#ALL_expression<-log2(ALL_expression+1)
####临床
load("G:/免疫亚型分型/result6/GSE34942/GSE34942_clin.info.rda")
GSE34942_clin.info$sample_type<-"GSE34942"
load("G:/免疫亚型分型/result6/GSE15459/GSE15459_clin.info.rda")
GSE15459_clin.info$sample_type<-"GSE15459"
clinical<-rbind(GSE34942_clin.info,GSE15459_clin.info)
clinical$sample<-rownames(clinical)

###临床样本和表达应该对应
ALL_expression<-ALL_expression[,clinical$sample]

###去批次
#install.packages("FactoMineR")
library(FactoMineR)
#install.packages("factoextra")
library(factoextra)
pre.pca <- PCA(t(ALL_expression),graph = FALSE)
p<-fviz_pca_ind(pre.pca,
             geom= c("point", "text"),
             col.ind = clinical$sample_type,
             addEllipses = TRUE,
             legend.title="Group"  )
p
ggsave("pre_PCA.pdf",p,width = 7,height = 6)
####WIDTH=7.HEIGHT=6
###去批次
library(sva)
expr_combat <- ComBat(dat = ALL_expression, batch = clinical$sample_type)
###ComBat_seq函数要求输入的数据为matrix，如果是数据框，则会返回一个不知道是啥的list，所以这里要注意
###batch和group参数都是一个向量，可以直接设置向量，但是要注意与count样本顺序保持一致
combat.pca <- PCA(t(expr_combat),graph = FALSE)
p<-fviz_pca_ind(combat.pca,col.ind = clinical$sample_type,
             geom = c("point", "text"),
             addEllipses = TRUE,
             legend.title="Group" )
p
ggsave("ComBat_PCA.pdf",p,width = 7,height = 6)
####WIDTH=7.HEIGHT=6
load("G:/免疫亚型分型/MOVICS/Specific_markers.rda")
expr_combat<-log2(expr_combat+1)
marker.exp<-expr_combat[rownames(expr_combat)%in%templates$probe,]
#NTP预测
library(MOVICS)
getwd()
yau.ntp.pred <- runNTP(expr       = marker.exp,
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
                     surv.info        = clinical,
                     convt.time       = "d", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR YAU") 

save(yau.ntp.pred, file ="yau.ntp.pred.Rdata")
table(yau.ntp.pred$clust.res)
agree.yau <- compAgree(moic.res  = yau.ntp.pred,
                       subt2comp = GSE15459_clin.info[, "Subtype", drop = FALSE],
                       doPlot    = TRUE,
                       fig.name  = "YAU PREDICTEDMOIC WITH Subtype")
rm(list=ls())
#########NTP########
#########NTP########
setwd("H:/免疫亚型分型/外部验证数据/STAD/GSE34942")
load("GSE34942.marker.exp.rda")
load("GSE34942_clin.info.rda")
load("H:/免疫亚型分型/MOVICS/Specific_markers.rda")
#NTP预测
library(MOVICS)
yau.ntp.pred <- runNTP(expr       = GSE34942.marker.exp,
                       templates  = templates, # the template has been already prepared in runMarker()
                       scaleFlag  = TRUE, # scale input data (by default)
                       centerFlag = TRUE, # center input data (by default)
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR YAU1")
save(yau.ntp.pred, file ="yau.ntp.pred.Rdata")

ntp.res<-yau.ntp.pred$ntp.res
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = GSE34942_clin.info,
                     convt.time       = "m", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR YAU") 
table(yau.ntp.pred$clust.res)
###检查预测的亚型和已有亚型分类之间的一致性。
agree.yau <- compAgree(moic.res  = yau.ntp.pred,
                       subt2comp = GSE34942_clin.info[, "Subtype", drop = FALSE],
                       doPlot    = TRUE,
                       fig.name  = "YAU PREDICTEDMOIC WITH Subtype")

rm(list=ls())
#BiocManager::install("GEOquery")
library(GEOquery)
library(limma)
#library(affy)
load("H:/免疫亚型分型/外部验证数据/STAD/GSE34942/GSE34942.marker.exp.rda")
load("H:/免疫亚型分型/外部验证数据/STAD/GSE15459/GSE15459.marker.exp.rda")
load("H:/免疫亚型分型/外部验证数据/STAD/GSE34942/GSE34942_clin.info.rda")
load("H:/免疫亚型分型/外部验证数据/STAD/GSE15459/GSE15459_clin.info.rda")
GEO_exp<-cbind(GSE34942.marker.exp,GSE15459.marker.exp)
GEO_clinical<-rbind(GSE34942_clin.info,GSE15459_clin.info)
#GSE15459.marker.exp<-scale(GSE15459.marker.exp)
load("H:/免疫亚型分型/MOVICS/Specific_markers.rda")
#NTP预测
library(MOVICS)
setwd("H:/免疫亚型分型/外部验证数据/STAD/Total_GEO")
GEO_exp<-scale(GEO_exp)
yau.ntp.pred <- runNTP(expr       = GEO_exp,
                       templates  = templates, # the template has been already prepared in runMarker()
                       scaleFlag  =T, # scale input data (by default)
                       centerFlag = T, # center input data (by default)
                       distance = "cosine",
                       ####用"cosine", "pearson", "spearman", "kendall"
                       ###方法分别试一下
                       doPlot     = TRUE, # to generate heatmap
                       fig.name   = "NTP HEATMAP FOR YAU1",
                       height = 5,
                       width=5)
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = GEO_clinical,
                     convt.time       = "m", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR YAU") 
save(yau.ntp.pred, file ="yau.ntp.pred.Rdata")

agree.yau <- compAgree(moic.res  = yau.ntp.pred,
                       subt2comp = GSE15459_clin.info[, "Subtype", drop = FALSE],
                       doPlot    = TRUE,
                       fig.name  = "YAU PREDICTEDMOIC WITH Subtype")




