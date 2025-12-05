options(stringsAsFactors = FALSE)
###清空环境
rm(list=ls())
#library(affy)
#读取验证集数据GSE91061
setwd("G:/免疫亚型分型/result7/GSE91061")
######GSE91061的平台ID为NCBI entrezID
GSE91061_exp<-read.table("GSE91061_fpkm.txt",
                         header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
####用clusterProfiler进行ID转换
library(clusterProfiler)
ID<-GSE91061_exp$X
idsTOsymbols = bitr(ID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
index<-match(GSE91061_exp$X,idsTOsymbols$ENTREZID)
GSE91061_exp$Symbol<-idsTOsymbols$SYMBOL[index]
GSE91061_exp<-GSE91061_exp[!is.na(GSE91061_exp$Symbol),]
rownames(GSE91061_exp)<-GSE91061_exp$Symbol
GSE91061_exp<-GSE91061_exp[,-c(1,ncol(GSE91061_exp))]
GSE91061_exp<-log2(GSE91061_exp+1)
####治疗之后的样本表达谱预测免疫亚型
pro_exp<-GSE91061_exp[,grep("On",colnames(GSE91061_exp))]
colnames(pro_exp)<-gsub("_.*","",colnames(pro_exp))

###合并整理临床数据
ALL_clinical<-read.table("clinical.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
PD1_clinical<-read.table("G:/免疫亚型分型/result7/GSE91061/Riaz_PD-1/clinical.txt",
                         header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
PD1_CTLA4_clinical<-read.table("G:/免疫亚型分型/result7/GSE91061/Riaz_PD-1_CTLA-4/clinical.txt",
                         header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
clinical<-rbind(PD1_clinical,PD1_CTLA4_clinical)

index<-match(clinical$X,ALL_clinical$Patient)
clinical$RES<-ALL_clinical$Response[index]
####提取表达样本
pro_exp<-pro_exp[,colnames(pro_exp)%in%clinical$X]
index<-match(colnames(pro_exp),clinical$X)
clinical_result<-clinical[index,]

clinical_result$PFS<-clinical_result$PFS/30
clinical_result$fustat<-clinical_result$OS.Event
clinical_result$futime<-clinical_result$OS
table(clinical_result$RES)
clinical_result<-clinical_result[-which(clinical_result$RES%in%"NE"),]
clinical_result$response<-"NDB"
clinical_result$response[clinical_result$RES%in%"CR"|
                           clinical_result$RES%in%"PR"]<-"DCB"
clinical_result$response[clinical_result$RES%in%"SD"|
                           clinical_result$PFS>6]<-"DCB"
table(clinical_result$response)
rownames(clinical_result)<-clinical_result$X
write.table(clinical_result,"G:/免疫亚型分型/result7/GSE91061/clinical_result.txt",row.names=F,sep="\t",quote=FALSE)

pro_exp<-pro_exp[,clinical_result$X]
which(clinical_result$X%in%PD1_clinical$X)
load("G:/免疫亚型分型/MOVICS/Specific_markers.rda")
###marker基因表达谱
marker_exp<-pro_exp[templates$probe,]
marker_exp<-na.omit(marker_exp)
####共有42个maker基因的表达
#Clinical_outcome<-clinical[clinical$Patient%in%colnames(marker_exp),]

#NTP预测
library(MOVICS)
getwd()
yau.ntp.pred <- runNTP(expr       = marker_exp,
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

clust.res<-yau.ntp.pred$clust.res
clust.res$clust<-paste0("CS",clust.res$clust)
index<-match(rownames(clinical_result),clust.res$samID)
clinical_result$cluster<-clust.res$clust[index]
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = clinical_result,
                     convt.time       = "d", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR YAU") 

save(yau.ntp.pred, file ="yau.ntp.pred.Rdata")
######免疫治疗反应结果比例图
#####免疫治疗反应比例分析
library(ggplot2)
GSE91061_Response<-data.frame(Response=clinical_result$response,
                             Cluster=clinical_result$cluster)
fisher_res<-table(GSE91061_Response)
fisher_test<-fisher.test(fisher_res)
fisher_test$p.value
p<-ggplot(data = GSE91061_Response) + 
  ggtitle(paste("Fisher’s exact test P=",fisher_test$p.value,sep = ""))+
  geom_bar(mapping = aes(x = Cluster, fill = Response), position = "fill")+
  scale_fill_manual(values=c("#f29c71","#21498d"))+
  theme_bw()+theme(panel.grid=element_blank())
p
ggsave("Response.pdf",p,width = 5.5,height = 6)

