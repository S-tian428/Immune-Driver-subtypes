options(stringsAsFactors = FALSE)
###清空环境
rm(list=ls())
#library(affy)
#读取验证集数据GSE91061
setwd("G:/免疫亚型分型/result7/GSE91061/Riaz_PD-1")
######GSE91061的平台ID为NCBI entrezID
####所有样本的基因表达谱
GSE91061_exp<-read.table("G:/免疫亚型分型/result7/GSE91061/GSE91061_fpkm.txt",
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

###只接受抗PD-1治疗患者的临床数据
clinical<-read.table("clinical.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
index<-which(colnames(pro_exp)%in%clinical$X)
###只接受PD-1治疗的样本表达
PD1_expression<-pro_exp[,index]
index<-match(colnames(PD1_expression),clinical$X)
clinical<-clinical[index,]

#####所有患者临床数据
All_clinical<-read.table("G:/免疫亚型分型/result7/GSE91061/clinical.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
PD1_clinical<-All_clinical[All_clinical$Patient%in%colnames(PD1_expression),c(1,3)]
index<-match(colnames(PD1_expression),PD1_clinical$Patient)
PD1_clinical<-PD1_clinical[index,]
PD1_clinical_result<-cbind(clinical,RES=PD1_clinical$Response)

PD1_clinical_result$PFS<-PD1_clinical_result$PFS/30
PD1_clinical_result$fustat<-PD1_clinical_result$OS.Event
PD1_clinical_result$futime<-PD1_clinical_result$OS
table(PD1_clinical_result$RES)
PD1_clinical_result$response<-"NDB"
PD1_clinical_result$response[PD1_clinical_result$RES%in%"CR"|
                            PD1_clinical_result$RES%in%"PR"]<-"DCB"
PD1_clinical_result$response[PD1_clinical_result$RES%in%"SD"|
                               PD1_clinical_result$PFS>6]<-"DCB"
table(PD1_clinical_result$response)
rownames(PD1_clinical_result)<-PD1_clinical_result$X

load("G:/免疫亚型分型/MOVICS/Specific_markers.rda")
###marker基因表达谱
marker_exp<-PD1_expression[templates$probe,]
marker_exp<-na.omit(marker_exp)
####共有42个maker基因的表达
#Clinical_outcome<-clinical[clinical$Patient%in%colnames(marker_exp),]

#load("H:/免疫亚型分型/MOVICS/STAD.multiomics.rda")
#STAD.surv<-STAD_multiomics$surv.info
#save(GSE91061_clin.info,file = "GSE91061_clin.info.rda")
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
index<-match(rownames(PD1_clinical_result),clust.res$samID)
PD1_clinical_result$cluster<-clust.res$clust[index]
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = PD1_clinical_result,
                     convt.time       = "d", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR YAU") 

save(yau.ntp.pred, file ="yau.ntp.pred.Rdata")
######免疫治疗反应结果比例图
####将NDB分为两组
#####免疫治疗反应比例分析
library(ggplot2)
PD1_Response<-data.frame(Response=PD1_clinical_result$response,
                              Cluster=PD1_clinical_result$cluster)
fisher_res<-table(PD1_Response)
fisher_test<-fisher.test(fisher_res)
fisher_test$p.value
p<-ggplot(data = PD1_Response) + 
  ggtitle(paste("Fisher’s exact test P=",fisher_test$p.value,sep = ""))+
  geom_bar(mapping = aes(x = Cluster, fill = Response), position = "fill")+
  scale_fill_manual(values=c("#f29c71","#21498d"))+
  theme_bw()+theme(panel.grid=element_blank())
p
ggsave("Response.pdf",p,width = 5.5,height = 6)

