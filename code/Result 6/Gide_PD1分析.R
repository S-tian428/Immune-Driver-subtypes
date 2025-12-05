options(stringsAsFactors = FALSE)
###清空环境
rm(list=ls())
#library(affy)
#读取验证集数据Gide
setwd("G:/免疫亚型分型/result7/Gide_PD1")
######Gide的平台ID为NCBI entrezID
Gide_exp<-read.table("Expression.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
####用clusterProfiler进行ID转换
library(clusterProfiler)
ID<-Gide_exp$Enrenze.ID
idsTOsymbols = bitr(ID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
index<-match(Gide_exp$Enrenze.ID,idsTOsymbols$ENTREZID)
Gide_exp$Symbol<-idsTOsymbols$SYMBOL[index]
Gide_exp<-Gide_exp[!is.na(Gide_exp$Symbol),]
rownames(Gide_exp)<-Gide_exp$Symbol
Gide_exp<-Gide_exp[,-c(1,ncol(Gide_exp))]
write.table(Gide_exp,"Gide_exp.txt",sep="\t",quote=FALSE)


#Gide_exp<-log2(Gide_exp+1)
####治疗之后的样本表达谱预测免疫亚型
###临床数据
clinical<-read.table("clinical.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
clinical$X<-paste0("X",clinical$X)
# match(colnames(pro_exp),clinical)
index<-match(colnames(Gide_exp),clinical$X)
clinical<-clinical[index,]
clinical$PFS<-clinical$PFS/30
clinical$fustat<-clinical$OS.Event
clinical$futime<-clinical$OS
table(clinical$Best.RECIST.response)
clinical$response<-"NDB"
clinical$response[clinical$Best.RECIST.response%in%"CR"|
                    clinical$Best.RECIST.response%in%"PR"]<-"DCB"
clinical$response[clinical$Best.RECIST.response%in%"SD"|
                    clinical$PFS>6]<-"DCB"

rownames(clinical)<-clinical$X
write.table(clinical,"clinical_result.txt",sep="\t",quote=FALSE)


table(clinical$response)
table(clinical$Response)
###每个样本都有生存数据
####和marker基因匹配上的基因名保留
load("G:/免疫亚型分型/MOVICS/Specific_markers.rda")
###marker基因表达谱
marker_exp<-Gide_exp[templates$probe,]
marker_exp<-na.omit(marker_exp)
####共有44个maker基因的表达
#save(Gide_clin.info,file = "Gide_clin.info.rda")
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
index<-match(rownames(clinical),clust.res$samID)
clinical$cluster<-clust.res$clust[index]
surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = clinical,
                     convt.time       = "d", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR YAU") 


save(yau.ntp.pred, file ="yau.ntp.pred.Rdata")
######免疫治疗反应结果比例图
#####免疫治疗反应比例分析
library(ggplot2)
Gide_Response<-data.frame(Response=clinical$response,
                          Cluster=clinical$cluster)
fisher_res<-table(Gide_Response)
fisher_test<-fisher.test(fisher_res)
fisher_test$p.value
p<-ggplot(data = Gide_Response) + 
  ggtitle(paste("Fisher’s exact test P=",fisher_test$p.value,sep = ""))+
  geom_bar(mapping = aes(x = Cluster, fill = Response), position = "fill")+
  scale_fill_manual(values=c("#f29c71","#21498d"))+
  theme_bw()+theme(panel.grid=element_blank())
p
ggsave("Response.pdf",p,width = 5.5,height = 6)

