options(stringsAsFactors = FALSE)
###清空环境
rm(list=ls())
#library(affy)
#读取验证集数据Gide
setwd("G:/免疫亚型分型/result7/Nathanson_CTLA-4")
######Gide的平台ID为NCBI entrezID
Nathanson_exp<-read.table("Expression.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
####用clusterProfiler进行ID转换
library(clusterProfiler)
ID<-Nathanson_exp$Entrez
idsTOsymbols = bitr(ID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
index<-match(Nathanson_exp$Entrez,idsTOsymbols$ENTREZID)
Nathanson_exp$Symbol<-idsTOsymbols$SYMBOL[index]
Nathanson_exp<-Nathanson_exp[!is.na(Nathanson_exp$Symbol),]
rownames(Nathanson_exp)<-Nathanson_exp$Symbol
Nathanson_exp<-Nathanson_exp[,-c(1,ncol(Nathanson_exp))]
write.table(Nathanson_exp,"Nathanson_exp.txt",sep="\t",quote=FALSE)

#Nathanson_exp<-log2(Nathanson_exp+1)
####治疗之后的样本表达谱预测免疫亚型
###临床数据
clinical<-read.table("clinical.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
Nathanson_exp<-Nathanson_exp[,clinical$Patient]
clinical$fustat<-clinical$OS.Event
clinical$futime<-clinical$OS
clinical$response[clinical$Response==0]<-"NDB"
clinical$response[clinical$Response==1]<-"DCB"
table(clinical$response)
rownames(clinical)<-clinical$Patient
write.table(clinical,"clinical_result.txt",sep="\t",quote=FALSE,row.names = F)

###每个样本都有生存数据
####和marker基因匹配上的基因名保留
load("G:/免疫亚型分型/MOVICS/Specific_markers.rda")
###marker基因表达谱
marker_exp<-Nathanson_exp[templates$probe,]
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
Nathanson_Response<-data.frame(Response=clinical$response,
                          Cluster=clinical$cluster)
fisher_res<-table(Nathanson_Response)
fisher_res
fisher_test<-fisher.test(fisher_res)
fisher_test$p.value
p<-ggplot(data = Nathanson_Response) + 
  ggtitle(paste("Fisher’s exact test P=",fisher_test$p.value,sep = ""))+
  geom_bar(mapping = aes(x = Cluster, fill = Response), position = "fill")+
  scale_fill_manual(values=c("#21498d","#e71d36"))+
  theme_bw()+theme(panel.grid=element_blank())
p
ggsave("Response.pdf",p,width = 5.5,height = 6)

