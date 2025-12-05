options(stringsAsFactors = FALSE)
###清空环境
rm(list=ls())
library(GEOquery)
#读取验证集数据GSE78220
setwd("G:/免疫亚型分型/result7/GSE78220")
GSE78220 = getGEO('GSE78220', destdir=".", AnnotGPL = F, getGPL = F) 
###临床数据
clinical<-pData(GSE78220[[1]])
GSE78220_exp<-read.table("GSE78220_PatientFPKM.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
rownames(GSE78220_exp)<-GSE78220_exp[,1]
GSE78220_exp<-GSE78220_exp[,-1]
GSE78220_exp<-log2(GSE78220_exp+1)
colnames(GSE78220_exp)<-gsub("\\..*","",colnames(GSE78220_exp))
index<-match(colnames(GSE78220_exp),clinical$title)
colnames(GSE78220_exp)<-clinical$geo_accession[index]
write.table(GSE78220_exp,"GSE78220_exp_result.txt",sep="\t",quote=FALSE)

####和marker基因匹配上的基因名保留
load("G:/免疫亚型分型/MOVICS/Specific_markers.rda")
###marker基因表达谱
marker_exp<-GSE78220_exp[templates$probe,]
marker_exp<-na.omit(marker_exp)
####共有42个maker基因的表达
#筛选出有生存的样本
GSE78220.marker.exp<-marker_exp[,colnames(marker_exp)%in%clinical$geo_accession]
save(GSE78220.marker.exp,file="GSE78220.marker.exp.rda")#192个样本
Clinical_outcome<-clinical[clinical$geo_accession%in%colnames(GSE78220.marker.exp),]

#load("H:/免疫亚型分型/MOVICS/STAD.multiomics.rda")
#STAD.surv<-STAD_multiomics$surv.info
GSE78220_clin.info<-data.frame(Response=Clinical_outcome$characteristics_ch1.1,
                               futime=Clinical_outcome$characteristics_ch1.6,
                               fustat=Clinical_outcome$`vital status:ch1`)
rownames(GSE78220_clin.info)<-Clinical_outcome$geo_accession
GSE78220_clin.info$Response
GSE78220_clin.info$Response<-gsub("anti-pd-1 response: ","",GSE78220_clin.info$Response)
table(GSE78220_clin.info$Response)
GSE78220_clin.info$futime<-gsub("overall survival \\(days): ","",GSE78220_clin.info$futime)
GSE78220_clin.info<-GSE78220_clin.info[-7,]
GSE78220_clin.info$futime<-as.numeric(GSE78220_clin.info$futime)
GSE78220_clin.info$fustat[GSE78220_clin.info$fusta%in%"Dead"]=1
GSE78220_clin.info$fustat[GSE78220_clin.info$fusta%in%"Alive"]=0
GSE78220_clin.info$fustat<-as.numeric(GSE78220_clin.info$fustat)
#NTP预测
library(MOVICS)
getwd()
yau.ntp.pred <- runNTP(expr       = GSE78220.marker.exp,
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
str(GSE78220_clin.info)

surv.yau <- compSurv(moic.res         = yau.ntp.pred,
                     surv.info        = GSE78220_clin.info,
                     convt.time       = "d", # switch to month
                     surv.median.line = "hv", # switch to both
                     fig.name         = "KAPLAN-MEIER CURVE OF NTP FOR YAU") 

save(yau.ntp.pred, file ="yau.ntp.pred.Rdata")
table(yau.ntp.pred$clust.res)
index<-match(rownames(GSE78220_clin.info),yau.ntp.pred$clust.res$samID)
GSE78220_clin.info$cluster<-yau.ntp.pred$clust.res$clust[index]
GSE78220_clin.info$RES<-"NDB"
GSE78220_clin.info$RES[GSE78220_clin.info$Response%in%"Complete Response"|
                         GSE78220_clin.info$Response%in%"Partial Response"]<-"DCB"
save(GSE78220_clin.info,file = "GSE78220_clin.info.rda")


#####免疫治疗反应比例分析
library(ggplot2)
GSE78220_Response<-data.frame(Response=GSE78220_clin.info$RES,
                               Cluster=GSE78220_clin.info$cluster)
fisher_res<-table(GSE78220_Response)
fisher_test<-fisher.test(fisher_res)
fisher_test$p.value
p<-ggplot(data = GSE78220_Response) + 
  ggtitle(paste("Fisher’s exact test P=",fisher_test$p.value,sep = ""))+
  geom_bar(mapping = aes(x = Cluster, fill = Response ), position = "fill")+
  scale_fill_manual(values=c("#f29c71","#21498d"))+
  theme_bw()+theme(panel.grid=element_blank())
p
ggsave("Response.pdf",p,width = 5.5,height = 6)

