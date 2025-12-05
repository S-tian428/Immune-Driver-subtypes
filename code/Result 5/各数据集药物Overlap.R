rm(list=ls())
setwd("G:/免疫亚型分型/result6/Drug_Venn/")
GSE26253<-read.table("G:/免疫亚型分型/result6/GSE26253/Oncopredict/Wilcox_test.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
GSE26253_CS1_drug<-GSE26253$Drug[GSE26253$response_cluster%in%"CS1"]
GSE26253_CS2_drug<-GSE26253$Drug[GSE26253$response_cluster%in%"CS2"]

write.table(GSE26253_CS1_drug,"CS1/GSE26253.txt",
            sep="\t",quote=FALSE,row.names=F,col.names = F)
write.table(GSE26253_CS2_drug,"CS2/GSE26253.txt",
            sep="\t",quote=FALSE,row.names=F,col.names = F)
GSE15459<-read.table("G:/免疫亚型分型/result6/GSE15459/Oncopredict/Wilcox_test.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
GSE15459_CS1_drug<-GSE15459$Drug[GSE15459$response_cluster%in%"CS1"]
GSE15459_CS2_drug<-GSE15459$Drug[GSE15459$response_cluster%in%"CS2"]
write.table(GSE15459_CS1_drug,"CS1/GSE15459.txt",
            sep="\t",quote=FALSE,row.names=F,col.names = F)
write.table(GSE15459_CS2_drug,"CS2/GSE15459.txt",
            sep="\t",quote=FALSE,row.names=F,col.names = F)

GSE84437<-read.table("G:/免疫亚型分型/result6/GSE84437/Oncopredict/Wilcox_test.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
GSE84437_CS1_drug<-GSE84437$Drug[GSE84437$response_cluster%in%"CS1"]
GSE84437_CS2_drug<-GSE84437$Drug[GSE84437$response_cluster%in%"CS2"]
write.table(GSE84437_CS1_drug,"CS1/GSE84437.txt",
            sep="\t",quote=FALSE,row.names=F,col.names = F)
write.table(GSE84437_CS2_drug,"CS2/GSE84437.txt",
            sep="\t",quote=FALSE,row.names=F,col.names = F)
TCGA<-read.table("G:/免疫亚型分型/result5/Oncopredict/Wilcox_test.txt",
                 header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
TCGA_CS1_drug<-TCGA$Drug[TCGA$response_cluster%in%"CS1"]
TCGA_CS2_drug<-TCGA$Drug[TCGA$response_cluster%in%"CS2"]
write.table(TCGA_CS1_drug,"CS1/TCGA.txt",
            sep="\t",quote=FALSE,row.names=F,col.names = F)
write.table(TCGA_CS2_drug,"CS2/TCGA.txt",
            sep="\t",quote=FALSE,row.names=F,col.names = F)
CS1_drug <- Reduce(intersect, list(
            GSE26253_CS1_drug,GSE84437_CS1_drug,
            GSE15459_CS1_drug,TCGA_CS1_drug))
write.table(CS1_drug,"CS1_intersect_drug.txt",
            sep="\t",quote=FALSE,row.names=F,col.names = F)
CS2_drug <- Reduce(intersect, list(
  GSE26253_CS2_drug,GSE84437_CS2_drug,
  GSE15459_CS2_drug,TCGA_CS2_drug))
write.table(CS2_drug,"CS2_intersect_drug.txt",
            sep="\t",quote=FALSE,row.names=F,col.names = F)
#####CS1和CS2 IC50差值的平均值排序
CS1_DIFF<-data.frame(TCGA=TCGA$Difference[TCGA$Drug%in%CS1_drug],
                     GSE15459=GSE15459$Difference[GSE15459$Drug%in%CS1_drug],
                     GSE26253=GSE26253$Difference[GSE26253$Drug%in%CS1_drug],
                     GSE84437=GSE84437$Difference[GSE84437$Drug%in%CS1_drug])
rownames(CS1_DIFF)<-CS1_drug 
a<-apply(CS1_DIFF, 1, mean)
CS1_DIFF<-cbind(CS1_DIFF,mean_diff=a)
CS1_DIFF<-CS1_DIFF[order(CS1_DIFF$mean_diff,decreasing = T),]
write.table(CS1_DIFF,"CS1_DIFF.txt",
            sep="\t",quote=FALSE,row.names=T,col.names = T)


#####CS1和CS2 IC50差值的平均值排序
CS2_DIFF<-data.frame(TCGA=abs(TCGA$Difference[TCGA$Drug%in%CS2_drug]),
                     GSE15459=abs(GSE15459$Difference[GSE15459$Drug%in%CS2_drug]),
                     GSE26253=abs(GSE26253$Difference[GSE26253$Drug%in%CS2_drug]),
                     GSE84437=abs(GSE84437$Difference[GSE84437$Drug%in%CS2_drug]))
rownames(CS2_DIFF)<-CS2_drug 
b<-apply(CS2_DIFF, 1, mean)
CS2_DIFF<-cbind(CS2_DIFF,mean_diff=b)
CS2_DIFF<-CS2_DIFF[order(CS2_DIFF$mean_diff,decreasing = T),]
write.table(CS2_DIFF,"CS2_DIFF.txt",
            sep="\t",quote=FALSE,row.names=T,col.names = T)
