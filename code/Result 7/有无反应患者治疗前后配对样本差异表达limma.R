options(stringsAsFactors = FALSE)
###清空环境
rm(list=ls())
library(limma)
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
write.table(GSE91061_exp,"G:/免疫亚型分型/result8/GSE91061_exp.txt",sep="\t",quote=FALSE)
###临床数据
Clinical<-read.table("G:/免疫亚型分型/result7/GSE91061/clinical_result.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
####Presonse  Pre vs Pro
Presonse_patients<-Clinical$X[Clinical$response%in%"DCB"]
sample<-gsub("_.*","",colnames(GSE91061_exp))
index<-which(sample%in%Presonse_patients)
Presonse_exp<-GSE91061_exp[,index]
###配对样本差异表达分析
Presonse_paired<-data.frame(sample=colnames(Presonse_exp))
Presonse_paired$Group<-"Pro"
Presonse_paired$Group[grep("Pre",colnames(Presonse_exp))]<-"Pre"
Presonse_paired$subject<-gsub("_.*","",Presonse_paired$sample)
individuals<-factor(Presonse_paired$subject)
treatment<-factor(Presonse_paired$Group)
design_paried <- model.matrix(~ individuals + treatment)
fit2 <- lmFit(Presonse_exp,design_paried)
fit2 <- eBayes(fit2)
tempOutput<-topTable(fit2, adjust="BH", coef=1, n=Inf)
####biomarker
####加载CS1 biomarker
CS1_biomarker<-read.table("G:/免疫亚型分型/result5/CS1_biomarker.txt",
                          header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
####加载CS2 biomarker
CS2_biomarker<-read.table("G:/免疫亚型分型/result5/CS2_biomarker.txt",
                          header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
biomarkers<-c(CS1_biomarker$Gene.Symbol,CS2_biomarker$Gene.Symbol)
biomarker_Response_PreOn<-tempOutput[biomarkers,]
biomarker_Response_PreOn<-na.omit(biomarker_Response_PreOn)
write.table(biomarker_Response_PreOn,"G:/免疫亚型分型/result8/Biomarker_ICB/Response_PreOn_limma.txt",sep="\t",quote=FALSE)


####NOPresonse  Pre vs Pro
NOPresonse_patients<-Clinical$X[Clinical$response%in%"NDB"]
index<-which(sample%in%NOPresonse_patients)
NOPresonse_exp<-GSE91061_exp[,index]
###配对样本差异表达分析
NOPresonse_paired<-data.frame(sample=colnames(NOPresonse_exp))
NOPresonse_paired$Group<-"Pro"
NOPresonse_paired$Group[grep("Pre",colnames(NOPresonse_exp))]<-"Pre"
NOPresonse_paired$subject<-gsub("_.*","",NOPresonse_paired$sample)
individuals<-factor(NOPresonse_paired$subject)
treatment<-factor(NOPresonse_paired$Group)
design_paried <- model.matrix(~ individuals + treatment)
fit2 <- lmFit(NOPresonse_exp,design_paried)
fit2 <- eBayes(fit2)
tempOutput<-topTable(fit2, adjust="BH", coef=1, n=Inf)
####biomarker
biomarker_NOPresonse_PreOn<-tempOutput[biomarkers,]
biomarker_NOPresonse_PreOn<-na.omit(biomarker_NOPresonse_PreOn)
write.table(biomarker_NOPresonse_PreOn,"G:/免疫亚型分型/result8/Biomarker_ICB/NOResponse_PreOn_limma.txt",sep="\t",quote=FALSE)
#####
gene1<-biomarker_NOPresonse_PreOn[biomarker_NOPresonse_PreOn$adj.P.Val>=0.05,]
gene2<-biomarker_Response_PreOn[biomarker_Response_PreOn$adj.P.Val<0.05,]
signature_gene<-intersect(rownames(gene1),rownames(gene2))
write.table(signature_gene,"G:/免疫亚型分型/result8/Biomarker_ICB/NORes_Response_signature.txt",sep="\t",quote=FALSE,col.names = F)
