options(stringsAsFactors = FALSE)
###清空环境
rm(list=ls())
library(limma)
setwd("G:/免疫亚型分型/result8")
GSE91061_exp<-read.table("GSE91061_exp.txt",
                         header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
GSE91061_pro<-GSE91061_exp[,grep("On",colnames(GSE91061_exp))]
colnames(GSE91061_pro)<-gsub("_.*","",colnames(GSE91061_pro))

Clinical<-read.table("G:/免疫亚型分型/result7/GSE91061/clinical_result.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)

res_sample<-Clinical$X[Clinical$response%in%"DCB"]
non_sample<-Clinical$X[Clinical$response%in%"NDB"]
#####
expression<-GSE91061_pro[,c(non_sample,res_sample)]
library(limma)
#####获得肿瘤样本和正常样本标筿
group1 <- rep('NoResponse',length(non_sample))
group2 <- rep('Response', length(res_sample))
#######热图中样本的标签及差异表达分析中的标筿
grouplist<-c(group1,group2)
#######差异表达分析
design <- model.matrix(~0+factor(grouplist))
colnames(design)=levels(factor(grouplist))#改design列名为分组信恿
rownames(design)=colnames(expression)
contrast.matrix <- makeContrasts("Response-NoResponse", levels=design)
fit <- lmFit(expression,design)
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效枿
fit2 <- eBayes(fit2)
tempOutput<-topTable(fit2, adjust="BH", coef=1, n=Inf)
#write.table(tempOutput,"Response_limma.txt",sep="\t",quote=FALSE)
pro_sigLimma<-tempOutput[
                           tempOutput$P.Value<0.05,]
####pre
GSE91061_pre<-GSE91061_exp[,grep("Pre",colnames(GSE91061_exp))]
colnames(GSE91061_pre)<-gsub("_.*","",colnames(GSE91061_pre))
expression<-GSE91061_pre[,c(non_sample,res_sample)]
#####获得肿瘤样本和正常样本标筿
#######热图中样本的标签及差异表达分析中的标筿
#######差异表达分析
design <- model.matrix(~0+factor(grouplist))
colnames(design)=levels(factor(grouplist))#改design列名为分组信恿
rownames(design)=colnames(expression)
contrast.matrix <- makeContrasts("Response-NoResponse", levels=design)
fit <- lmFit(expression,design)
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效枿
fit2 <- eBayes(fit2)
tempOutput<-topTable(fit2, adjust="BH", coef=1, n=Inf)
pre_nosig_Limma<-tempOutput[tempOutput$P.Value>=0.05,]

#####
pre_Pro_gene<-intersect(rownames(pre_nosig_Limma),rownames(pro_sigLimma))
write.table(pre_Pro_gene,"G:/免疫亚型分型/result8/Biomarker_ICB/Pre_Pro_gene.txt",sep="\t",quote=FALSE,row.names = F)
###基因limma结果
CS_biomarker<-read.table("G:/免疫亚型分型/result8/Biomarker_ICB/Sig_Cox_result.txt",
                         header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
CS_biomarker<-rownames(CS_biomarker)
pre_nosig_Limma[CS_biomarker,]
pro_sigLimma[CS_biomarker,]
