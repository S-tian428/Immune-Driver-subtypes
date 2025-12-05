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
#######CS1
CS1_biomarker<-read.table("G:/免疫亚型分型/result5/CS1_biomarker.txt",
                          header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
#######CS2
CS2_biomarker<-read.table("G:/免疫亚型分型/result5/CS2_biomarker.txt",
                          header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
Signiture<-c(CS1_biomarker$Gene.Symbol,CS2_biomarker$Gene.Symbol)
Signiture_limma<-tempOutput[Signiture,]
Signiture_limma<-na.omit(Signiture_limma)
write.table(Signiture_limma,"G:/免疫亚型分型/result8/Biomarker_ICB/On_NOResRse_limma.txt",sep="\t",quote=FALSE)


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
Signiture_limma1<-tempOutput[Signiture,]
Signiture_limma1<-na.omit(Signiture_limma1)
write.table(Signiture_limma,"G:/免疫亚型分型/result8/Biomarker_ICB/Pre_NOResRse_limma.txt",sep="\t",quote=FALSE)
#####
gene_On<-Signiture_limma[Signiture_limma$P.Value<0.05,]
gene_Pre<-Signiture_limma1[Signiture_limma1$P.Value>=0.05,]
signature_gene<-intersect(rownames(gene_On),rownames(gene_Pre))
write.table(signature_gene,"G:/免疫亚型分型/result8/Biomarker_ICB/Pre_On_signature.txt",sep="\t",quote=FALSE,col.names = F)
