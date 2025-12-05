rm(list=ls())
#library(GSVA)
CS1_biomarker<-read.table("G:/免疫亚型分型/result5/CS1_biomarker.txt",
                          header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
####加载CS2 biomarker
CS2_biomarker<-read.table("G:/免疫亚型分型/result5/CS2_biomarker.txt",
                          header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
biomarker<-c(CS1_biomarker$Gene.Symbol,CS2_biomarker$Gene.Symbol)

Res_NoRes_gene<-read.table("G:/免疫亚型分型/result8/Biomarker_ICB/Response_NoResponse_gene.txt",
                          header=F,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
Res_NoRes_gene<-Res_NoRes_gene$V2

Pre_Pro_gene<-read.table("G:/免疫亚型分型/result8/Biomarker_ICB/Pre_Pro_gene.txt",
                         header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
Pre_Pro_gene<-Pre_Pro_gene$x
#####重要的CS_biomarker
Res_NoRes_biomarker<-intersect(Res_NoRes_gene,biomarker)
Pre_Pro_biomarker<-intersect(Pre_Pro_gene,biomarker)
CS_biomarker<-union(Res_NoRes_biomarker,Pre_Pro_biomarker)
write.table(CS_biomarker,"G:/免疫亚型分型/result8/Biomarker_ICB/ICB_CS_biomarker.txt",sep="\t",quote=FALSE,row.names = F)
Signature_gene<-union(Res_NoRes_gene,Pre_Pro_gene)


setwd("G:/免疫亚型分型/result8")
GSE91061_exp<-read.table("GSE91061_exp.txt",
                         header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
Clinical<-read.table("G:/免疫亚型分型/result7/GSE91061/clinical_result.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
sample<-gsub("_.*","",colnames(GSE91061_exp))
index<-which(sample%in%Clinical$X)
GSE91061_exp<-GSE91061_exp[,index]

#####Pre
GSE91061_pre<-GSE91061_exp[,grep("Pre",colnames(GSE91061_exp))]
colnames(GSE91061_pre)<-gsub("_.*","",colnames(GSE91061_pre))
signature_pre_exp<-t(GSE91061_pre[rownames(GSE91061_pre)%in%Signature_gene,])


GSE91061_pro<-GSE91061_exp[,grep("On",colnames(GSE91061_exp))]
colnames(GSE91061_pro)<-gsub("_.*","",colnames(GSE91061_pro))
signature_pro_exp<-t(GSE91061_pro[rownames(GSE91061_pro)%in%Signature_gene,])

biomarker_PreExp<-t(GSE91061_pre[rownames(GSE91061_pre)%in%CS_biomarker,])
biomarker_ProExp<-t(GSE91061_pro[rownames(GSE91061_pro)%in%CS_biomarker,])
####
######计算相关性
###治疗前
library(psych)
pre_cortest_psy <- corr.test(biomarker_PreExp,signature_pre_exp,method = "pearson")
pre_P<-as.data.frame(pre_cortest_psy$p)
pre_P$CS.biomarker<-rownames(pre_P)
###宽数据变长数据
library(tidyr)
pre_P_long <- gather(pre_P, key = "signature", value = "P",
                    -`CS.biomarker`)
pre_P_long<-na.omit(pre_P_long)
pre_P_sig<-pre_P_long[pre_P_long$P<0.05,]

write.table(pre_P_sig,"G:/免疫亚型分型/result8/Biomarker_ICB/pre_P_sig.txt",sep="\t",quote=FALSE,row.names = F)
pre_R<-as.data.frame(pre_cortest_psy$r)
pre_R$CS.biomarker<-rownames(pre_R)
###宽数据变长数据
pre_R_long <- gather(pre_R, key = "signature", value = "R",
                     -`CS.biomarker`)
pre_R_long<-na.omit(pre_R_long)
pre_R_sig<-pre_R_long[pre_R_long$R>=0.5,]
write.table(pre_R_sig,"G:/免疫亚型分型/result8/Biomarker_ICB/pre_R_sig.txt",sep="\t",quote=FALSE,row.names = F)
Pre_Corgene<-intersect(pre_P_sig$signature,pre_R_sig$signature)

###治疗后
library(psych)
pro_cortest_psy <- corr.test(biomarker_ProExp,signature_pro_exp,method = "pearson")
pro_P<-as.data.frame(pro_cortest_psy$p)
pro_P$CS.biomarker<-rownames(pro_P)
###宽数据变长数据
library(tidyr)
pro_P_long <- gather(pro_P, key = "signature", value = "P",
                     -`CS.biomarker`)
pro_P_long<-na.omit(pro_P_long)
pro_P_sig<-pro_P_long[pro_P_long$P<0.05,]
write.table(pro_P_sig,"G:/免疫亚型分型/result8/Biomarker_ICB/pro_P_sig.txt",sep="\t",quote=FALSE,row.names = F)
pro_R<-as.data.frame(pro_cortest_psy$r)
pro_R$CS.biomarker<-rownames(pro_R)
###宽数据变长数据
pro_R_long <- gather(pro_R, key = "signature", value = "R",
                     -`CS.biomarker`)
pro_R_long<-na.omit(pro_R_long)
pro_R_sig<-pro_R_long[pro_R_long$R>=0.5,]
write.table(pro_R_sig,"G:/免疫亚型分型/result8/Biomarker_ICB/pro_R_sig.txt",sep="\t",quote=FALSE,row.names = F)
pro_Corgene<-intersect(pro_P_sig$signature,pro_R_sig$signature)
Signiture_final<-intersect(Pre_Corgene,pro_Corgene)
write.table(Signiture_final,"G:/免疫亚型分型/result8/Biomarker_ICB/Signiture_final.txt",sep="\t",quote=FALSE,row.names = F)
