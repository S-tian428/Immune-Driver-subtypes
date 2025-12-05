options(stringsAsFactors = FALSE)
###清空环境
rm(list=ls())
setwd("G:/免疫亚型分型/result8")
GSE91061_exp<-read.table("GSE91061_exp.txt",
                         header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
Clinical<-read.table("G:/免疫亚型分型/result7/GSE91061/clinical_result.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
sample<-gsub("_.*","",colnames(GSE91061_exp))
index<-which(sample%in%Clinical$X)
GSE91061_exp<-GSE91061_exp[,index]
GSE91061_exp_On<-GSE91061_exp[,grep("On",colnames(GSE91061_exp))]
colnames(GSE91061_exp_On)<-gsub("_.*","",colnames(GSE91061_exp_On))
Signature_gene<-read.table("G:/免疫亚型分型/result8/Biomarker_ICB/ICB_CS_biomarker.txt",
                           header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
Signature_gene<-Signature_gene$x
Signature_EXP<-GSE91061_exp_On[Signature_gene,]
Signature_EXP<-data.frame(t(Signature_EXP))
index<-match(Clinical$X,rownames(Signature_EXP))
clinical_exp<-cbind(Clinical,Signature_EXP[index,])
#####但cox回归分析
library(stringr)
library(caret)
library(survminer)
library(survival)
i=17
colnames(clinical_exp)[1:10]
Cox_result<-NULL
for (i in 17:ncol(clinical_exp)) {
  cox_genes <- coxph(Surv(OS, OS.Event) ~ clinical_exp[,i], data = clinical_exp)
  coef <- coef(cox_genes) #回归系数
  SE <- sqrt(diag(vcov(cox_genes))) #标准误
  HR <- exp(coef) #风险比
  cox_need <- cbind(HR = HR,
                    LowerCI = exp(coef - qnorm(.975, 0, 1) * SE),
                    upperCI = exp(coef + qnorm(.975, 0, 1) * SE),
                    pvalue = 1 - pchisq((coef/SE)^2, 1))
  rownames(cox_need)<-colnames(clinical_exp)[i]
  
  Cox_result<-rbind(Cox_result,cox_need)
}
Cox_result<-as.data.frame(Cox_result)
Sig_Cox_result<-Cox_result[Cox_result$pvalue<0.05,]
write.table(Sig_Cox_result,"G:/免疫亚型分型/result8/Biomarker_ICB/Sig_Cox_result.txt",sep="\t",quote=FALSE,row.names = T)

