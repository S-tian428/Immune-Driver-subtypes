##免疫细胞浸润
#####Cibersort
rm(list=ls())
setwd("G:/免疫亚型分型/result8/Biomarker_ICB/Cibersort/")
#STAD.FPKM<-STAD_multiomics$STAD.FPKM
source('Cibersort.R')
###tpm
result <- CIBERSORT("LM22.txt","GSE91061_exp.txt", perm = 1000, QN = T)
save(result,file = "CIBERSORT_result.rda")
load("CIBERSORT_result.rda")
####整理成矩阵形式
cibersort_result<-as.data.frame(result[,-(23:25)])
Clinical<-read.table("G:/免疫亚型分型/result7/GSE91061/clinical_result.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
sample<-gsub("_.*","",rownames(cibersort_result))
index<-which(sample%in%Clinical$X)
cibersort_result<-cibersort_result[index,]
#####Pre
cibersort_pro<-cibersort_result[grep("On",rownames(cibersort_result)),]
rownames(cibersort_pro)<-gsub("_.*","",rownames(cibersort_pro))
cibersort_pro<-as.matrix(cibersort_pro)

GSE91061_exp<-read.table("GSE91061_exp.txt",
                         header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
sample<-gsub("_.*","",colnames(GSE91061_exp))
index<-which(sample%in%Clinical$X)
GSE91061_exp<-GSE91061_exp[,index]
GSE91061_pro<-GSE91061_exp[,grep("On",colnames(GSE91061_exp))]
colnames(GSE91061_pro)<-gsub("_.*","",colnames(GSE91061_pro))
CS_biomarker<-read.table("G:/免疫亚型分型/result8/Biomarker_ICB/Sig_Cox_result.txt",
                         header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
####加载CS2 biomarker
CS_biomarker<-rownames(CS_biomarker)
CS_pro_exp<-t(GSE91061_pro[rownames(GSE91061_pro)%in%CS_biomarker,])
library(ggcorrplot)
#install.packages("ggthemes")
library(ggthemes)
###治疗后
library(psych)
pro_cortest_psy <- corr.test(CS_pro_exp,cibersort_pro,method = "pearson")

cmt<-pro_cortest_psy$r
pmt<-pro_cortest_psy$p.adj
p<-ggcorrplot(cmt,method = "circle",outline.color = "white",
              ggtheme = theme_bw(),colors = c("#839EDB", "white", "#FF8D8D"),lab = T,lab_size=2,
              p.mat=pmt,insig="pch",pch.col = "red", pch.cex = 3, tl.cex = 12)
####单个细胞与基因相关性
pro_Imm_CS_exp<-data.frame(CS_pro_exp,cibersort_pro)
colnames(pro_Imm_CS_exp)
library(ggstatsplot)
p<-ggscatterstats(
  
  data = pro_Imm_CS_exp,
  
  x = IGHM,
  
  y = T.cells.CD8,
  
  type = "p",
  
  conf.level = 0.99,
  
  # marginal=F,
  
  messages = TRUE
  
)
p
ggsave("G:/免疫亚型分型/result8/Biomarker_ICB/免疫细胞相关性/IGHM_T.cells.CD8.pdf",
       width = 6,height = 6)
p<-ggscatterstats(
  
  data = pro_Imm_CS_exp,
  
  x = IGHM,
  
  y = T.cells.CD4.memory.activated,
  
  type = "p",
  
  conf.level = 0.99,
  
  # marginal=F,
  
  messages = TRUE
  
)
p
ggsave("G:/免疫亚型分型/result8/Biomarker_ICB/免疫细胞相关性/IGHM_T.cells.CD4.memory.activated.pdf",
       width = 6,height = 6)

p<-ggscatterstats(
  
  data = pro_Imm_CS_exp,
  
  x = IGHM,
  
  y = T.cells.CD8,
  
  type = "p",
  
  conf.level = 0.99,
  
  # marginal=F,
  
  messages = TRUE
  
)
p
ggsave("G:/免疫亚型分型/result8/Biomarker_ICB/免疫细胞相关性/IGHM_T.cells.CD8.pdf",
       width = 6,height = 6)
####CCL19
p<-ggscatterstats(
  
  data = pro_Imm_CS_exp,
  
  x = CCL19,
  
  y = T.cells.CD4.memory.activated,
  
  type = "p",
  
  conf.level = 0.99,
  
  # marginal=F,
  
  messages = TRUE
  
)
p
ggsave("G:/免疫亚型分型/result8/Biomarker_ICB/免疫细胞相关性/CCL19_T.cells.CD4.memory.activated.pdf",
       width = 6,height = 6)
p<-ggscatterstats(
  
  data = pro_Imm_CS_exp,
  
  x = CCL19,
  
  y = T.cells.CD8,
  
  type = "p",
  
  conf.level = 0.99,
  
  # marginal=F,
  
  messages = TRUE
  
)
p
ggsave("G:/免疫亚型分型/result8/Biomarker_ICB/免疫细胞相关性/CCL19_T.cells.CD8.pdf",
       width = 6,height = 6)
p<-ggscatterstats(
  
  data = pro_Imm_CS_exp,
  
  x = CCL19,
  
  y = B.cells.naive,
  
  type = "p",
  
  conf.level = 0.99,
  
  # marginal=F,
  
  messages = TRUE
  
)
p
ggsave("G:/免疫亚型分型/result8/Biomarker_ICB/免疫细胞相关性/CCL19_B.cells.naive.pdf",
       width = 6,height = 6)
