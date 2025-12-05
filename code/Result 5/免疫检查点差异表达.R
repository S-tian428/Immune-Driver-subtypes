##免疫细胞浸润
#####Cibersort
rm(list=ls())
setwd("I:/免疫亚型分型/result6/GSE26253")
load("GSE26253_NormExp.rda")
load("yau.ntp.pred.Rdata")
clust.res<-yau.ntp.pred$clust.res
clust.res$clust<-paste0("CS",clust.res$clust)
clust.res<-clust.res[order(clust.res$clust),]

gene_exp<-gene_exp[,clust.res$samID]
library(limma)
group1 <- clust.res$clust[clust.res$clust%in%"CS1"]
group2 <- clust.res$clust[clust.res$clust%in%"CS2"]
grouplist<-c(group1,group2)
#######差异表达分析
design <- model.matrix(~0+factor(grouplist))
colnames(design)=levels(factor(grouplist))#改design列名为分组信息
rownames(design)=colnames(gene_exp)
contrast.matrix <- makeContrasts("CS2-CS1", levels=design)
fit <- lmFit(gene_exp,design)
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
fit2 <- eBayes(fit2)
tempOutput<-topTable(fit2, adjust="BH", coef=1, n=Inf)
limma_result<-na.omit(tempOutput)
Sig_limma<-limma_result[limma_result$adj.P.Val<0.05,]
Sig_limma<-Sig_limma[order(Sig_limma$logFC),]
ICIs<-c("HLA-DRB5","HLA-DPB1","HLA-DOA","HLA-DPA1","HLA-DQA1","CD27")
ICIs_limma<-Sig_limma[ICIs,]
