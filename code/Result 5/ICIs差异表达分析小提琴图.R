##免疫细胞浸润
#####Cibersort
rm(list=ls())
#library(GEOquery)
setwd("H:/免疫亚型分型/result6/GSE62254")
load("yau.ntp.pred.Rdata")
load("GSE62254_NormExp.rda")
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
ICIs_limma<-na.omit(ICIs_limma)
Sig_genes<-rownames(ICIs_limma)
###绘制小提琴图
Sig_exp<-gene_exp[Sig_genes,]
Sig_exp<-as.data.frame(t(Sig_exp))
Sig_exp$Cluster<-clust.res$clust
library(ggplot2)
library(ggpubr)
i=1
Sig_exp$Cluster<-as.factor(Sig_exp$Cluster)
setwd("H:/免疫亚型分型/result6/GSE62254/ICIs/")
for (i in 1:length(Sig_genes)) {
  gene<-colnames(Sig_exp)[i]
  p<-ggplot(Sig_exp, aes(x=Cluster, y=Sig_exp[,i],color=Cluster))+geom_violin(linewidth=1)+
    scale_color_brewer(palette="Dark2")
  p
  p<-p+stat_summary(fun=mean, geom="point", size=5,color="darkgrey")+
    theme_bw()+ylab(gene)+xlab("")
  p<-p+stat_compare_means(method = "wilcox.test")
  p
  myfilename <-paste0(gene,".pdf")
  #library(ggplot2)
  ggsave(myfilename,p,width=3.8,height = 3.5)
  
  print(i)
}
