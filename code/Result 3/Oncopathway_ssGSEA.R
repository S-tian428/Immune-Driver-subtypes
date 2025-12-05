rm(list=ls())
library(GSVA)
setwd("G:/免疫亚型分型/result3/oncopathway")
gene_set<-read.table("oncogenic pathways.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
list<- split(as.matrix(gene_set)[,1], gene_set[,2])
####表达
load("G:/免疫亚型分型/MOVICS/STAD.multiomics.rda")
Expression<-as.matrix(STAD_multiomics$STAD.FPKM)

##We can calculate GSVA enrichment scores as follows:
library(GSVA)
#expr--Expression
#cellMarker--list
gsvaPar <- ssgseaParam(exprData = Expression, 
                       geneSets = list,
                       normalize = TRUE)
gsva_data <- gsva(gsvaPar, verbose = FALSE)
gsva_matrix<-data.frame(t(gsva_data))
###cluster
load("G:/免疫亚型分型/MOVICS/cmoic.brca.rda")
Sample_cluster<-cmoic.brca$clust.res
Sample_cluster$clust<-paste0("CS",Sample_cluster$clust)
index<-match(rownames(gsva_matrix),Sample_cluster$samID)
gsva_matrix$Cluster<-Sample_cluster$clust[index]
library(ggplot2)
library(ggpubr)
library(forcats)
require(ggsci)
i=1
gsva_matrix$Cluster<-as.factor(gsva_matrix$Cluster)
for (i in 1:10) {
  gene<-colnames(gsva_matrix)[i]
  p<-ggplot(gsva_matrix, aes(x=Cluster, y=gsva_matrix[,i],color=Cluster))+geom_violin(linewidth=1)+
    scale_color_brewer(palette="Dark2")
  p<-p+stat_summary(fun=mean, geom="point", size=5,color="darkgrey")+
    theme_classic(base_size = 15)+ylab("")+xlab("")+ggtitle(gene)
  p<-p+stat_compare_means(method = "wilcox.test")
  p
  myfilename <-paste(gene,"_boxplot.pdf",sep="")
  #library(ggplot2)
  ggsave(myfilename,p,width=4,height = 4)
  
  print(i)
}






