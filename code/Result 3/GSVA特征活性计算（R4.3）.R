rm(list=ls())
library(GSVA)
setwd("G:/免疫亚型分型/result3/oncopathway")
gmt1<-clusterProfiler::read.gmt("G:/免疫亚型分型/result3/GSVA/Immune signature.gmt")
index<-which(gmt1$gene%in%"")
gmt1<-gmt1[-index,]
list<-split(gmt1$gene, gmt1$term)
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
write.table(gsva_matrix,"G:/免疫亚型分型/result3/GSVA/Immune signature_GSVA.txt",
            sep="\t",quote=FALSE)
##########
gmt2<-clusterProfiler::read.gmt("G:/免疫亚型分型/result3/GSVA/活性相关性/Immune signature1.gmt")
index<-which(gmt2$gene%in%"")
gmt2<-gmt2[-index,]
list<-split(gmt2$gene, gmt2$term)

#expr--Expression
#cellMarker--list
gsvaPar <- ssgseaParam(exprData = Expression, 
                       geneSets = list,
                       normalize = TRUE)
gsva_data <- gsva(gsvaPar, verbose = FALSE)
gsva_matrix<-data.frame(t(gsva_data))
write.table(gsva_matrix,"G:/免疫亚型分型/result3/GSVA/活性相关性/Immune signature1_GSVA.txt",
            sep="\t",quote=FALSE)
