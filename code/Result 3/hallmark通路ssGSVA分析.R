rm(list=ls())
#install.packages("tidyverse")
library(data.table)###读文件
library(dplyr)#####去重复
library(tidyverse)####列变行名
library(GSVA)
library(GSEABase)
#install.packages('R.utils')
load("H:/免疫亚型分型/MOVICS/STAD.multiomics.rda")
STAD_fpkm<-STAD_multiomics$STAD.FPKM
c2gmt<-clusterProfiler::read.gmt("H:/免疫亚型分型/result3/hallmark.pathway/h.all.v2023.2.Hs.symbols.gmt")
genesets = split(c2gmt$gene, c2gmt$term)
TCGA_exp<-as.matrix(STAD_fpkm)
gs.exp <- gsva(TCGA_exp,genesets)
write.table(gs.exp,"H:/免疫亚型分型/result3/hallmark.pathway/hallmark_result.txt",sep="\t",quote=FALSE)






