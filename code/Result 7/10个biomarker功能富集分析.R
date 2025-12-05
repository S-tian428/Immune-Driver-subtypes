rm(list=ls())
CS_biomarker<-read.table("G:/免疫亚型分型/result8/Biomarker_ICB/Sig_Cox_result.txt",
                         header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
####加载CS2 biomarker
CS_biomarker<-rownames(CS_biomarker)
library(clusterProfiler)
library(org.Hs.eg.db)
CS_ids = bitr(CS_biomarker, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#sce.markers = merge(sce.markers, ids, by.x='V1', by.y='SYMBOL')
#####GO分析 
CS_go<-enrichGO(CS_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.05,keyType = 'ENTREZID')
#使用simplify 对GO富集分析结果进行精简
CS_go<- simplify(CS_go)
CS_result<-CS_go@result
CS_KEG<- enrichKEGG(gene=CS_ids$ENTREZID, organism='hsa', pvalueCutoff=0.05)
CS_kegg<-CS_KEG@result
