rm(list=ls())
signatureGSVA<-read.table("G:/免疫亚型分型/result3/GSVA/活性相关性/Immune signature_GSVA.txt",
                             header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)

###cluster
load("G:/免疫亚型分型/MOVICS/cmoic.brca.rda")
Sample_cluster<-cmoic.brca$clust.res
Sample_cluster$clust<-paste0("CS",Sample_cluster$clust)
CS1_sample<-Sample_cluster$samID[Sample_cluster$clust%in%"CS1"]
CS2_sample<-Sample_cluster$samID[Sample_cluster$clust%in%"CS2"]
CS1_GSVA<-signatureGSVA[CS1_sample,]
CS2_GSVA<-signatureGSVA[CS2_sample,]
####相关性
library(ggcorrplot)
library(ggthemes)
CS1_corr<-round(cor(as.matrix(CS1_GSVA)),3)
# 计算相关性的P值矩阵
CS1_p.mat <- cor_pmat(as.matrix(CS1_GSVA))
p1<-ggcorrplot(CS1_corr, method = "square", type = "upper", ggtheme = ggplot2::theme_bw(), title = "", 
           show.legend = TRUE, legend.title = "Corr", show.diag = T, 
           colors = c("#b2df8a", "white", "#0072b5"), outline.color = "white", 
           hc.order = F,  lab = T, lab_col = "black", 
           lab_size = 2.5, p.mat = CS1_p.mat, sig.level = 0.05, insig = "blank", pch = 4, pch.col = "black", pch.cex = 5, tl.cex = 12, 
           tl.col = "black", tl.srt = 45, digits = 2)
p1
ggsave("G:/免疫亚型分型/result3/GSVA/活性相关性/CS1_signature_cor.pdf",p1,height=10,width=10)
CS2_corr<-round(cor(as.matrix(CS2_GSVA)),3)
# 计算相关性的P值矩阵
CS2_p.mat <- cor_pmat(as.matrix(CS2_GSVA))
p2<-ggcorrplot(CS2_corr, method = "square", type = "lower", ggtheme = ggplot2::theme_bw(), title = "", 
           show.legend = TRUE, legend.title = "Corr", show.diag = T, 
           colors = c("#ffff99", "white", "#f79c4a"), outline.color = "white", 
            hc.method = "complete", lab = T, lab_col = "black", 
           lab_size = 2.5, p.mat = CS2_p.mat, sig.level = 0.05, insig = "blank", pch = 4, pch.col = "black", pch.cex = 5, tl.cex = 12, 
           tl.col = "black", tl.srt = 45, digits = 2)
p2
ggsave("G:/免疫亚型分型/result3/GSVA/活性相关性/CS2_signature_cor.pdf",p2,height=10,width=10)











