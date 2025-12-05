rm(list=ls())
#library(GSVA)
CS_biomarker<-read.table("G:/免疫亚型分型/result8/Biomarker_ICB/Sig_Cox_result.txt",
                          header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
####加载CS2 biomarker
CS_biomarker<-rownames(CS_biomarker)

ICIs<-c("PDCD1","CTLA4")
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
ICI_pre_exp<-t(GSE91061_pre[rownames(GSE91061_pre)%in%ICIs,])
CS_pre_exp<-t(GSE91061_pre[rownames(GSE91061_pre)%in%CS_biomarker,])


GSE91061_pro<-GSE91061_exp[,grep("On",colnames(GSE91061_exp))]
colnames(GSE91061_pro)<-gsub("_.*","",colnames(GSE91061_pro))
ICI_pro_exp<-t(GSE91061_pro[rownames(GSE91061_pro)%in%ICIs,])
CS_pro_exp<-t(GSE91061_pro[rownames(GSE91061_pro)%in%CS_biomarker,])

####
######计算相关性
###治疗前
library(psych)
pre_cortest_psy <- corr.test(ICI_pre_exp,CS_pre_exp,method = "pearson")
pre_P<-as.data.frame(pre_cortest_psy$p)
pre_P$ICIs<-rownames(pre_P)
###宽数据变长数据
library(tidyr)
pre_P_long <- gather(pre_P, key = "signature", value = "P",
                    -`ICIs`)
pre_P_long<-na.omit(pre_P_long)
pre_P_sig<-pre_P_long[pre_P_long$P<0.05,]
#write.table(pre_P_sig,"G:/免疫亚型分型/result8/Biomarker_ICB/pre_P_sig.txt",sep="\t",quote=FALSE,row.names = F)
pre_R<-as.data.frame(pre_cortest_psy$r)
pre_R$ICIs<-rownames(pre_R)
###宽数据变长数据
pre_R_long <- gather(pre_R, key = "signature", value = "R",
                     -`ICIs`)
pre_R_long<-na.omit(pre_R_long)
pre_R_sig<-pre_R_long[pre_R_long$R>=0.5,]
#write.table(pre_R_sig,"G:/免疫亚型分型/result8/Biomarker_ICB/pre_R_sig.txt",sep="\t",quote=FALSE,row.names = F)

###治疗后
library(psych)
pro_cortest_psy <- corr.test(ICI_pro_exp,CS_pro_exp,method = "pearson")
pro_P<-as.data.frame(pro_cortest_psy$p)
pro_P$ICIs<-rownames(pro_P)
###宽数据变长数据
library(tidyr)
pro_P_long <- gather(pro_P, key = "signature", value = "P",
                     -`ICIs`)
pro_P_long<-na.omit(pro_P_long)
pro_P_sig<-pro_P_long[pro_P_long$P<0.05,]
#write.table(pro_P_sig,"G:/免疫亚型分型/result8/Biomarker_ICB/pro_P_sig.txt",sep="\t",quote=FALSE,row.names = F)
pro_R<-as.data.frame(pro_cortest_psy$r)
pro_R$ICIs<-rownames(pro_R)
###宽数据变长数据
pro_R_long <- gather(pro_R, key = "signature", value = "R",
                     -`ICIs`)
pro_R_long<-na.omit(pro_R_long)
pro_R_sig<-pro_R_long[pro_R_long$R>=0.3,]
#####表达相关性图
library(ggstatsplot)
pre_ICIs_CS_exp<-data.frame(ICI_pre_exp,CS_pre_exp)
colnames(pre_ICIs_CS_exp)
p<-ggscatterstats(
  
  data = pre_ICIs_CS_exp,
  
  x = IGHM,
  
  y = CTLA4,
  
  type = "p",
  
  conf.level = 0.99,
  
  # marginal=F,
  
  messages = TRUE
  
)
ggsave("G:/免疫亚型分型/result8/Biomarker_ICB/表达相关性/ICIs/CTLA4_IGHM_Pre.pdf",
       width = 6,height = 6)

p<-ggscatterstats(
  
  data = pre_ICIs_CS_exp,
  
  x = CCL19,
  
  y = CTLA4,
  
  type = "p",
  
  conf.level = 0.99,
  
  # marginal=F,
  
  messages = TRUE
  
)
ggsave("G:/免疫亚型分型/result8/Biomarker_ICB/表达相关性/ICIs/CTLA4_CCL19_Pre.pdf",
       width = 6,height = 6)


p<-ggscatterstats(
  
  data = pre_ICIs_CS_exp,
  
  x = MFAP4,
  
  y = PDCD1,
  
  type = "p",
  
  conf.level = 0.99,
  
  # marginal=F,
  
  messages = TRUE
  
)
p
ggsave("G:/免疫亚型分型/result8/Biomarker_ICB/表达相关性/ICIs/MFAP4_PDCD1_Pre.pdf",
       width = 6,height = 6)

library(ggcorrplot)
library(ggthemes)
library(psych)

cmt<-pro_cortest_psy$r
pmt<-pro_cortest_psy$p.adj
p<-ggcorrplot(cmt,method = "circle",outline.color = "white",
           ggtheme = theme_bw(),colors = c("#839EDB", "white", "#FF8D8D"),lab = T,lab_size=2,
           p.mat=pmt,insig="pch",pch.col = "red", pch.cex = 3, tl.cex = 12)
ggsave("G:/免疫亚型分型/result8/Biomarker_ICB/表达相关性/ICIs/corrplot_Pro.pdf",
       width = 8,height = 6)

cmt<-pre_cortest_psy$r
pmt<-pre_cortest_psy$p.adj
p<-ggcorrplot(cmt,method = "circle",outline.color = "white",
           ggtheme = theme_bw(),colors = c("#839EDB", "white", "#FF8D8D"),lab = T,lab_size=2,
           p.mat=pmt,insig="pch",pch.col = "red", pch.cex = 3, tl.cex = 12)
ggsave("G:/免疫亚型分型/result8/Biomarker_ICB/表达相关性/ICIs/corrplot_Pre.pdf",
       width = 8,height = 6)



pro_ICIs_CS_exp<-data.frame(ICI_pro_exp,CS_pro_exp)
p<-ggscatterstats(
  
  data = pro_ICIs_CS_exp,
  
  x = IGHM,
  
  y = CTLA4,
  
  type = "p",
  
  conf.level = 0.99,
  
  # marginal=F,
  
  messages = TRUE
  
)
p
ggsave("G:/免疫亚型分型/result8/Biomarker_ICB/表达相关性/ICIs/CTLA4_IGHM_Pro.pdf",
       width = 6,height = 6)

p<-ggscatterstats(
  
  data = pro_ICIs_CS_exp,
  
  x = CCL19,
  
  y = CTLA4,
  
  type = "p",
  
  conf.level = 0.99,
  
  # marginal=F,
  
  messages = TRUE
  
)
ggsave("G:/免疫亚型分型/result8/Biomarker_ICB/表达相关性/ICIs/CTLA4_CCL19_Pro.pdf",
       width = 6,height = 6)
######计算相关性
###治疗前
library(psych)
pre_cortest_psy <- corr.test(CS_pre_exp,CS_pre_exp,method = "pearson")
pre_P<-as.data.frame(pre_cortest_psy$p)
pre_P$CS<-rownames(pre_P)
###宽数据变长数据
library(tidyr)
pre_P_long <- gather(pre_P, key = "signature", value = "P",
                     -`CS`)
pre_P_long<-na.omit(pre_P_long)
pre_P_sig<-pre_P_long[pre_P_long$P<0.05,]

write.table(pre_P_sig,"G:/免疫亚型分型/result8/Biomarker_ICB/pre_P_sig.txt",sep="\t",quote=FALSE,row.names = F)
pre_R<-as.data.frame(pre_cortest_psy$r)
pre_R$CS<-rownames(pre_R)
###宽数据变长数据
pre_R_long <- gather(pre_R, key = "signature", value = "R",
                     -`CS`)
pre_R_long<-na.omit(pre_R_long)
pre_R_sig<-pre_R_long[pre_R_long$R>=0.5,]
####表达相关性热图
library(ggcorrplot)
#install.packages("ggthemes")
library(ggthemes)
corr <- round(cor(CS_pre_exp), 1)
p.mat <- cor_pmat(CS_pre_exp)
ggcorrplot(corr, method = "circle", type = "upper", ggtheme = ggplot2::theme_bw(), title = "", 
           show.legend = TRUE, legend.title = "Corr", show.diag = T, 
           colors = c("#839EDB", "white", "#FF8D8D"), outline.color = "white", 
           hc.order = T, hc.method = "complete", lab = T, lab_col = "black", 
           lab_size = 2, p.mat = p.mat, sig.level = 0.05,  pch = 4, pch.col = "black", pch.cex = 5, tl.cex = 12, 
           tl.col = "black", tl.srt = 45, digits = 2)


corr <- round(cor(CS_pro_exp), 1)
p.mat <- cor_pmat(CS_pro_exp)
ggcorrplot(corr, method = "circle", type = "upper", ggtheme = ggplot2::theme_bw(), title = "", 
           show.legend = TRUE, legend.title = "Corr", show.diag = T, 
           colors = c("#839EDB", "white", "#FF8D8D"), outline.color = "white", 
           hc.order = T, hc.method = "complete", lab = T, lab_col = "black", 
           lab_size = 2, p.mat = p.mat, sig.level = 0.05,  pch = 4, pch.col = "black", pch.cex = 5, tl.cex = 12, 
           tl.col = "black", tl.srt = 45, digits = 2)


###治疗后
library(psych)
pro_cortest_psy <- corr.test(CS_pro_exp,CS_pro_exp,method = "pearson")
pro_P<-as.data.frame(pro_cortest_psy$p)
pro_P$CS<-rownames(pro_P)
###宽数据变长数据
library(tidyr)
pro_P_long <- gather(pro_P, key = "signature", value = "P",
                     -`CS`)
pro_P_long<-na.omit(pro_P_long)
pro_P_sig<-pro_P_long[pro_P_long$P<0.05,]
write.table(pro_P_sig,"G:/免疫亚型分型/result8/Biomarker_ICB/pro_P_sig.txt",sep="\t",quote=FALSE,row.names = F)
pro_R<-as.data.frame(pro_cortest_psy$r)
pro_R$CS<-rownames(pro_R)
###宽数据变长数据
pro_R_long <- gather(pro_R, key = "signature", value = "R",
                     -`CS`)
pro_R_long<-na.omit(pro_R_long)
pro_R_sig<-pro_R_long[pro_R_long$R>=0.5,]




write.table(pro_R_sig,"G:/免疫亚型分型/result8/Biomarker_ICB/pro_R_sig.txt",sep="\t",quote=FALSE,row.names = F)
pro_Corgene<-intersect(pro_P_sig$signature,pro_R_sig$signature)
Signiture_final<-intersect(Pre_Corgene,pro_Corgene)
write.table(Signiture_final,"G:/免疫亚型分型/result8/Biomarker_ICB/Signiture_final.txt",sep="\t",quote=FALSE,row.names = F)
