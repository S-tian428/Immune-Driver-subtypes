options(stringsAsFactors = FALSE)
###清空环境
rm(list=ls())
library(limma)
#读取验证集数据GSE91061
setwd("G:/免疫亚型分型/result7/GSE91061")
######GSE91061的平台ID为NCBI entrezID
GSE91061_exp<-read.table("GSE91061_fpkm.txt",
                         header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
####用clusterProfiler进行ID转换
library(clusterProfiler)
ID<-GSE91061_exp$X
idsTOsymbols = bitr(ID, fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db")
index<-match(GSE91061_exp$X,idsTOsymbols$ENTREZID)
GSE91061_exp$Symbol<-idsTOsymbols$SYMBOL[index]
GSE91061_exp<-GSE91061_exp[!is.na(GSE91061_exp$Symbol),]
rownames(GSE91061_exp)<-GSE91061_exp$Symbol
GSE91061_exp<-GSE91061_exp[,-c(1,ncol(GSE91061_exp))]
GSE91061_exp<-log2(GSE91061_exp+1)
write.table(GSE91061_exp,"G:/免疫亚型分型/result8/GSE91061_exp.txt",sep="\t",quote=FALSE)

####筛选配对样本
# load("G:/免疫亚型分型/result7/GSE91061/yau.ntp.pred.Rdata")
# cluster<-yau.ntp.pred$clust.res
# GSE91061_colnames<-gsub("_.*","",colnames(GSE91061_exp))
# index<-which(GSE91061_colnames%in%cluster$samID)
# GSE91061_exp<-GSE91061_exp[,index]
######两种亚型治疗前后样本差异表达分析
###CS1
GSE91061_paired<-data.frame(sample=colnames(GSE91061_exp))
GSE91061_paired$Group<-"Pro"
GSE91061_paired$Group[grep("Pre",colnames(GSE91061_exp))]<-"Pre"
GSE91061_paired$subject<-gsub("_.*","",GSE91061_paired$sample)
individuals<-factor(GSE91061_paired$subject)
treatment<-factor(GSE91061_paired$Group)
design_paried <- model.matrix(~ individuals + treatment)
fit2 <- lmFit(GSE91061_exp,design_paried)
fit2 <- eBayes(fit2)
tempOutput<-topTable(fit2, adjust="BH", coef=1, n=Inf)
write.table(tempOutput,"G:/免疫亚型分型/result8/ALL_limma.txt",sep="\t",quote=FALSE)
######绘制火山图
library(dplyr)
library(tibble)
data <- 
  tempOutput %>% 
  mutate(change = as.factor(ifelse(adj.P.Val < 0.05 & abs(logFC) > 1,
                                   ifelse(logFC > 1 ,'Up','Down'),'No Change'))) %>% 
  rownames_to_column('gene')
colnames(data)
# 普通火山图
library(ggplot2)
p <- ggplot(data = data, 
            aes(x = logFC, 
                y = -log10(P.Value))) +  # 设置x轴为logFC，y轴为-P.Value的对数值
  geom_point(alpha = 0.4, size = 3.5, 
             aes(color = change)) +      # 添加散点，根据change列着色
  ylab("-log10(Pvalue)") +               # 设置y轴标签
  scale_color_manual(values = c( "grey", "#f2672a")) +  # 设置颜色映射，蓝色表示下调，灰色表示稳定，红色表示上调
  geom_vline(xintercept = c(-1,1), lty = 4, col = "black", lwd = 0.8) +  # 添加垂直参考线，用于标记logFC阈值
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.8) +   # 添加水平参考线，用于标记-P.Value阈值
  theme_bw()  # 使用网格白底主题
p
####加载CS1 biomarker
CS1_biomarker<-read.table("G:/免疫亚型分型/result5/CS1_biomarker.txt",
                          header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
####加载CS2 biomarker
CS2_biomarker<-read.table("G:/免疫亚型分型/result5/CS2_biomarker.txt",
                          header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
biomarkers<-c(CS1_biomarker$Gene.Symbol,CS2_biomarker$Gene.Symbol)
# 标签添加
# 创建一个新的列label，并初始化为NA
data$label <- ""  
#data<-na.omit(data)
# 根据symbol的值，为特定基因添加标签信息
sig_data<-data[data$adj.P.Val<0.05&abs(data$logFC)>1,]
biomarker<-intersect(sig_data$gene,biomarkers)
index<-match(biomarker,data$gene)
data[index,]$label<-biomarker
library(ggrepel)
p1<-p +  # 基于普通火山图p
  geom_text_repel(data = data, aes(x = logFC, 
                                   y = -log10(P.Value), label = label),
                  max.overlaps = 2000)
p1
ggsave("G:/免疫亚型分型/result8/All_volcano.pdf",p1,width = 9.5,height = 8)
Biomarker_limma<-sig_data[sig_data$gene%in%biomarkers,]
write.table(Biomarker_limma,"G:/免疫亚型分型/result8/Biomarker_limma_Pre_On.txt",sep="\t",quote=FALSE,row.names = F)
