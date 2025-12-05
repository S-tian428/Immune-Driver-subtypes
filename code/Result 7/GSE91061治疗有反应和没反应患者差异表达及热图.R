options(stringsAsFactors = FALSE)
###清空环境
rm(list=ls())
library(limma)
setwd("G:/免疫亚型分型/result8")
GSE91061_exp<-read.table("GSE91061_exp.txt",
                         header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
GSE91061_pro<-GSE91061_exp[,grep("On",colnames(GSE91061_exp))]
colnames(GSE91061_pro)<-gsub("_.*","",colnames(GSE91061_pro))

Clinical<-read.table("G:/免疫亚型分型/result7/GSE91061/clinical_result.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
# load("G:/免疫亚型分型/result7/GSE91061/yau.ntp.pred.Rdata")
# clust<-yau.ntp.pred$clust.res
# clust$clust<-paste0("CS",clust$clust)
# index<-match(Clinical$X,clust$samID)
# Clinical$cluster<-clust$clust[index]
res_sample<-Clinical$X[Clinical$response%in%"DCB"]
non_sample<-Clinical$X[Clinical$response%in%"NDB"]
#####
expression<-GSE91061_pro[,c(non_sample,res_sample)]
library(limma)
#####获得肿瘤样本和正常样本标筿
group1 <- rep('NoResponse',length(non_sample))
group2 <- rep('Response', length(res_sample))
#######热图中样本的标签及差异表达分析中的标筿
grouplist<-c(group1,group2)
#######差异表达分析
design <- model.matrix(~0+factor(grouplist))
colnames(design)=levels(factor(grouplist))#改design列名为分组信恿
rownames(design)=colnames(expression)
contrast.matrix <- makeContrasts("Response-NoResponse", levels=design)
fit <- lmFit(expression,design)
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效枿
fit2 <- eBayes(fit2)
tempOutput<-topTable(fit2, adjust="BH", coef=1, n=Inf)
write.table(tempOutput,"Response_limma.txt",sep="\t",quote=FALSE)
Sig_respongse_limma<-tempOutput[tempOutput$P.Value<0.05,]

#######CS1
CS1_biomarker<-read.table("G:/免疫亚型分型/result5/CS1_biomarker.txt",
                          header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
CS1_sig<-Sig_respongse_limma[CS1_biomarker$Gene.Symbol,]
CS1_sig<-na.omit(CS1_sig)
CS1_sig$Cluster<-"CS1"
#######CS2
CS2_biomarker<-read.table("G:/免疫亚型分型/result5/CS2_biomarker.txt",
                          header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
CS2_sig<-Sig_respongse_limma[CS2_biomarker$Gene.Symbol,]
CS2_sig<-na.omit(CS2_sig)
CS2_sig$Cluster<-"CS2"
Signiture_result<-rbind(CS1_sig,CS2_sig)
write.table(Signiture_result,"Signiture_result.txt",sep="\t",quote=FALSE)

#######热图
Sig_expression<-expression[rownames(Signiture_result),]
#构建分组信息
annotation_col<-data.frame(grouplist)
rownames(annotation_col) <- colnames(expression)
annotation_row<-data.frame(Signiture_result$Cluster)
rownames(annotation_row)<-rownames(Signiture_result)
colnames(annotation_row)<-"cluster"
rowcolor <- c("#ffff99","#b15928") 
names(rowcolor) <- c("CS1","CS2") #类型颜色

colcolor <- c("#2ec4b6","#e71d36") 
names(colcolor) <- c("NoResponse","Response") #类型颜色

ann_colors<-list(cluster=rowcolor,grouplist=colcolor)
library(pheatmap)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
p<-pheatmap(as.matrix(Sig_expression)
            ,color = c(colorRampPalette(colors = c("#3b95e0","white"))(length(bk)/2),
                       colorRampPalette(colors = c("white","#f78014"))(length(bk)/2)),
            # legend_breaks=seq(-8,8,2),
            breaks=bk,
            cluster_rows = F,cluster_cols = F,fontsize_row=5,
            scale = 'row',show_colnames = F,show_rownames = T,
            annotation_col = annotation_col,
            annotation_row = annotation_row,
            annotation_colors = ann_colors)
p
library(ggpubr)
ggsave("Biomarker_heatmap.pdf",p,width = 8,height = 3)
getwd()
