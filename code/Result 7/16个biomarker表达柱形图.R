rm(list=ls())
#library(GSVA)
CS1_biomarker<-read.table("G:/免疫亚型分型/result5/CS1_biomarker.txt",
                          header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
####加载CS2 biomarker
CS2_biomarker<-read.table("G:/免疫亚型分型/result5/CS2_biomarker.txt",
                          header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
biomarker<-c(CS1_biomarker$Gene.Symbol,CS2_biomarker$Gene.Symbol)

Res_NoRes_gene<-read.table("G:/免疫亚型分型/result8/Biomarker_ICB/Response_NoResponse_gene.txt",
                           header=F,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
Res_NoRes_gene<-Res_NoRes_gene$V2

Pre_Pro_gene<-read.table("G:/免疫亚型分型/result8/Biomarker_ICB/Pre_Pro_gene.txt",
                         header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
Pre_Pro_gene<-Pre_Pro_gene$x
#####重要的CS_biomarker
Res_NoRes_biomarker<-intersect(Res_NoRes_gene,biomarker)
Pre_Pro_biomarker<-intersect(Pre_Pro_gene,biomarker)
CS_biomarker<-union(Res_NoRes_biomarker,Pre_Pro_biomarker)
write.table(CS_biomarker,"G:/免疫亚型分型/result8/Biomarker_ICB/ICB_CS_biomarker.txt",sep="\t",quote=FALSE,row.names = F)

####单个基因表达柱形图
setwd("G:/免疫亚型分型/result8")
GSE91061_exp<-read.table("GSE91061_exp.txt",
                         header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
Clinical<-read.table("G:/免疫亚型分型/result7/GSE91061/clinical_result.txt",
                     header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
res_sample<-Clinical$X[Clinical$response%in%"DCB"]
non_sample<-Clinical$X[Clinical$response%in%"NDB"]
sample<-gsub("_.*","",colnames(GSE91061_exp))
index<-which(sample%in%Clinical$X)
GSE91061_exp<-GSE91061_exp[,index]
#####Pre/On样本柱形图
####pre
GSE91061_pre<-GSE91061_exp[,grep("Pre",colnames(GSE91061_exp))]
colnames(GSE91061_pre)<-gsub("_.*","",colnames(GSE91061_pre))
pre_expression<-GSE91061_pre[,c(non_sample,res_sample)]
pre_expression<-as.data.frame(t(pre_expression))
pre_biomarker_expression<-pre_expression[,Pre_Pro_biomarker]
index<-match(rownames(pre_biomarker_expression),Clinical$X)
pre_biomarker_expression$response<-Clinical$response[index]
######表达柱形图
library(ggpubr)
j=1
setwd("G:/免疫亚型分型/result8/Biomarker_ICB/Pre_Pro_biomarker/Pre")
for (j in 1:(ncol(pre_biomarker_expression)-1)) {
  p = ggplot(pre_biomarker_expression, aes(x = response, y = pre_biomarker_expression[,j], fill=response)) + 
    labs(x="", y="log2(FPKM+1)",title = colnames(pre_biomarker_expression)[j]) +
    theme_bw(base_size = 9) + ###去除背景颜色
    #geom_text(aes(label=Hugo_Symbol),check_overlap = TRUE,size=6) +
    geom_boxplot(width=0.5, outlier.colour = NA) +
    # 不显示离群点
    geom_jitter(width = 0.1) +
    theme(axis.text = element_text(size = 20))+ 
    #theme(axis.title.x=element_text(vjust=2, size=20,face = "bold"))+
    # 添加散点
    # geom_line(aes(group=Hugo_Symbol),
    #           color="black", alpha=1,linetype=2,
    #           linewidth=0.8) +
    theme(axis.title.y=element_text(vjust=2, size=20))+
    scale_fill_manual(values = c('#fdd692','#ec7357'))+
    stat_compare_means( method = "wilcox.test")
  p
  ggsave(paste0(colnames(pre_biomarker_expression)[j],".pdf"),p,height = 6,width = 6)
}
####PRO
####pre
GSE91061_pro<-GSE91061_exp[,grep("On",colnames(GSE91061_exp))]
colnames(GSE91061_pro)<-gsub("_.*","",colnames(GSE91061_pro))
pro_expression<-GSE91061_pro[,c(non_sample,res_sample)]
pro_expression<-as.data.frame(t(pro_expression))
pro_biomarker_expression<-pro_expression[,Pre_Pro_biomarker]
index<-match(rownames(pro_biomarker_expression),Clinical$X)
pro_biomarker_expression$response<-Clinical$response[index]
######表达柱形图
library(ggpubr)
j=1
setwd("G:/免疫亚型分型/result8/Biomarker_ICB/Pre_Pro_biomarker/Pro")
for (j in 1:(ncol(pro_biomarker_expression)-1)) {
  p = ggplot(pro_biomarker_expression, aes(x = response, y = pro_biomarker_expression[,j], fill=response)) + 
    labs(x="", y="log2(FPKM+1)",title = colnames(pro_biomarker_expression)[j]) +
    theme_bw(base_size = 9) + ###去除背景颜色
    #geom_text(aes(label=Hugo_Symbol),check_overlap = TRUE,size=6) +
    geom_boxplot(width=0.5, outlier.colour = NA) +
    # 不显示离群点
    geom_jitter(width = 0.1) +
    theme(axis.text = element_text(size = 20))+ 
    #theme(axis.title.x=element_text(vjust=2, size=20,face = "bold"))+
    # 添加散点
    # geom_line(aes(group=Hugo_Symbol),
    #           color="black", alpha=1,linetype=2,
    #           linewidth=0.8) +
    theme(axis.title.y=element_text(vjust=2, size=20))+
    scale_fill_manual(values = c('#fdd692','#ec7357'))+
    stat_compare_means( method = "wilcox.test")
  p
  ggsave(paste0(colnames(pro_biomarker_expression)[j],".pdf"),p,height = 6,width = 6)
}
#####Presponse
###配对样本差异表达分析
sample<-gsub("_.*","",colnames(GSE91061_exp))
index<-which(sample%in%res_sample)
Presonse_exp<-GSE91061_exp[,index]
Presonse_exp<-data.frame(t(Presonse_exp))
Presonse_biomarker_expression<-Presonse_exp[,Res_NoRes_biomarker]
Presonse_biomarker_expression$sample<-gsub("_.*","",rownames(Presonse_biomarker_expression))
Presonse_biomarker_expression$Group<-"Pro"
Presonse_biomarker_expression$Group[grep("Pre",rownames(Presonse_biomarker_expression))]<-"Pre"

j=1
setwd("G:/免疫亚型分型/result8/Biomarker_ICB/Res_NoRes_biomarker/Res")
for (j in 1:(ncol(Presonse_biomarker_expression)-2)) {
  p = ggplot(Presonse_biomarker_expression, aes(x = Group, y = Presonse_biomarker_expression[,j], fill=Group)) + 
    labs(x="", y="log2(FPKM+1)",title = colnames(Presonse_biomarker_expression)[j]) +
    theme_bw(base_size = 9) + ###去除背景颜色
    #geom_text(aes(label=Hugo_Symbol),check_overlap = TRUE,size=6) +
    geom_boxplot(outlier.size = 0.5) +
    geom_line(aes(group = sample), color = "grey80", size = 1) +
    geom_point(size = 2) +
    theme(axis.text = element_text(size = 20))+ 
    #theme(axis.title.x=element_text(vjust=2, size=20,face = "bold"))+
    # 添加散点
    # geom_line(aes(group=Hugo_Symbol),
    #           color="black", alpha=1,linetype=2,
    #           linewidth=0.8) +
    theme(axis.title.y=element_text(vjust=2, size=20))+
    scale_fill_manual(values = c('#fdd692','#ec7357'))+
    stat_compare_means( method = "wilcox.test")
  p
  ggsave(paste0(colnames(Presonse_biomarker_expression)[j],".pdf"),p,height = 6,width = 6)
}


#####NonPresponse
###配对样本差异表达分析
index<-which(sample%in%non_sample)
NoPresonse_exp<-GSE91061_exp[,index]
NoPresonse_exp<-data.frame(t(NoPresonse_exp))
NoPresonse_biomarker_expression<-NoPresonse_exp[,Res_NoRes_biomarker]
NoPresonse_biomarker_expression$sample<-gsub("_.*","",rownames(NoPresonse_biomarker_expression))
NoPresonse_biomarker_expression$Group<-"Pro"
NoPresonse_biomarker_expression$Group[grep("Pre",rownames(NoPresonse_biomarker_expression))]<-"Pre"

j=1
setwd("G:/免疫亚型分型/result8/Biomarker_ICB/Res_NoRes_biomarker/NoRes")
for (j in 1:(ncol(NoPresonse_biomarker_expression)-2)) {
  p = ggplot(NoPresonse_biomarker_expression, aes(x = Group, y = NoPresonse_biomarker_expression[,j], fill=Group)) + 
    labs(x="", y="log2(FPKM+1)",title = colnames(NoPresonse_biomarker_expression)[j]) +
    theme_bw(base_size = 9) + ###去除背景颜色
    #geom_text(aes(label=Hugo_Symbol),check_overlap = TRUE,size=6) +
    geom_boxplot(outlier.size = 0.5) +
    geom_line(aes(group = sample), color = "grey80", size = 1) +
    geom_point(size = 2) +
    theme(axis.text = element_text(size = 20))+ 
    #theme(axis.title.x=element_text(vjust=2, size=20,face = "bold"))+
    # 添加散点
    # geom_line(aes(group=Hugo_Symbol),
    #           color="black", alpha=1,linetype=2,
    #           linewidth=0.8) +
    theme(axis.title.y=element_text(vjust=2, size=20))+
    scale_fill_manual(values = c('#fdd692','#ec7357'))+
    stat_compare_means( method = "wilcox.test")
  p
  ggsave(paste0(colnames(NoPresonse_biomarker_expression)[j],".pdf"),p,height = 6,width = 6)
}


