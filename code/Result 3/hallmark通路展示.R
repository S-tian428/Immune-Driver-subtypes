rm(list = ls())
###33种癌症样本的hallmark结果
ssGSEA_hallmark<-read.table("H:/免疫亚型分型/result3/hallmark.pathway/hallmark_result.txt",
                            header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
#ssGSEA_hallmark<-as.data.frame(t(ssGSEA_hallmark))
colnames(ssGSEA_hallmark)<-gsub("\\.","-",colnames(ssGSEA_hallmark))
load("H:/免疫亚型分型/MOVICS/STAD.multiomics.rda")
surv.info<-STAD_multiomics$surv.info

# index<-match(rownames(ssGSEA_hallmark),rownames(surv.info))
# ssGSEA_hallmark$Cluster<-surv.info$Cluster[index]
#####分突变样本和野生型样本
###突变型样本
#####cs1
CS1_sample<-rownames(surv.info)[surv.info$Cluster%in%"CS1"]
ssGSEA_hallmark_CS1<-ssGSEA_hallmark[,CS1_sample]
#####野生型样本
CS2_sample<-rownames(surv.info)[surv.info$Cluster%in%"CS2"]
ssGSEA_hallmark_CS2<-ssGSEA_hallmark[,CS2_sample]
exprSet<-cbind(ssGSEA_hallmark_CS1,ssGSEA_hallmark_CS2)
####野生型和突变型样本数要从矩阵中提取
group1 <- rep('CS1',length(CS1_sample))
group2 <- rep('CS2', length(CS2_sample))
grouplist<-c(group1,group2)

design <- model.matrix(~0+factor(grouplist))
colnames(design)=levels(factor(grouplist))#??design????Ϊ??????Ϣ
rownames(design)=colnames(exprSet)
library(limma)
contrast.matrix <- makeContrasts("CS2-CS1", levels=design)
fit <- lmFit(exprSet,design)
fit2 <- contrasts.fit(fit, contrast.matrix) ##??һ??????Ҫ?????ҿ??????п???Ч??
fit2 <- eBayes(fit2)
tempOutput<-topTable(fit2, adjust="BH", coef=1, n=Inf)
result<-na.omit(tempOutput)
result_sig<-result[result$adj.P.Val<0.05,]

result_sig<-result_sig[order(result_sig$logFC,decreasing = T),]
write.table(result,"H:/免疫亚型分型/result3/hallmark.pathway/limma_result.txt",sep="\t",quote=FALSE,row.names = T)
write.table(result_sig,"H:/免疫亚型分型/result3/hallmark.pathway/Limma_result_sig.txt",sep="\t",quote=FALSE,row.names = T)

#####绘制上下调P值热图
result_sig$mutation<-NA
result_sig$wild<-NA
result_sig$mutation[result_sig$logFC>0]=-log(result_sig$adj.P.Val[result_sig$logFC>0])
result_sig$wild[result_sig$logFC<0]=-log(result_sig$adj.P.Val[result_sig$logFC<0])
#####拿前6个基因做柱形图
sig_pathway<-rownames(result_sig)[1:6]

sig_pathway_mut<-t(ssGSEA_hallmark_mut)[,sig_pathway]
sig_pathway_mut<-data.frame(sig_pathway_mut)
sig_pathway_mut$Type<-"MUT"
###野生型
sig_pathway_wild<-t(ssGSEA_hallmark_wild)[,sig_pathway]
sig_pathway_wild<-data.frame(sig_pathway_wild)
sig_pathway_wild$Type<-"Wild"
####合并
sig_pathway<-rbind(sig_pathway_mut,sig_pathway_wild)
colnames(sig_pathway)


######格式
k=1
sig_box_result<-NULL
for(k in 1:6){
  sig_box<-data.frame(sig_pathway[,k],
                      sig_pathway[,7],
                      colnames(sig_pathway)[k])
  colnames(sig_box)<-c("Hallmark_score","Type","Pathway")
  sig_box_result<-rbind(sig_box_result,sig_box)
}
library(ggplot2)
library(ggpubr)
library(patchwork)
p <- ggboxplot(sig_box_result, x = "Type", y = "Hallmark_score",
               color = "Type", palette = "aaas",
               #add = "jitter",
               facet.by = "Pathway", 
               short.panel.labs = F)
p
p1 =p + stat_compare_means(label =  "p.signif", label.x = 1.5)
p1

p <- ggboxplot(sig_pathway, x = "Type", y = "HALLMARK_E2F_TARGETS",
               # 配色方案 ?ggboxplot
               color = "Type", palette = "aaas")
p
#  Add p-value
p1 = p + stat_compare_means() #default Wilcoxon

p2 = p + stat_compare_means(method = "t.test")
p1 + p2


p <- ggboxplot(sig_pathway, x = "supp", y = "len",
               color = "supp", palette = "jco",
               #add = "jitter",
               facet.by = "dose", 
               short.panel.labs = F)
p1 = p + stat_compare_means(label = "p.format")
# p + stat_compare_means(label =  "p.signif", label.x = 1.5)

p <- ggboxplot(ToothGrowth, x = "dose", y = "len",
               color = "supp", palette = "jco")
p2 = p + stat_compare_means(aes(group = supp))










logP_halmark<-data.frame(result_sig$mutation,result_sig$wild)
rownames(logP_halmark)<-rownames(result_sig)
logP_halmark<-as.matrix(logP_halmark)
p<-pheatmap(logP_halmark,
            na_col = "grey",
            show_colnames = F,cellwidth = 10,
            cluster_col = F,cluster_rows =F,
            col = colorRampPalette(c("#1f77b4", "white", "#b10026"))(50)
)
####将47个显著通路按照FC指排序
result<-result[order(result$logFC,decreasing = T),]
exprSet_result<-exprSet[rownames(result),]
####样本标签注释一下
annotation_col = data.frame(
  PatientsType = factor(c(group1, group2))
)
rownames(annotation_col) =colnames(exprSet)
#####设置注释样本信息的颜色
ann_colors = list(
  PatientsType = c(Mutation = "#b10026", Wild = "#3b4992")
)
data=as.matrix(exprSet_result)
library(RColorBrewer)
#install.packages("pheatmap")
library(pheatmap)
library(ggplot2)
#breaks
bk <- c(seq(-1,-0.1,by=0.01),seq(0,1,by=0.01))
# 做热图：
library(ComplexHeatmap)
pheatmap(data, annotation_col = annotation_col,
         use_raster = FALSE,scale="none",
         cluster_col = F,cluster_rows =F,
         show_colnames =F)
p<-pheatmap(data,annotation_col = annotation_col,
            annotation_colors = ann_colors,
            scale = "row",na_col = "grey",
            cluster_col = F,cluster_rows =F,
            show_colnames =F,
            col = colorRampPalette(c("#1f77b4", "white", "#b10026"))(50)
)
p
