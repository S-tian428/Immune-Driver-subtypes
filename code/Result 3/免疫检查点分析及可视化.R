##免疫细胞浸润
#####Cibersort
rm(list=ls())
load("I:/免疫亚型分型/MOVICS/STAD.multiomics.rda")
setwd("I:/免疫亚型分型/result3/ICGs/")
STAD.FPKM<-STAD_multiomics$STAD.FPKM
surv.info<-STAD_multiomics$surv.info
surv.info<-surv.info[order(surv.info$Cluster),]
STAD.FPKM<-STAD.FPKM[,rownames(surv.info)]
####CS1和CS2组间进行所有基因的进行差异表达分析
####limma差异表达分析
library(limma)
group1 <- surv.info$Cluster[surv.info$Cluster%in%"CS1"]
group2 <- surv.info$Cluster[surv.info$Cluster%in%"CS2"]
#######热图中样本的标签及差异表达分析中的标签
grouplist<-c(group1,group2)
#######差异表达分析
design <- model.matrix(~0+factor(grouplist))
colnames(design)=levels(factor(grouplist))#改design列名为分组信息
rownames(design)=colnames(STAD.FPKM)
contrast.matrix <- makeContrasts("CS2-CS1", levels=design)
fit <- lmFit(STAD.FPKM,design)
fit2 <- contrasts.fit(fit, contrast.matrix) ##这一步很重要，大家可以自行看看效果
fit2 <- eBayes(fit2)
tempOutput<-topTable(fit2, adjust="BH", coef=1, n=Inf)
limma_result<-na.omit(tempOutput)
Sig_limma<-limma_result[limma_result$adj.P.Val<0.05,]
Sig_limma<-Sig_limma[order(Sig_limma$logFC),]
write.table(Sig_limma,"Limma_signaficant.txt",sep = "\t",quote = F)

####一些免疫因子的差异情况
Sig_limma<-read.table("Limma_signaficant.txt",
                            header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)

###差异趋化因子
Chemokine<-read.table("I:/免疫亚型分型/result3/Chemokine/36860872.txt",
                            header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)

Chemokine<-Chemokine$gene

####差异免疫检查点基因
checkpoint_gene<-read.table("Immune checkpoint gene.txt",
                            header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
ICGs<-checkpoint_gene$Symbol
####ICGs&Chemokine(没有)
ICGs.Chemokine<-intersect(ICGs,Chemokine)

####只有一部分显著
index_Chemokine<-which(rownames(Sig_limma)%in%Chemokine)
Chemokine_limma<-Sig_limma[index_Chemokine,]
write.table(Chemokine_limma,"I:/免疫亚型分型/result3/Chemokine/Chemokine_limma.txt",sep="\t",quote=FALSE)


####只有一部分显著
index<-which(rownames(Sig_limma)%in%ICGs)
ICGs_limma<-Sig_limma[index,]
write.table(ICGs_limma,"ICGs_limma.txt",sep="\t",quote=FALSE)
#####绘制热图
ICG_exp<-STAD.FPKM[rownames(ICGs_limma),]
Chemokine_exp<-STAD.FPKM[rownames(Chemokine_limma),]
ICGs.Chemokine_exp<-rbind(ICG_exp,Chemokine_exp)

#构建分组信息
annotation_col<-data.frame(clust=surv.info$Cluster)
rownames(annotation_col) <- rownames(surv.info)

annotation_row<-data.frame(Genes=c(rep("Immune Checkpoint Genes",nrow(ICG_exp)),
                           rep("Chemokine Genes",nrow(Chemokine_exp))))
rownames(annotation_row)<-rownames(ICGs.Chemokine_exp)

rowcolor <- c("#ffff99","#b15928") 
names(rowcolor) <- c("Immune Checkpoint Genes","Chemokine Genes") #类型颜色

colcolor <- c("#2ec4b6","#e71d36") 
names(colcolor) <- c("CS1","CS2") #类型颜色

ann_colors<-list(Genes=rowcolor,clust=colcolor)
#mycol<-colorRampPalette(c("#17692c", "white", "#ff7f00"))(100)
#mycol<-colorRampPalette(c("#6a3d9a", "white", "#b15928"))(100)
library(pheatmap)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
p<-pheatmap(ICGs.Chemokine_exp
            ,color = c(colorRampPalette(colors = c("#3b95e0","white"))(length(bk)/2),
                       colorRampPalette(colors = c("white","#f78014"))(length(bk)/2)),
           # legend_breaks=seq(-8,8,2),
            breaks=bk,
            cluster_rows = F,cluster_cols = F,fontsize_row=5,
            scale = 'row',show_colnames = F,show_rownames = T,
            annotation_col = annotation_col,
            annotation_row = annotation_row,
            annotation_colors = ann_colors)
library(ggpubr)
ggsave("Chemokine.ICGs_heatmap.pdf",p,width = 7,height = 5)
#####单个基因表达小提琴图
Sig_genes<-rownames(ICGs_limma)[abs(ICGs_limma$logFC)>1]
Gene_exp<-STAD.FPKM[Sig_genes,]
Gene_exp<-as.data.frame(t(Gene_exp))
Gene_exp$Cluster<-surv.info$Cluster
####小提琴图
library(ggplot2)
i=6
Gene_exp$Cluster<-as.factor(Gene_exp$Cluster)
for (i in 1:length(Sig_genes)) {
  gene<-colnames(Gene_exp)[i]
  p<-ggplot(Gene_exp, aes(x=Cluster, y=Gene_exp[,i],color=Cluster))+geom_violin(linewidth=1)+
    scale_color_brewer(palette="Dark2")
  p
  p<-p+stat_summary(fun=mean, geom="point", size=5,color="darkgrey")+
    theme_classic(base_size = 15)+ylab(gene)+xlab("")
  p<-p+stat_compare_means(method = "wilcox.test")
  p
  myfilename <-paste0(gene,"_boxplot.pdf")
  #library(ggplot2)
  ggsave(myfilename,p,width=4.3,height = 3.7)
  
  print(i)
}

#####火山图
library(tibble)
library(ggrepel)
library(dplyr)
data <- 
  Chemokine_limma %>% 
  mutate(change = as.factor(ifelse(P.Value < 0.05 & abs(logFC) > 1,
                                   ifelse(logFC > 1 ,'Up','Down'),'No change'))) %>% 
  rownames_to_column('gene')
colnames(data)
p<-ggplot(data,aes(logFC, -log10(adj.P.Val), label = gene))+
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#999999")+
  geom_vline(xintercept = c(-1,1), linetype = "dashed", color = "#999999")+
  geom_point(aes(color = change),
             size = 2
             ) +
  theme_bw(base_size = 12)+
  ggsci::scale_color_jama() +
  theme(panel.grid = element_blank(),
        legend.position = 'right') +
  # 添加标签
  geom_text_repel(data = filter(data,abs(logFC)> 1),
                  max.overlaps = getOption("ggrepel.max.overlaps", default = 20),
                  aes(label = gene, 
                      color = change),
                  size = 2) +
  xlab("Log2FC")+
  ylab("-Log10(FDR q-value)")

ggsave("I:/免疫亚型分型/result3/Chemokine/volcano.pdf",p,width = 7,height = 3)

###差异免疫循环因子
CICs<-read.table("I:/免疫亚型分型/result3/Cancer Immunity Cycle/Cancer_Immunity_Cycle.txt",
                 header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)

CICs_positive<-CICs$GeneSymbol[CICs$Direction%in%"positive"]
CICs_negative<-CICs$GeneSymbol[CICs$Direction%in%"negative"]

####分成三组
####CICs_positive&CICs_negative
CICs_positive.negative<-intersect(CICs_positive,CICs_negative)
###CICs_positive
CICs_positive<-setdiff(CICs_positive,CICs_positive.negative)
###CICs_negative
CICs_negative<-setdiff(CICs_negative,CICs_positive.negative)
####只有一部分显著
index_CICs_positive<-which(rownames(Sig_limma)%in%CICs_positive)
CICs_positive_limma<-Sig_limma[index_CICs_positive,]
index_CICs_negative<-which(rownames(Sig_limma)%in%CICs_negative)
CICs_negative_limma<-Sig_limma[index_CICs_negative,]
index_CICs_positive.negative<-which(rownames(Sig_limma)%in%CICs_positive.negative)
CICs_positive.negative_limma<-Sig_limma[index_CICs_positive.negative,]
write.table(CICs_positive_limma,"I:/免疫亚型分型/result3/Cancer Immunity Cycle/Positive_limma.txt",sep="\t",quote=FALSE)
write.table(CICs_negative_limma,"I:/免疫亚型分型/result3/Cancer Immunity Cycle/Negative_limma.txt",sep="\t",quote=FALSE)
write.table(CICs_positive.negative_limma,"I:/免疫亚型分型/result3/Cancer Immunity Cycle/Positive.Negative_limma.txt",sep="\t",quote=FALSE)
####表达
CICs_pos_exp<-STAD.FPKM[rownames(CICs_positive_limma),]
CICs_neg_exp<-STAD.FPKM[rownames(CICs_negative_limma),]
CICs_pos.neg_exp<-STAD.FPKM[rownames(CICs_positive.negative_limma),]

CICs_exp<-rbind(CICs_pos_exp,CICs_neg_exp,CICs_pos.neg_exp)

#构建分组信息
annotation_col<-data.frame(clust=surv.info$Cluster)
rownames(annotation_col) <- rownames(surv.info)

annotation_row<-data.frame(Genes=c(rep("Positive",nrow(CICs_pos_exp)),
                                   rep("Negative",nrow(CICs_neg_exp)),
                                   rep("Positive & Negative",nrow(CICs_pos.neg_exp))))
rownames(annotation_row)<-rownames(CICs_exp)

rowcolor <- c("#17692c","#6a3d9a","#1f78b4") 
names(rowcolor) <- c("Positive","Negative","Positive & Negative") #类型颜色

colcolor <- c("#2ec4b6","#e71d36") 
names(colcolor) <- c("CS1","CS2") #类型颜色

ann_colors<-list(Genes=rowcolor,clust=colcolor)
#mycol<-colorRampPalette(c("#17692c", "white", "#ff7f00"))(100)
#mycol<-colorRampPalette(c("#6a3d9a", "white", "#b15928"))(100)
library(pheatmap)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))
p<-pheatmap(CICs_exp
            ,color = c(colorRampPalette(colors = c("#3b95e0","white"))(length(bk)/2),
                       colorRampPalette(colors = c("white","#f78014"))(length(bk)/2)),
            # legend_breaks=seq(-8,8,2),
            breaks=bk,
            cluster_rows = F,cluster_cols = F,fontsize_row=5,
            scale = 'row',show_colnames = F,show_rownames = T,
            annotation_col = annotation_col,
            annotation_row = annotation_row,
            annotation_colors = ann_colors)
library(ggpubr)
ggsave("I:/免疫亚型分型/result3/Cancer Immunity Cycle/CICs_heatmap.pdf",p,width = 7,height = 5)
#####单个基因表达小提琴图

