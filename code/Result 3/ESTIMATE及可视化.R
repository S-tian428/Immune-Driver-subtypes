rm(list=ls())
load("I:/免疫亚型分型/MOVICS/STAD.multiomics.rda")
setwd("I:/免疫亚型分型/result3/Estimate/")
surv.info<-STAD_multiomics$surv.info
STAD.FPKM<-STAD_multiomics$STAD.FPKM
write.table(STAD.FPKM,"STAD.FPKM.txt",sep="\t",quote=FALSE)
###########2、Estimate分析#####
##estimate##
library(utils)
#install.packages("estimate")
#install.packages("estimate",
#                 repos="http://r-forge.r-project.org",
#                dependencies=TRUE)
library(estimate)
#设置输入输出文件
# 手动添加行名作为一列
##将准备好的表达谱保存为txt格式，用genesymbol,改成id="GeneSymbol"即可
filterCommonGenes(input.f="STAD.FPKM.txt", output.f="STAD_expression.gct", id="GeneSymbol")
estimateScore(input.ds="STAD_expression.gct", output.ds="expression_estimate_score.gct",
              platform="illumina")
estimate_score <- read.table("expression_estimate_score.gct",  skip = 2, header = TRUE)
rownames(estimate_score)<-estimate_score[,1]
estimate_score<-estimate_score[,-(1:2)]
colnames(estimate_score)<-gsub("\\.","-",colnames(estimate_score))
estimate_score<-as.data.frame(t(estimate_score))
####计算肿瘤纯度
estimate_score$Purity<- cos(0.6049872018+0.0001467884 * estimate_score$ESTIMATEScore)
#####获得患者的亚型分组信息
index<-match(rownames(estimate_score),rownames(surv.info))
estimate_score$Cluster<-surv.info$Cluster[index]
#绘制柱形图
####################StromalScore
library(ggpubr)
p<-ggplot(estimate_score, aes(Cluster, StromalScore)) +
  geom_boxplot(aes(fill = Cluster), notch = FALSE, size = 0.4) +
  scale_fill_manual(values=c("#8DBC80FF","#17692CFF")) +
  guides(fill=guide_legend(title="Cluster")) +
  theme_bw()
p1<-p+
  stat_compare_means(method = "wilcox.test",label = "p.signif",
                     label.x.npc = "center",size=3.5)
p1
ggsave("StromalScore.pdf", p1,width=4.5,height = 4.3)
#####免疫得分
p<-ggplot(estimate_score, aes(Cluster, ImmuneScore)) +
  geom_boxplot(aes(fill = Cluster), notch = FALSE, size = 0.4) +
  scale_fill_manual(values=c("#CAB2D6FF","#6A3D9AFF")) +
  guides(fill=guide_legend(title="Cluster")) +
  theme_bw()
p1<-p+stat_compare_means(method = "wilcox.test",label = "p.signif",
                         label.x.npc = "center",size=3.5)
p1
ggsave("ImmuneScore.pdf", p1,width=4.5,height = 4.3)

####################ESTIMATEScore
library(ggpubr)
p<-ggplot(estimate_score, aes(Cluster, ESTIMATEScore)) +
  geom_boxplot(aes(fill = Cluster), notch = FALSE, size = 0.4) +
  scale_fill_manual(values=c("#8DBC80FF","#17692CFF")) +
  guides(fill=guide_legend(title="Cluster")) +
  theme_bw()
p1<-p+stat_compare_means(method = "wilcox.test",label = "p.signif",
                         label.x.npc = "center",size=3.5)
p1
ggsave("ESTIMATEScore.pdf", p1,width=4.5,height = 4.3)
#####免疫得分
p<-ggplot(estimate_score, aes(Cluster, Purity)) +
  geom_boxplot(aes(fill = Cluster), notch = FALSE, size = 0.4) +
  scale_fill_manual(values=c("#CAB2D6FF","#6A3D9AFF")) +
  guides(fill=guide_legend(title="Cluster")) +
  theme_bw()
p1<-p+stat_compare_means(method = "wilcox.test",label = "p.signif",
                         label.x.npc = "center",size=3.5)
p1
ggsave("Purity.pdf", p1,width=4.5,height = 4.3)

