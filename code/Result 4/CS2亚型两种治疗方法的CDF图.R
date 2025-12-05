rm(list=ls())
####TIDE输入前需要对表达谱进行归一化处理
load("G:/免疫亚型分型/MOVICS/STAD.multiomics.rda")
surv.info<-STAD_multiomics$surv.info
surv.info$patients<-gsub("-01A|-01B","",rownames(surv.info))
TICA_STAD <- read.table("G:/免疫亚型分型/result5/TICA/TCIA-ClinicalData.tsv", header=TRUE, sep="\t")
index<-match(TICA_STAD$barcode,surv.info$patients)
TICA_STAD$Cluster<-surv.info$Cluster[index]
TICA_STAD<-TICA_STAD[!is.na(TICA_STAD$Cluster),]
TICA_STAD<-TICA_STAD[,c(1,22,23,24,25,27)]
colnames(TICA_STAD)
library(ggplot2)
library(ggpubr)
pd1_pos_CS2<-TICA_STAD$ips_ctla4_neg_pd1_pos[TICA_STAD$Cluster%in%"CS2"]
CS2_ips1<-data.frame(IPS=pd1_pos_CS2,Catogary="anti-PD-1")
pd1_ctla4_CS2<-TICA_STAD$ips_ctla4_pos_pd1_pos[TICA_STAD$Cluster%in%"CS2"]
CS2_ips2<-data.frame(IPS=pd1_ctla4_CS2,Catogary="anti_PD-1&anti_CTLA-4")
CS2_ips<-rbind(CS2_ips1,CS2_ips2)
####绘制概率分布图
p<-ggplot(CS2_ips, aes(x=IPS, col=Catogary)) + 
  stat_ecdf(geom="smooth", se=F, size=1.2) + 
  theme_bw() 
ggsave("G:/免疫亚型分型/result5/TICA/CDF.pdf",p,width = 5,height = 3)
wt<- wilcox.test(IPS ~ Catogary, data = CS2_ips)
wt$p.value


