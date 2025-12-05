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
TICA_STAD$ips_ctla4_pos_pd1_pos
colnames(TICA_STAD)
library(ggplot2)
library(ggsci)
library(ggpubr)
p = ggplot(TICA_STAD, aes(x = Cluster, y = ips_ctla4_pos_pd1_pos, fill=Cluster)) + 
  labs(x="", y="IPS-CTLA4(+)/PD-1(+)") +
  theme_bw(base_size = 9) + ###去除背景颜色
  #geom_text(aes(label=Hugo_Symbol),check_overlap = TRUE,size=6) +
  geom_boxplot(width=0.5) +
  # 不显示离群点
  #geom_jitter(width = 0.1) +
  theme(axis.text = element_text(size = 20))+ 
  #theme(axis.title.x=element_text(vjust=2, size=20,face = "bold"))+
  # 添加散点
  # geom_line(aes(group=Hugo_Symbol),
  #           color="black", alpha=1,linetype=2,
  #           linewidth=0.8) +
  theme(axis.title.y=element_text(vjust=2, size=20))+
  scale_fill_manual(values = c('#cab2d6','#ffff99'))+
  stat_compare_means( method = "wilcox.test",label.x.npc = "center",size=7)
p
ggsave("G:/免疫亚型分型/result5/TICA/CTLA4+_PD-1+.pdf",p,height = 5,width = 5)
p = ggplot(TICA_STAD, aes(x = Cluster, y = ips_ctla4_neg_pd1_pos, fill=Cluster)) + 
  labs(x="", y="IPS-CTLA4(-)/PD-1(+)") +
  theme_bw(base_size = 9) + ###去除背景颜色
  #geom_text(aes(label=Hugo_Symbol),check_overlap = TRUE,size=6) +
  geom_boxplot(width=0.5) +
  # 不显示离群点
  #geom_jitter(width = 0.1) +
  theme(axis.text = element_text(size = 20))+ 
  #theme(axis.title.x=element_text(vjust=2, size=20,face = "bold"))+
  # 添加散点
  # geom_line(aes(group=Hugo_Symbol),
  #           color="black", alpha=1,linetype=2,
  #           linewidth=0.8) +
  theme(axis.title.y=element_text(vjust=2, size=20))+
  scale_fill_manual(values = c('#cab2d6','#ffff99'))+
  stat_compare_means( method = "wilcox.test",label.x.npc = "center",size=7)
p
ggsave("G:/免疫亚型分型/result5/TICA/CTLA4-_PD-1+.pdf",p,height = 5,width = 5)

p = ggplot(TICA_STAD, aes(x = Cluster, y = ips_ctla4_neg_pd1_neg, fill=Cluster)) + 
  labs(x="", y="IPS-CTLA4(-)/PD-1(-)") +
  theme_bw(base_size = 9) + ###去除背景颜色
  #geom_text(aes(label=Hugo_Symbol),check_overlap = TRUE,size=6) +
  geom_boxplot(width=0.5) +
  # 不显示离群点
  #geom_jitter(width = 0.1) +
  theme(axis.text = element_text(size = 20))+ 
  #theme(axis.title.x=element_text(vjust=2, size=20,face = "bold"))+
  # 添加散点
  # geom_line(aes(group=Hugo_Symbol),
  #           color="black", alpha=1,linetype=2,
  #           linewidth=0.8) +
  theme(axis.title.y=element_text(vjust=2, size=20))+
  scale_fill_manual(values = c('#cab2d6','#ffff99'))+
  stat_compare_means( method = "wilcox.test",label.x.npc = "center",size=7)
p
ggsave("G:/免疫亚型分型/result5/TICA/CTLA4-_PD-1-.pdf",p,height = 5,width = 5)

p = ggplot(TICA_STAD, aes(x = Cluster, y = ips_ctla4_pos_pd1_neg, fill=Cluster)) + 
  labs(x="", y="IPS-CTLA4(+)/PD-1(-)") +
  theme_bw(base_size = 9) + ###去除背景颜色
  #geom_text(aes(label=Hugo_Symbol),check_overlap = TRUE,size=6) +
  geom_boxplot(width=0.5) +
  # 不显示离群点
  #geom_jitter(width = 0.1) +
  theme(axis.text = element_text(size = 20))+ 
  #theme(axis.title.x=element_text(vjust=2, size=20,face = "bold"))+
  # 添加散点
  # geom_line(aes(group=Hugo_Symbol),
  #           color="black", alpha=1,linetype=2,
  #           linewidth=0.8) +
  theme(axis.title.y=element_text(vjust=2, size=20))+
  scale_fill_manual(values = c('#cab2d6','#ffff99'))+
  stat_compare_means( method = "wilcox.test",label.x.npc = "center",size=7)
p
ggsave("G:/免疫亚型分型/result5/TICA/CTLA4+_PD-1-.pdf",p,height = 5,width = 5)
