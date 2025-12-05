rm(list=ls())
####TIDE输入前需要对表达谱进行归一化处理
load("H:/免疫亚型分型/MOVICS/STAD.multiomics.rda")
STAD_fpkm<-STAD_multiomics$STAD.FPKM
surv.info<-STAD_multiomics$surv.info
normalized_data <-  t(apply(STAD_fpkm, 1, function(x)x-(mean(x))))
setwd("H:/免疫亚型分型/result5/TIDE/")
write.table(data.frame(ID=rownames(normalized_data),normalized_data,
                       check.names = F),"TIDE_input.txt",
            sep="\t",quote=FALSE,row.names = T)
TIDE_STAD <- read.csv(file="H:/免疫亚型分型/result5/TIDE/TIDE_RESULT.csv", header=TRUE, sep=",")
index<-match(TIDE_STAD$Patient,rownames(surv.info))
TIDE_STAD$Cluster<-surv.info$Cluster[index]
ggplot(TIDE_STAD,aes(x=Patient,TIDE))+
  geom_col(aes(fill=Cluster))+
  theme(axis.text.x = element_text(angle= 90 , vjust= .5 , hjust= 1 ))+ 
  scale_fill_manual(values=c("#00c16e","#7552cc"))#自定义颜色
library(ggplot2)
library(ggsci)
library(ggpubr)
p = ggplot(TIDE_STAD, aes(x = Cluster, y = Dysfunction, color=Cluster)) + 
  labs(x="", y="Mutation Number") +
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
  scale_color_manual(values = c('orangered3','#1f78b4'))+
  stat_compare_means( method = "wilcox.test",label = "p.signif",label.x.npc = "center",size=7)
p
