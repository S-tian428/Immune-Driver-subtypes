##免疫细胞浸润
#####Cibersort
rm(list=ls())
library(tibble)
library(dplyr)
library(tidyr)
load("I:/免疫亚型分型/MOVICS/STAD.multiomics.rda")
setwd("I:/免疫亚型分型/result3/MCPcounter/")
surv.info<-STAD_multiomics$surv.info
library(devtools)
#install_github("ebecht/MCPcounter",ref="master", subdir="Source")
library(MCPcounter)
probesets <- data.table::fread("probesets.txt",data.table = F,header = F)
genes=read.table("genes.txt",
                 sep="\t",stringsAsFactors=FALSE,header=TRUE,colClasses="character",check.names=FALSE)
probesets=read.table("probesets.txt",
                     sep="\t",stringsAsFactors=FALSE,header=F,colClasses="character",check.names=FALSE)
STAD.FPKM=STAD_multiomics$STAD.FPKM
results<- MCPcounter.estimate(STAD.FPKM,featuresType=c("HUGO_symbols")[1],
                              probesets=probesets,
                              genes=genes)
results<-t(results)
results<-as.data.frame(results)
results$Cluster<-surv.info$Cluster
#画图
library(ggplot2)
library(ggpubr)
library(ggsci)
library(lemon)
#设置分组信息
###宽数据变长数据
library(tidyr)
library(tidyr)
results <- gather(results, key = "ImmCells", value = "value",
                    -"Cluster")
p<-ggplot(results, aes(x = Cluster, y = value, fill = Cluster)) +
  geom_boxplot(position = position_dodge(0.8)) +
  theme_bw(base_size = 9) + ###去除背景颜色
  theme(panel.grid=element_blank()) +
  facet_wrap(. ~ ImmCells,scales = 'free',ncol = 5)+
  scale_fill_nejm() +
  stat_compare_means( method = "wilcox.test",label = "p.signif",label.x.npc = "center",size=2.3)
ggsave("MCPcounter.pdf",p,width = 7,height = 4)
