rm(list = ls())
#BiocManager::install("BSgenome.Hsapiens.UCSC.hg38")
library(BSgenome.Hsapiens.UCSC.hg38)
df <- data.frame(chromName = seqnames(BSgenome.Hsapiens.UCSC.hg38), # 染色体名字
                 chromlength = seqlengths(BSgenome.Hsapiens.UCSC.hg38)# 每条染色体长度
)
df$chromNum <- 1:length(df$chromName) # 染色体名字纯数字版，为了和scores.gistic里面的对应

# 我们用的原始CNV文件是没有性染色体信息的
df <- df[1:22,] # 只要前22条染色体信息
df$chromlengthCumsum <- cumsum(as.numeric(df$chromlength)) # 染色体累加长度
# 得到每条染色体从0开始的起始坐标
df$chormStartPosFrom0 <- c(0,df$chromlengthCumsum[-nrow(df)])
# 计算每条染色体中间位置坐标，用来最后加文字
tmp_middle <- diff(c(0,df$chromlengthCumsum)) / 2
df$chromMidelePosFrom0 <- df$chormStartPosFrom0 + tmp_middle
# 如果你不知道用哪个函数读取，多试几次就知道了！
####CS1
scores <- read.table("G:/免疫亚型分型/result2/Gistic2.0/CS1_results/scores.gistic",sep="\t",header=T,stringsAsFactors = F)
chromID <- scores$Chromosome
scores$StartPos <- scores$Start + df$chormStartPosFrom0[chromID]
scores$EndPos <- scores$End + df$chormStartPosFrom0[chromID]
###把Del的G.score变成负数
range(scores$G.score)
scores[scores$Type == "Del", "G.score"] <- scores[scores$Type == "Del", "G.score"] * -1
range(scores$G.score)
#install.packages("ggplot2")
library(ggplot2)
library(ggsci)
df$ypos <- rep(c(0.2,0.25),11)
CS1_plot<-ggplot(scores, aes(StartPos, G.score))+
  geom_area(aes(group=Type, fill=factor(Type,levels = c("Del","Amp"))))+
  scale_fill_lancet(guide=guide_legend(reverse = T),name="Type")+
  geom_vline(data = df ,mapping=aes(xintercept=chromlengthCumsum),linetype=2)+
  geom_text(data = df,aes(x=chromMidelePosFrom0,y=ypos,label=chromName))+
  scale_x_continuous(expand = c(0,-1000),limits = c(0,2.9e9),name = NULL,labels = NULL)+
  ylim(-0.3,0.3)+
  theme_minimal()+
  theme(legend.position = "top",
        axis.text.y = element_text(color = "black",size = 14),
        axis.title.y = element_text(color = "black",size = 16)
  )
ggsave("G:/免疫亚型分型/result2/Gistic2.0/Gscore_CS1.pdf",CS1_plot,width = 10,height = 5.5)

####frequency
scores[scores$Type == "Del", "frequency"] <- scores[scores$Type == "Del", "frequency"] * -1
CS1_frequency_plot<-ggplot(scores, aes(StartPos, frequency))+
  geom_area(aes(group=Type, fill=factor(Type,levels = c("Del","Amp"))))+
  scale_fill_lancet(guide=guide_legend(reverse = T),name="Type")+
  geom_vline(data = df ,mapping=aes(xintercept=chromlengthCumsum),linetype=2)+
  geom_text(data = df,aes(x=chromMidelePosFrom0,y=ypos,label=chromName))+
  scale_x_continuous(expand = c(0,-1000),limits = c(0,2.9e9),name = NULL,labels = NULL)+
  ylim(-0.3,0.3)+
  theme_minimal()+
  theme(legend.position = "top",
        axis.text.y = element_text(color = "black",size = 14),
        axis.title.y = element_text(color = "black",size = 16)
  )
ggsave("G:/免疫亚型分型/result2/Gistic2.0/frequency_CS1.pdf",CS1_frequency_plot,width = 11,height = 5)

####CS2
scores <- read.table("G:/免疫亚型分型/result2/Gistic2.0/CS2_results/scores.gistic",sep="\t",header=T,stringsAsFactors = F)
chromID <- scores$Chromosome
scores$StartPos <- scores$Start + df$chormStartPosFrom0[chromID]
scores$EndPos <- scores$End + df$chormStartPosFrom0[chromID]
###把Del的G.score变成负数
range(scores$G.score)
scores[scores$Type == "Del", "G.score"] <- scores[scores$Type == "Del", "G.score"] * -1
range(scores$G.score)
#install.packages("ggplot2")
library(ggplot2)
library(ggsci)
df$ypos <- rep(c(0.2,0.25),11)
CS2_plot<-ggplot(scores, aes(StartPos, G.score))+
  geom_area(aes(group=Type, fill=factor(Type,levels = c("Del","Amp"))))+
  scale_fill_lancet(guide=guide_legend(reverse = T),name="Type")+
  geom_vline(data = df ,mapping=aes(xintercept=chromlengthCumsum),linetype=2)+
  geom_text(data = df,aes(x=chromMidelePosFrom0,y=ypos,label=chromName))+
  scale_x_continuous(expand = c(0,-1000),limits = c(0,2.9e9),name = NULL,labels = NULL)+
  ylim(-0.3,0.3)+
  theme_minimal()+
  theme(legend.position = "top",
        axis.text.y = element_text(color = "black",size = 14),
        axis.title.y = element_text(color = "black",size = 16)
  )
ggsave("G:/免疫亚型分型/result2/Gistic2.0/Gscore_CS2.pdf",CS2_plot,width = 10,height = 5.5)

####frequency
scores[scores$Type == "Del", "frequency"] <- scores[scores$Type == "Del", "frequency"] * -1
CS2_frequency_plot<-ggplot(scores, aes(StartPos, frequency))+
  geom_area(aes(group=Type, fill=factor(Type,levels = c("Del","Amp"))))+
  scale_fill_lancet(guide=guide_legend(reverse = T),name="Type")+
  geom_vline(data = df ,mapping=aes(xintercept=chromlengthCumsum),linetype=2)+
  geom_text(data = df,aes(x=chromMidelePosFrom0,y=ypos,label=chromName))+
  scale_x_continuous(expand = c(0,-1000),limits = c(0,2.9e9),name = NULL,labels = NULL)+
  ylim(-0.3,0.3)+
  theme_minimal()+
  theme(legend.position = "top",
        axis.text.y = element_text(color = "black",size = 14),
        axis.title.y = element_text(color = "black",size = 16)
  )
ggsave("G:/免疫亚型分型/result2/Gistic2.0/frequency_CS2.pdf",CS2_frequency_plot,width = 11,height = 5)
