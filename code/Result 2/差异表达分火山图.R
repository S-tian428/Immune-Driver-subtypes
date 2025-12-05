rm(list=ls())
setwd("H:/免疫亚型分型/result5")
####绘制火山图
#####火山图
library(tibble)
library(ggrepel)
library(dplyr)
limma_result<-read.table("limma_result.txt",header = T,sep="\t",stringsAsFactors=FALSE,quote = "")
colnames(limma_result)
limma_result<-na.omit(limma_result)
data <- 
  limma_result %>% 
  mutate(change = as.factor(ifelse(adj.P.Val < 0.05 & abs(logFC) > 1,
                                   ifelse(logFC > 1 ,'Up','Down'),'No Change'))) %>% 
  rownames_to_column('gene')
colnames(data)
# 普通火山图
p <- ggplot(data = data, 
            aes(x = logFC, 
                y = -log10(P.Value))) +  # 设置x轴为logFC，y轴为-P.Value的对数值
  geom_point(alpha = 0.4, size = 3.5, 
             aes(color = change)) +      # 添加散点，根据change列着色
  ylab("-log10(Pvalue)") +               # 设置y轴标签
  scale_color_manual(values = c("#128acc", "grey", "#f2672a")) +  # 设置颜色映射，蓝色表示下调，灰色表示稳定，红色表示上调
  geom_vline(xintercept = c(-1,1), lty = 4, col = "black", lwd = 0.8) +  # 添加垂直参考线，用于标记logFC阈值
  geom_hline(yintercept = -log10(0.05), lty = 4, col = "black", lwd = 0.8) +   # 添加水平参考线，用于标记-P.Value阈值
  theme_bw()  # 使用网格白底主题
p
# 标签添加

# 按需求添加标签

# 创建一个新的列label，并初始化为NA
data$label <- ""  
#data<-na.omit(data)
# 根据symbol的值，为特定基因添加标签信息
data<-data[order(data$logFC,decreasing = T),]
data$label[(nrow(data)-9):nrow(data)] <- data$gene[(nrow(data)-9):nrow(data)]
data$label[1:10] <- data$gene[1:10]

p1<-p +  # 基于普通火山图p
  geom_text_repel(data = data, aes(x = logFC, 
                                   y = -log10(P.Value), label = label),
                  max.overlaps = 100)
p1
ggsave("H:/免疫亚型分型/result2/Diffgene_volcano.pdf",p1,width = 7,height = 5)

# # 在普通火山图p的基础上，添加标签，并使用geom_label_repel函数进行标签的绘制
# p0 <- p + geom_label_repel(data = data, aes(label = label),
#                            size = 4,                           # 设置标签大小
#                            box.padding = unit(0.5, "lines"),   # 设置标签内边距
#                            point.padding = unit(0.8, "lines"), # 设置标签与点的距离
#                            segment.color = "black",            # 设置标签边界线颜色
#                            show.legend = FALSE,                # 不显示图例
#                            max.overlaps = 10000)               # 设置标签重叠的最大次数
# p0
