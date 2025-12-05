rm(list=ls())
setwd("H:/免疫亚型分型/result5")
limma_result<-read.table("limma_result.txt",header = T,sep="\t",stringsAsFactors=FALSE,quote = "")
CS1_up<-rownames(limma_result)[limma_result$logFC<(-1) & limma_result$adj.P.Val<0.05]
CS2_up<-rownames(limma_result)[limma_result$logFC>1 & limma_result$adj.P.Val<0.05]
####分别对上下调基因进行GO和KEGG富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
####CS1_up基因
######基因名和ID转换
CS1_ids = bitr(CS1_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#sce.markers = merge(sce.markers, ids, by.x='V1', by.y='SYMBOL')
#####GO分析 
CS1_go<-enrichGO(CS1_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                qvalueCutoff = 0.05,keyType = 'ENTREZID')
#使用simplify 对GO富集分析结果进行精简
CS1_go<- simplify(CS1_go)
CS1_result<-CS1_go@result
# CS1_result<-CS1_result[order(CS1_result$Count,decreasing = T),]
###BP取前十个
# CS1_BP<-CS1_result[CS1_result$ONTOLOGY%in%"BP",]
CS1_result<-CS1_result[1:10,]
# CS1_MF<-CS1_result[CS1_result$ONTOLOGY%in%"MF",]
# CS1_CC<-CS1_result[CS1_result$ONTOLOGY%in%"CC",]
# CS1_ggplot<-rbind(CS1_BP,CS1_MF,CS1_CC)
CS1_result$log10adj=-log10(CS1_result$p.adjust)
colnames(CS1_result)
  #横向柱状图#
library(scales)
# scale_colour_gradient2(low = muted("red"),
#                        
#                        mid = "white",
#                        
#                        high = muted("blue"),
#                        
#                        midpoint = 110)
#纵向柱状图#
#横向柱状图#
library(RColorBrewer)
cols <- rev(brewer.pal(11, 'YlOrRd'))
cols
[1] "#313695" "#4575B4" "#74ADD1" "#ABD9E9" "#E0F3F8" "#FFFFBF"
[7] "#FEE090" "#FDAE61" "#F46D43" "#D73027" "#A50026"
[1] "#800026" "#BD0026" "#E31A1C" "#FC4E2A" "#FD8D3C" "#FEB24C"
[7] "#FED976" "#FFEDA0" "#FFFFCC"
CS1_GOplot<-ggplot(CS1_result, 
       aes(x=Description,y=Count)) +  #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(aes(fill=log10adj),stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_gradient(low ="#ABD9E9",
                        high ="#FEB24C")+ #柱状图填充颜色
  #facet_grid(.~ONTOLOGY, scale = 'free_x', space = 'free_x')+
  xlab("GO term") + #x轴标签
  ylab("Gene Number") +  #y轴标签
  labs(title = "CS1 GO Terms")+ #设置标题
  theme_bw() + 
  theme(axis.text.x=element_text(family="sans",face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 )) #对字体样式、颜色、还有横坐标角度（）

CS1_GOplot
#help(theme) #查阅这个函数其他具体格式
ggsave("H:/免疫亚型分型/result2/CS1_GOplot.pdf",CS1_GOplot,width = 7,height=8)

####CS2_up基因
######基因名和ID转换
CS2_ids = bitr(CS2_up, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#sce.markers = merge(sce.markers, ids, by.x='V1', by.y='SYMBOL')
#####GO分析 
CS2_go<-enrichGO(CS2_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.05,keyType = 'ENTREZID')
#使用simplify 对GO富集分析结果进行精简
CS2_go<- simplify(CS2_go)
CS2_result<-CS2_go@result

CS2_result<-CS2_result[1:10,]
# CS1_MF<-CS1_result[CS1_result$ONTOLOGY%in%"MF",]
# CS1_CC<-CS1_result[CS1_result$ONTOLOGY%in%"CC",]
# CS1_ggplot<-rbind(CS1_BP,CS1_MF,CS1_CC)
CS2_result$log10adj=-log10(CS2_result$p.adjust)
colnames(CS2_result)
#横向柱状图#
library(scales)
# scale_colour_gradient2(low = muted("red"),
#                        
#                        mid = "white",
#                        
#                        high = muted("blue"),
#                        
#                        midpoint = 110)
#纵向柱状图#
#横向柱状图#
library(RColorBrewer)
cols <- rev(brewer.pal(11, 'RdYlBu'))
cols
[1] "#313695" "#4575B4" "#74ADD1" "#ABD9E9" "#E0F3F8" "#FFFFBF"
[7] "#FEE090" "#FDAE61" "#F46D43" "#D73027" "#A50026"
CS2_GOplot<-ggplot(CS2_result, 
                   aes(x=Description,y=Count)) +  #x、y轴定义；根据ONTOLOGY填充颜色
  geom_bar(aes(fill=log10adj),stat="identity", width=0.8) +  #柱状图宽度
  scale_fill_gradient(low ="#ABD9E9",
                      high ="#FEB24C")+ #柱状图填充颜色
  #facet_grid(.~ONTOLOGY, scale = 'free_x', space = 'free_x')+
  xlab("GO term") + #x轴标签
  ylab("Gene Number") +  #y轴标签
  labs(title = "CS2 GO Terms")+ #设置标题
  theme_bw() + 
  theme(axis.text.x=element_text(family="sans",face = "bold", color="gray50",angle = 70,vjust = 1, hjust = 1 )) #对字体样式、颜色、还有横坐标角度（）

CS2_GOplot
#help(theme) #查阅这个函数其他具体格式
ggsave("H:/免疫亚型分型/result2/CS2_GOplot.pdf",CS2_GOplot,width = 10,height=12)




