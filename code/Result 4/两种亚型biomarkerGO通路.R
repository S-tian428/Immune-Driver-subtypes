#############GO 富集分析
#####模块1
rm(list=ls())
#library(stringr)
#library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
load("H:/免疫亚型分型/MOVICS/marker.up.rda")
templates<-marker.up$templates
CS1_biomarker<-templates$probe[templates$class%in%"CS1"]
CS2_biomarker<-templates$probe[templates$class%in%"CS2"]
#dat <- dat %>% dplyr::group_by(group_type) %>% dplyr::do(head(., n = 13)) 
######cytoscape

######基因名和ID转换
CS1_ids = bitr(CS1_biomarker, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#sce.markers = merge(sce.markers, ids, by.x='V1', by.y='SYMBOL')
#####GO分析 BP
CS1_go<-enrichGO(CS1_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                qvalueCutoff = 0.05,keyType = 'ENTREZID')
#使用simplify 对GO富集分析结果进行精简
CS1_go<- simplify(CS1_go)
CS1_go<-CS1_go@result
CS1_go<-CS1_go[order(CS1_go$Count,decreasing = T),]
CS1_go<-CS1_go[CS1_go$ONTOLOGY%in%"BP",]

write.table(CS1_go,"H:/免疫亚型分型/result5/CS1_biomarker_GO.txt",sep="\t",quote=FALSE,row.names = F)
####可视化
#######BP
####取前10个通路画图
CS1_go<-CS1_go[1:10,]
####提取有用的两列
CS1_go<-data.frame(Description=CS1_go$Description,geneID=CS1_go$geneID)
####数据转换，将一行分成多行
library(tidyr)
CS1_go <-CS1_go %>% as_tibble() %>% 
  separate_rows(geneID, sep = "/")
i=1
####id和基因名转换
CS1_go$gene_symbol<-NA
for (i in 1:nrow(CS1_go)) {
  
  index<-which(CS1_ids$ENTREZID%in%CS1_go$geneID[i])
  if(length(index)>0){
    CS1_go$gene_symbol[i]<-CS1_ids$SYMBOL[index]
  } 
  
}
#####长数据变宽数据
CS1_wide<-data.frame(spread(CS1_go,gene_symbol,geneID))
rownames(CS1_wide)<-CS1_wide$Description
CS1_wide<-CS1_wide[,-1]
####若基因富集到通路则赋值1，否则赋值0
CS1_wide[!is.na(CS1_wide)]<-1
CS1_wide[is.na(CS1_wide)]<-0
####将字符型变成数值型
a<-apply(CS1_wide,2,as.numeric)
rownames(a)<-rownames(CS1_wide)
CS1_wide<-t(a)
b<-apply(CS1_wide,1,sum)
CS1_wide<-as.data.frame(CS1_wide)
CS1_wide$logFC<-b
####加入差异表达的FC值
#expression_result<-read.table("F:/肝损伤再灌注/差异表达分析/mRNA/Sig_mRNA_result.txt",
#                              header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
#index1<-match(rownames(go_result_wide),expression_result$genes)
#go_result_wide$logFC<-expression_result$logFC[index1]
#go_result_wide<-as.data.frame(a)
#go_result_wide$logFC<-0
#str(go_result_wide)
CS1_wide<-as.matrix(CS1_wide)
####画图
#install.packages("GOplot")
library(GOplot)
CS1_p<-GOChord(CS1_wide, 
           space = 0.02, #GO term处间隔大小设置
           #limit = c(3, 5),#第一个数值为至少分配给一个基因的Goterm数，第二个数值为至少分配给一个GOterm的基因数
           gene.order = 'logFC', lfc.max=max(b),####lfc的最大值默认只有3，需要自己设置
           gene.space = 0.25, gene.size = 5,#基因排序，间隔，名字大小设置
           #lfc.col=c('firebrick3', 'white','royalblue3'),
           lfc.col=c('#FF8C00','white','#8B4500'),##上调下调颜色设置
           #ribbon.col=colorRampPalette(c("blue", "red"))(length(EC$process)),#GO term 颜色设置
           ribbon.col=brewer.pal(10, "Set3")#GO term 颜色设置
)
CS1_p
#go_result_wide<-as.matrix(go_result_wide)
ggsave("H:/免疫亚型分型/result5/CS1_GOChord.pdf",CS1_p,width=10,height=11)


######基因名和ID转换
CS2_ids = bitr(CS2_biomarker, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
#sce.markers = merge(sce.markers, ids, by.x='V1', by.y='SYMBOL')
#####GO分析 BP
CS2_go<-enrichGO(CS2_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont='ALL',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.05,keyType = 'ENTREZID')
#使用simplify 对GO富集分析结果进行精简
CS2_go<- simplify(CS2_go)
CS2_go<-CS2_go@result
CS2_go<-CS2_go[order(CS2_go$Count,decreasing = T),]
CS2_go<-CS2_go[CS2_go$ONTOLOGY%in%"BP",]
write.table(CS2_go,"H:/免疫亚型分型/result5/CS2_biomarker_GO.txt",sep="\t",quote=FALSE,row.names = F)
####可视化
#######BP
####取前10个通路画图
CS2_go<-CS2_go[1:10,]
####提取有用的两列
CS2_go<-data.frame(Description=CS2_go$Description,geneID=CS2_go$geneID)
####数据转换，将一行分成多行
library(tidyr)
CS2_go <-CS2_go %>% as_tibble() %>% 
  separate_rows(geneID, sep = "/")
i=1
####id和基因名转换
CS2_go$gene_symbol<-NA
for (i in 1:nrow(CS2_go)) {
  
  index<-which(CS2_ids$ENTREZID%in%CS2_go$geneID[i])
  if(length(index)>0){
    CS2_go$gene_symbol[i]<-CS2_ids$SYMBOL[index]
  } 
  
}
#####长数据变宽数据
CS2_wide<-data.frame(spread(CS2_go,gene_symbol,geneID))
rownames(CS2_wide)<-CS2_wide$Description
CS2_wide<-CS2_wide[,-1]
####若基因富集到通路则赋值1，否则赋值0
CS2_wide[!is.na(CS2_wide)]<-1
CS2_wide[is.na(CS2_wide)]<-0
####将字符型变成数值型
a<-apply(CS2_wide,2,as.numeric)
rownames(a)<-rownames(CS2_wide)
CS2_wide<-t(a)
b<-apply(CS2_wide,1,sum)
CS2_wide<-as.data.frame(CS2_wide)
CS2_wide$logFC<-b
####加入差异表达的FC值
#expression_result<-read.table("F:/肝损伤再灌注/差异表达分析/mRNA/Sig_mRNA_result.txt",
#                              header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
#index1<-match(rownames(go_result_wide),expression_result$genes)
#go_result_wide$logFC<-expression_result$logFC[index1]
#go_result_wide<-as.data.frame(a)
#go_result_wide$logFC<-0
#str(go_result_wide)
CS2_wide<-as.matrix(CS2_wide)
####画图
#install.packages("GOplot")
library(GOplot)
CS2_p<-GOChord(CS2_wide, 
               space = 0.02, #GO term处间隔大小设置
               #limit = c(3, 5),#第一个数值为至少分配给一个基因的Goterm数，第二个数值为至少分配给一个GOterm的基因数
               gene.order = 'logFC', lfc.max=max(b),####lfc的最大值默认只有3，需要自己设置
               gene.space = 0.25, gene.size = 5,#基因排序，间隔，名字大小设置
               #lfc.col=c('firebrick3', 'white','royalblue3'),
               lfc.col=c('#FF8C00','white','#8B4500'),##上调下调颜色设置
               #ribbon.col=colorRampPalette(c("blue", "red"))(length(EC$process)),#GO term 颜色设置
               ribbon.col=brewer.pal(10, "Set3")#GO term 颜色设置
)
CS2_p
#go_result_wide<-as.matrix(go_result_wide)
ggsave("H:/免疫亚型分型/result5/CS2_GOChord.pdf",CS2_p,width=10,height=11)
