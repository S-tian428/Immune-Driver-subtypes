#############GO 富集分析
#####模块1
rm(list=ls())
#library(stringr)
#library(enrichplot)
library(clusterProfiler)
library(org.Hs.eg.db)
load("G:/免疫亚型分型/MOVICS/marker.up.rda")
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
write.table(CS1_go,"G:/免疫亚型分型/result5/CS1_biomarker_GO.txt",sep="\t",quote=FALSE,row.names = F)
CS1_go_plot<-CS1_go[1:10,]
CS1_go_plot$subtypes<-"CS1"
CS1_go_plot$Count=-CS1_go_plot$Count
#CS1_go<-CS1_go[CS1_go$ONTOLOGY%in%"BP",]
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
write.table(CS2_go,"G:/免疫亚型分型/result5/CS2_biomarker_GO.txt",sep="\t",quote=FALSE,row.names = F)
CS2_go_plot<-CS2_go[1:10,]
CS2_go_plot$subtypes<-"CS2"
go_PLOT<-rbind(CS1_go_plot,CS2_go_plot)
go_PLOT$Description[15]<-"adaptive immune response based on ..."
colnames(go_PLOT)
go_PLOT$Count
library(ggplot2)
p<-ggplot(go_PLOT, aes(x = reorder(Description,Count), Count,fill=subtypes)) + 
  geom_bar(stat = 'identity',alpha = 0.7) + 
  coord_flip() + 
  theme_bw() + #去除背景色
  labs(x = "")+
  scale_fill_manual(values = c("#0072b5","#f79c4a"))#设置颜色
ggsave("G:/免疫亚型分型/result6/biomarker_GO.pdf",p,height = 3,width = 6)
