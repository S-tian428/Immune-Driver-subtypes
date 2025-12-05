rm(list = ls())
setwd("I:/免疫亚型分型/MOVICS/")#输入
load("cmoic.brca.rda")
load("STAD.multiomics.rda")
surv.info<-STAD_multiomics$surv.info
#####样本临床信息
STAD_clinical<-fread("I:/免疫亚型分型/原始多组学数据/TCGA-STAD.GDC_phenotype.tsv.gz",header = T, sep = '\t',data.table = F)
index1<-match(rownames(surv.info),STAD_clinical$submitter_id.samples)
surv.info$age<-STAD_clinical$age_at_initial_pathologic_diagnosis[index1]
surv.info$age
clust.res<-cmoic.brca$clust.res
clust.res$clust<-paste0("CS",clust.res$clust)
index<-match(rownames(surv.info),clust.res$samID)
surv.info$Cluster<-clust.res$clust[index]
####组间年龄差异
my_comparisons <-list(c("CS1", "CS2"))
library(ggplot2)
library(ggpubr)
# 使用ggplot2绘制小提琴图
colnames(surv.info)
plot<-ggplot(surv.info, aes(x = Cluster, y = age, fill = Cluster)) +
  geom_violin(trim = FALSE) +  # 绘制小提琴图，trim = FALSE保留完整形状
  geom_boxplot(width = 0.1, fill = "white") +  # 在小提琴图内部绘制箱式图，设置宽度和填充色
  theme_minimal() + # 使用简洁的主题
 # labs(title = "distribution of age among different groups", x = "group", y = "age", fill = "Group") +  # 设置轴标签和图例标题
  scale_fill_brewer(palette = "Pastel1")  # 使用漂亮的颜色填充（可选）

# 添加显著性标记
p_with_signif <- plot + stat_compare_means(comparisons = my_comparisons, # 定义你想比较哪些组，例如list(c("Group1", "Group2"))
                                           # 显示显著性标记和p值
                                           method = "wilcox.test", # 使用t检验进行比较
                                           tip.length =0, # 显著性标记线的长度
                                          # size = 5
                                          ) # 显著性标记的大小
p_with_signif
ggsave("I:/免疫亚型分型/result2/Age.pdf",p_with_signif,width = 4.5,height = 4)

####stage_M
table(surv.info$stageM)
StageM_surv.info<-data.frame(stageM=surv.info$stageM,
                             Cluster=surv.info$Cluster)

StageM_surv.info<-StageM_surv.info[StageM_surv.info$stageM%in%
                                     c("M0","M1"),]
fisher_res<-table(StageM_surv.info)
fisher_test<-fisher.test(fisher_res)
fisher_test$p.value
p<-ggplot(data = StageM_surv.info) + 
  ggtitle(paste("Fisher’s exact test P=",fisher_test$p.value,sep = ""))+
  geom_bar(mapping = aes(x = Cluster, fill = stageM), position = "fill")+
  scale_fill_manual(values=c("#21498d","#e71d36"))+theme_bw()+theme(panel.grid=element_blank())
p
ggsave("I:/免疫亚型分型/result2/StageM.pdf",p,width = 5.5,height = 6)

####stage_N
table(surv.info$stageN)
surv.info$stageN<-gsub("N3.*","N3",surv.info$stageN)

StageN_surv.info<-data.frame(stageN=surv.info$stageN,
                             Cluster=surv.info$Cluster)

StageN_surv.info<-StageN_surv.info[StageN_surv.info$stageN%in%
                                     c("N0","N1","N2","N3"),]
fisher_res<-table(StageN_surv.info)
fisher_test<-fisher.test(fisher_res)
fisher_test$p.value
p<-ggplot(data = StageN_surv.info) + 
  ggtitle(paste("Fisher’s exact test P=",fisher_test$p.value,sep = ""))+
  geom_bar(mapping = aes(x = Cluster, fill = stageN), position = "fill")+
  scale_fill_manual(values=c("#21498d","#b3cde3",'#F29C71',"#e71d36"))+theme_bw()+theme(panel.grid=element_blank())
p
ggsave("I:/免疫亚型分型/result2/StageN.pdf",p,width = 5.5,height = 6)
####stage_T
table(surv.info$stageT)
surv.info$stageT<-gsub("T1.*","T1",surv.info$stageT)
surv.info$stageT<-gsub("T2.*","T2",surv.info$stageT)
surv.info$stageT<-gsub("T4.*","T4",surv.info$stageT)
StageT_surv.info<-data.frame(stageT=surv.info$stageT,
                             Cluster=surv.info$Cluster)
fisher_res<-table(StageT_surv.info)
fisher_test<-fisher.test(fisher_res)
fisher_test$p.value
p<-ggplot(data = StageT_surv.info) + 
  ggtitle(paste("Fisher’s exact test P=",fisher_test$p.value,sep = ""))+
  geom_bar(mapping = aes(x = Cluster, fill = stageT), position = "fill")+
  scale_fill_manual(values=c("#21498d","#b3cde3",'#F29C71',"#e71d36"))+theme_bw()+theme(panel.grid=element_blank())
p
ggsave("I:/免疫亚型分型/result2/StageT.pdf",p,width = 5.5,height = 6)

####grade
table(surv.info$grade)
Grade_surv.info<-data.frame(grade=surv.info$grade,
                             Cluster=surv.info$Cluster)

Grade_surv.info<-Grade_surv.info[Grade_surv.info$grade%in%
                                     c("G1","G2","G3"),]
fisher_res<-table(Grade_surv.info)
fisher_res
fisher_test<-fisher.test(fisher_res)
fisher_test$p.value
p<-ggplot(data = Grade_surv.info) + 
  ggtitle(paste("Fisher’s exact test P=",fisher_test$p.value,sep = ""))+
  geom_bar(mapping = aes(x = Cluster, fill = grade), position = "fill")+
  scale_fill_manual(values=c("#21498d",'#F29C71',"#e71d36"))+theme_bw()+theme(panel.grid=element_blank())
p
ggsave("I:/免疫亚型分型/result2/Grade.pdf",p,width = 5.5,height = 6)

####stage
table(surv.info$stage)
surv.info$stage<-gsub("iv","IV",surv.info$stage)
surv.info$stage<-gsub("iii.*","III",surv.info$stage)
surv.info$stage<-gsub("ii.*","II",surv.info$stage)
surv.info$stage<-gsub("i.*","I",surv.info$stage)
Stage_surv.info<-data.frame(stage=surv.info$stage,
                            Cluster=surv.info$Cluster)
table(Stage_surv.info$stage)
Stage_surv.info<-Stage_surv.info[Stage_surv.info$stage%in%
                                   c("stage I","stage II",
                                     "stage III","stage IV"),]
fisher_res<-table(Stage_surv.info)
fisher_res
fisher_test<-fisher.test(fisher_res)
fisher_test$p.value
p<-ggplot(data = Stage_surv.info) + 
  ggtitle(paste("Fisher’s exact test P=",fisher_test$p.value,sep = ""))+
  geom_bar(mapping = aes(x = Cluster, fill = stage), position = "fill")+
  scale_fill_manual(values=c("#21498d","#b3cde3",'#F29C71',"#e71d36"))+theme_bw()+theme(panel.grid=element_blank())
p
ggsave("I:/免疫亚型分型/result2/Stage.pdf",p,width = 5.5,height = 6)

####response
table(surv.info$response)
surv.info$response<-gsub("Complete Remission/Response","Complete Response",surv.info$response)
surv.info$response<-gsub("Partial Remission/Response","Partial Response",surv.info$response)

Response_surv.info<-data.frame(response=surv.info$response,
                            Cluster=surv.info$Cluster)
Response_surv.info<-Response_surv.info[Response_surv.info$response%in%
                                    c("Complete Response","Partial Response",
                                    "Progressive Disease","Stable Disease"),]
fisher_res<-table(Response_surv.info)
fisher_res<-t(fisher_res)
fisher_test<-fisher.test(fisher_res)
fisher_test$p.value
p<-ggplot(data = Response_surv.info) + 
  ggtitle(paste("Fisher’s exact test P=",fisher_test$p.value,sep = ""))+
  geom_bar(mapping = aes(x = Cluster, fill = response), position = "fill")+
  scale_fill_manual(values=c("#21498d","#b3cde3",'#F29C71',"#e71d36"))+theme_bw()+theme(panel.grid=element_blank())
p
ggsave("I:/免疫亚型分型/result2/Response.pdf",p,width = 5.5,height = 6)
###Gender
table(surv.info$gender)
#汇总M分期信息为表格
Gender_surv.info<-data.frame(gender=surv.info$gender,
                             Cluster=surv.info$Cluster)

fisher_res<-table(Gender_surv.info)
fisher_test<-fisher.test(fisher_res)
fisher_test$p.value
p<-ggplot(data = Gender_surv.info) + 
  ggtitle(paste("Fisher’s exact test P=",fisher_test$p.value,sep = ""))+
  geom_bar(mapping = aes(x = Cluster, fill = gender), position = "fill")+
  scale_fill_manual(values=c("#21498d","#e71d36"))+theme_bw()+theme(panel.grid=element_blank())
p

ggsave("I:/免疫亚型分型/result2/Gender.pdf",p,width = 5.5,height = 6)

STAD_multiomics$surv.info<-surv.info
save(STAD_multiomics,file = "STAD.multiomics.rda")
getwd()
