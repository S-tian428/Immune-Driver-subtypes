rm(list = ls())
load("H:/免疫亚型分型/MOVICS/STAD.multiomics.rda")
surv.info<-STAD_multiomics$surv.info
STAD.mut<-STAD_multiomics$STAD.maf
load("H:/免疫亚型分型/MOVICS/cmoic.brca.rda")
index<-match(rownames(surv.info),cmoic.brca$clust.res$samID)
surv.info$Cluster<-cmoic.brca$clust.res$clust[index]
surv.info$Cluster<-paste0("CS",surv.info$Cluster)
STAD_multiomics$surv.info<-surv.info
save(STAD_multiomics,file = "H:/免疫亚型分型/MOVICS/STAD.multiomics.rda")
table(STAD.mut$Variant_Type)
####提取突变频率最高的Top基因做柱形图比较
mutations<-data.frame(Hugo_Symbol=STAD.mut$Hugo_Symbol,
                      Tumor_Sample_Barcode=STAD.mut$Tumor_Sample_Barcode)
#####计算基因在每个样本中的突变频率
#mutations<-unique(mutations)
####确定样本亚型
index<-match(mutations$Tumor_Sample_Barcode,rownames(surv.info))
mutations$Cluster<-surv.info$Cluster[index]
mutations$num<-1
CS1_mutation<-mutations[mutations$Cluster%in%"CS1",]
CS1_mutation<-unique(CS1_mutation)
CS1_mut<- aggregate(CS1_mutation$num,by= list(CS1_mutation$Hugo_Symbol,
                      CS1_mutation$Cluster), sum)
colnames(CS1_mut)<-c("Hugo_Symbol","Cluster","nMut")
CS1_mut$freq<-CS1_mut$nMut/length(rownames(surv.info)[surv.info$Cluster%in%"CS1"])

CS2_mutation<-mutations[mutations$Cluster%in%"CS2",]
CS2_mutation<-unique(CS2_mutation)
CS2_mut<- aggregate(CS2_mutation$num,by= list(CS2_mutation$Hugo_Symbol,
                                              CS2_mutation$Cluster), sum)
colnames(CS2_mut)<-c("Hugo_Symbol","Cluster","nMut")
CS2_mut$freq<-CS2_mut$nMut/length(rownames(surv.info)[surv.info$Cluster%in%"CS2"])
mutation_result<-rbind(CS1_mut,CS2_mut)
Sig_gene<-read.table("H:/免疫亚型分型/result2/INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION.txt",
                     header = T,sep="\t",stringsAsFactors=FALSE,quote = "")
Sig_gene<-Sig_gene[Sig_gene$padj<0.05,]
library(stringr)
Sig_gene$nMUT<-as.numeric(str_split_fixed(Sig_gene$TMB, " ", 2)[,1])
Sig_gene<-Sig_gene[order(Sig_gene$nMUT,decreasing = T),]
Sig_gene<-Sig_gene$Gene..Mutated.
mutation_result<-mutation_result[mutation_result$Hugo_Symbol%in%Sig_gene,]
rownames(mutation_result)<-NULL
# mutation_result$nMut<-1
# mutation_result1<- aggregate(mutation_result$nMut,by= list(mutation_result$Hugo_Symbol,
#     mutation_result$subtype), sum)
colnames(mutation_result)
table(mutation_result$Hugo_Symbol)
library(ggplot2)
library(ggsci)
library(ggpubr)
p = ggplot(mutation_result, aes(x = Cluster, y = freq, color=Cluster)) + 
  labs(x="", y="Mutation Number") +
  theme_bw(base_size = 9) + ###去除背景颜色
  geom_text(aes(label=Hugo_Symbol),check_overlap = TRUE,size=6) +
  geom_boxplot(width=0.5, outlier.colour = NA) +
  # 不显示离群点
 geom_jitter(width = 0.1) +
 theme(axis.text = element_text(size = 20))+ 
 #theme(axis.title.x=element_text(vjust=2, size=20,face = "bold"))+
  # 添加散点
  geom_line(aes(group=Hugo_Symbol),
            color="black", alpha=1,linetype=2,
            linewidth=0.8) +
  theme(axis.title.y=element_text(vjust=2, size=20))+
  scale_color_manual(values = c('orangered3','#1f78b4'))+
  stat_compare_means( method = "wilcox.test",label = "p.signif",label.x.npc = "center",size=7)
p
ggsave("Siggene_nMut_boxplot1.pdf",p,width = 6,height = 6)
ggsave("Siggene_freq_boxplot1.pdf",p,width = 6,height = 6)

