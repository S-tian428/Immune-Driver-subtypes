rm(list = ls())
library(data.table)
####免疫相关的lncRNA
###免疫通路相关的lncRNA
lncRNA_Pathway<-read.table("H:/免疫亚型分型/result1/ImmlncRNA/STAD/Immune_pathway_sig.txt",
                           header = T,sep = "\t",stringsAsFactors = F,quote = "")
PathwaylncRNA<-unique(lncRNA_Pathway$lncRNA.Symbol)
###免疫细胞相关的lncRNA
lncRNA_cells<-read.table("H:/免疫亚型分型/result1/ImmlncRNA/STAD/Immune_cell_sig.txt",
                         header = T,sep = "\t",stringsAsFactors = F,quote = "")
#lncRNA_cells<-lncRNA_cells[abs(lncRNA_cells$R.Value)>0.5,]
CelllncRNA<-unique(lncRNA_cells$lncRNA.Symbol)
Imm_lncRNA<-union(PathwaylncRNA,CelllncRNA)
######差异表达lncRNA
lncRNA_Exp<-read.table("H:/免疫亚型分型/result1/ImmlncRNA/STAD/Diff_expression_sig.txt",
                       header = T,sep = "\t",stringsAsFactors = F,quote = "",fill = T)
EXPlncRNA<-unique(lncRNA_Exp$lncRNA.Symbol)
###差异表达的免疫相关的lncRNA 608个
Imm_lncRNA<-Reduce(intersect, list(Imm_lncRNA, EXPlncRNA))
#####展示其中哪些是通路相关的lncRNA,哪些是细胞相关的lncRNA
Imm_CelllncRNA<-intersect(CelllncRNA,Imm_lncRNA)
Imm_PathwaylncRNA<-intersect(PathwaylncRNA,Imm_lncRNA)
####绘制韦恩图
###差异表达免疫通路相关的lncRNA
####差异表达免疫细胞相关的lncRNA
lncRNA<-list("ImmCell-lncRNA"=Imm_CelllncRNA,
             "ImmPathway-lncRNA"=Imm_PathwaylncRNA)
library(VennDiagram)
setwd("H:/免疫亚型分型/result1/")
venn.plot <-venn.diagram(
  lncRNA,
  category.names = c("Immune cells-related lncRNA",
                     "Immune pathways-related lncRNA"),
  filename = NULL,
  cex=2,
  #lty = "dotted",
  compression = "lzw",
  lwd = 2,
  col = c("#440154ff",  "#21908dff"),
  fill = c("#440154ff", "#21908dff"),
  alpha = 0.30,
  cat.cex=2,
  )
dev.off()
pdf('Imm-pathway-cell-lncRNA_Venn.pdf', width = 5, height = 5)
grid.draw(venn.plot)
dev.off()
#####设计哪些免疫通路和免疫细胞
Immpathway<-lncRNA_Pathway[lncRNA_Pathway$lncRNA.Symbol%in%Imm_PathwaylncRNA,3:4]
Immpathway<-data.frame(table(Immpathway$Immune.Pathway))
colnames(Immpathway)<-c("Pathway","lncRNA Number")
library(ggplot2)
p <- ggplot(Immpathway,aes(x=`lncRNA Number`,
                     y=factor(Pathway)))+
  geom_point(color="#9acd32",size=5)+
  geom_segment(aes(x=0,xend=`lncRNA Number`,
                   y=Pathway,yend=Pathway
                  ),color="#9acd32",size=1
              )+
  labs(y=NULL)+
  theme_bw()
p
ggsave("Immune Pathway-lncRNA.pdf",p,width=5,height = 5)
####免疫细胞
Immcell<-lncRNA_cells[lncRNA_cells$lncRNA.Symbol%in%Imm_CelllncRNA,]
Immcell<-Immcell[order(Immcell$Immune.Cell),]
####数据太大了，做成附表
write.table(Immcell,"TableS1_Immcell-lncRNA.txt",
            sep="\t",quote=FALSE,row.names = F)


####免疫相关的miRNA
###免疫通路相关的miRNA
miRNA_Pathway<-read.table("H:/免疫亚型分型/result1/ImmMiRNA/STAD/Immune_pathway_sig.txt",
                          header = T,sep = "\t",stringsAsFactors = F,quote = "")
PathwaymiRNA<-unique(miRNA_Pathway$miRNA.Symbol)
###免疫细胞相关的miRNA
miRNA_cells<-read.table("H:/免疫亚型分型/result1/ImmMiRNA/STAD/Immune_cell_sig.txt",
                        header = T,sep = "\t",stringsAsFactors = F,quote = "")
CellmiRNA<-unique(miRNA_cells$miRNA.Symbol)
Imm_miRNA<-union(PathwaymiRNA,CellmiRNA)
######差异表达miRNA
miRNA_Exp<-read.table("H:/免疫亚型分型/result1/ImmMiRNA/STAD/Diff_expression_sig.txt",
                      header = T,sep = "\t",stringsAsFactors = F,quote = "",fill = T)
EXPmiRNA<-unique(miRNA_Exp$miRNA.Symbol)
#Imm_miRNA_union<- Reduce(union, list(PathwaymiRNA, CellmiRNA,EXPmiRNA))
#write.table(Imm_miRNA_union,"H:/免疫亚型分型/ImmMiRNA/STAD/mirBase_ID_union.txt",
#            sep = "\t",row.names = F,quote = F)
####取三组miRNA的交集：193个
Imm_miRNA<- Reduce(intersect, list(Imm_miRNA,EXPmiRNA))
#####展示其中哪些是通路相关的miRNA,哪些是细胞相关的miRNA
Imm_CellmiRNA<-intersect(CellmiRNA,Imm_miRNA)
Imm_PathwaymiRNA<-intersect(PathwaymiRNA,Imm_miRNA)
####绘制韦恩图
###差异表达免疫通路相关的lncRNA
####差异表达免疫细胞相关的lncRNA
miRNA<-list("ImmCell-miRNA"=Imm_CellmiRNA,
             "ImmPathway-miRNA"=Imm_PathwaymiRNA)
library(VennDiagram)
setwd("H:/免疫亚型分型/result1/")
venn.plot <-venn.diagram(
  miRNA,
  category.names = c("Immune cells-related miRNA",
                     "Immune pathways-related miRNA"),
  filename = NULL,
  cex=2,
  #lty = "dotted",
  compression = "lzw",
  lwd = 2,
  col = c("#440154ff",  "#21908dff"),
  fill = c("#440154ff", "#21908dff"),
  alpha = 0.30,
  cat.cex=2,
)
dev.off()
pdf('Imm-pathway-cell-miRNA_Venn.pdf', width = 5, height = 5)
grid.draw(venn.plot)
dev.off()
#####设计哪些免疫通路和免疫细胞
Immpathway<-miRNA_Pathway[miRNA_Pathway$miRNA.Symbol%in%Imm_PathwaymiRNA,3:4]
Immpathway<-data.frame(table(Immpathway$Immune.Pathway))
colnames(Immpathway)<-c("Pathway","miRNA Number")
library(ggplot2)
p <- ggplot(Immpathway,aes(x=`miRNA Number`,
                           y=factor(Pathway)))+
  geom_point(color="#ffa500",size=5)+
  geom_segment(aes(x=0,xend=`miRNA Number`,
                   y=Pathway,yend=Pathway
  ),color="#ffa500",size=1
  )+
  labs(y=NULL)+
  theme_bw()
p
ggsave("Immune Pathway-miRNA.pdf",p,width=5,height = 5)

####免疫细胞
Immcell<-miRNA_cells[miRNA_cells$miRNA.Symbol%in%Imm_CellmiRNA,]
unique(Immcell$Immune.Cell)
Immcell<-Immcell[order(Immcell$Immune.Cell),]
####数据太大了，做成附表
write.table(Immcell,"TableS2_Immcell-miRNA.txt",
            sep="\t",quote=FALSE,row.names = F)

####这些miRNA对应的是miRBase ID,和表达谱中的miRNA名字对不上
###需要利用ensemble BioMart工具转化
####获得miRNA对应的mirBase ID
miRNA_mirBase<-read.table("H:/免疫亚型分型/ImmMiRNA/STAD/miRNA_mirBase.txt",
                          header = F,sep = "\t",stringsAsFactors = F,quote = "")
miRNA_mirBase$V7<-gsub("Name=","",miRNA_mirBase$V7)
index<-which(miRNA_mirBase$V7%in%Imm_miRNA)
miRNA_mirBase<-miRNA_mirBase[index,]
miRBase_id<-unique(miRNA_mirBase$V8)
###得到过度ID，可以输入到ensemble中
miRBase_id<-gsub("Derives_from=","",miRBase_id)
write.table(miRBase_id,"H:/免疫亚型分型/ImmMiRNA/STAD/ensemble_input.txt",
            sep = "\t",row.names = F,quote = F)
###读入ensemble准化后的miRNA名
ensemble_mirRNA<-read.table("H:/免疫亚型分型/result1/ImmMiRNA/STAD/ensemble_out1.txt",
                            header = T,sep = "\t",stringsAsFactors = F,quote = "")
###184个
Imm_miRNA<-ensemble_mirRNA$Gene.name

####免疫相关的gene
Imm_gene<-read.table("G:/免疫亚型分型/result1/ImmGene/Immport.txt",
                     header = T,sep = "\t",stringsAsFactors = F,quote = "")
Immgene<-unique(Imm_gene$Symbol)
####差异免疫基因：196个
####自己识别差异免疫基因
limma_result<-read.table("G:/免疫亚型分型/Expression/FPKM_limma.txt",
                         header = T,sep = "\t",stringsAsFactors = F,quote = "")
###筛选条件
limma_result<-limma_result[abs(limma_result$logFC)>1 & limma_result$adj.P.Val<0.05,]
####196
Immgene<-intersect(Immgene,rownames(limma_result))
Imm_limma<-limma_result[Immgene,]
write.table(Imm_limma,"G:/免疫亚型分型/result1/ImmGene/Immgene_siglimma.txt",
            sep = "\t",row.names = T,quote = F)
#
Imm_gene_result<-Imm_gene[Imm_gene$Symbol%in%Immgene,]
Immpathway<-data.frame(table(Imm_gene_result$Category))
colnames(Immpathway)<-c("Pathway","mRNA Number")
library(ggplot2)
p <- ggplot(Immpathway,aes(x=`mRNA Number`,
                           y=factor(Pathway)))+
  geom_point(color="#8b5a2b",size=5)+
  geom_segment(aes(x=0,xend=`mRNA Number`,
                   y=Pathway,yend=Pathway
  ),color="#8b5a2b",size=1
  )+
  labs(y=NULL)+
  theme_bw()
p
ggsave("Immune Pathway-mRNA.pdf",p,width=5,height = 5)
###表达热图
STAD_FPKM<-read.table("H:/免疫亚型分型/Expression/STAD_FPKM.txt",
                      header = T,sep = "\t",stringsAsFactors = F,quote = "")
colnames(STAD_FPKM)<-gsub("\\.","-",colnames(STAD_FPKM))
limma_result_Order<-limma_result[Immgene,]
limma_result_Order<-limma_result_Order[order(limma_result_Order$logFC,decreasing = T),]
Imm_FPKM<-STAD_FPKM[rownames(limma_result_Order),]
z1<-grep("-11.$",colnames(Imm_FPKM),value = T)
z2<-setdiff(colnames(Imm_FPKM),z1)
Imm_FPKM<-Imm_FPKM[,c(z1,z2)]
Imm_FPKM<-as.matrix(Imm_FPKM)
bk <- c(seq(-1,-0.1,by=0.01),seq(0,1,by=0.01))
library(ComplexHeatmap)
library(RColorBrewer)
# 做热图：

annotation_col = data.frame(c(rep("Normal",length(z1)),
                              rep("Tomor",length(z2))))
row.names(annotation_col) <- colnames(Imm_FPKM)
colnames(annotation_col)<-"SampleType"
annotation_colors =list(SampleType=c("Normal"="#27408a",
                      "Tomor"="#228b22")) 
#c(colorRampPalette(colors = c("darkgreen","gray99"))(length(bk)/2),colorRampPalette(colors = c("grey99","red"))(length(bk)/2))
p<-pheatmap(Imm_FPKM,
            annotation_col = annotation_col,
            annotation_colors=annotation_colors,
            border_color = NA,na_col = "grey",
            cluster_col = F,cluster_rows =F,
            show_rownames = F,
            show_colnames = F,
            scale="row",
            #color = colorRampPalette(brewer.pal(9,"Blues")[6:1],
            color = c("blue","white","red"),
            # color=mycol,
            #display_numbers = P_wide,
            fontsize_number = 12, number_color = "black",
            # legend_breaks=seq(-8,8,2),
           #breaks=bk
           )
pdf("mRNA_heatmap.pdf",width =6,height =5)
p
dev.off()



#####最终使用的免疫相关基因集
Imm_signiture_result<-c(Imm_lncRNA,Imm_miRNA,Immgene)

#####STAD driver gene
driver_gene<-read.table("H:/免疫亚型分型/NCG/NCG_driver_STAD.txt",
                        header = T,sep = "\t",stringsAsFactors = F,quote = "")
####110个
driver_gene<-unique(driver_gene$symbol)


####MOVICS
#library(MOVICS)
####基因表达谱
STAD_FPKM<-read.table("H:/免疫亚型分型/Expression/STAD_FPKM.txt",
                      header = T,sep = "\t",stringsAsFactors = F,quote = "")
colnames(STAD_FPKM)<-gsub("\\.","-",colnames(STAD_FPKM))
####提取免疫相关lncRNA,miRNA,mRNA表达谱
signiture.exp<-STAD_FPKM[Imm_signiture_result,]
signiture.exp<-na.omit(signiture.exp)

driver.exp<-STAD_FPKM[driver_gene,]
driver.exp<-na.omit(driver.exp)

######甲基化数据
STAD_methylation<-fread("H:/免疫亚型分型/原始多组学数据/TCGA-STAD.methylation450.tsv.gz",header = T, sep = '\t',data.table = F)
#####将甲基化位点对应到基因上
gene_id<-fread("H:/免疫亚型分型/原始多组学数据/illuminaMethyl450_hg38_GDC",header = T, sep = '\t',data.table = F,fill=T)
#####按照基因将一行拆分成多行
library(tidyr)
gene_id <-gene_id %>% as_tibble() %>% 
  separate_rows(gene, sep = ",")

Immgene_id<-gene_id[gene_id$gene%in%Imm_signiture_result,1:2]
driver_id<-gene_id[gene_id$gene%in%driver_gene,1:2]
###319个基因匹配到甲基化探针
###Imm gene 492
#unique(Immgene_id$gene)
###driver 107
#unique(driver_id$gene)
###提取免疫基因甲基化谱
Immgene_methylation=merge(Immgene_id,STAD_methylation,by.y ="Composite Element REF",by.x = "#id" )
#View(head(STAD_methylation))
###基因的平均甲基化水平
Immgene_methylation<-aggregate(Immgene_methylation[,3:ncol(Immgene_methylation)],by=list(Immgene_methylation$gene),mean)
rownames(Immgene_methylation)<-Immgene_methylation$Group.1
Immgene_methylation<-Immgene_methylation[,-1]

#####提取driver gene 甲基化
Drivergene_methylation=merge(driver_id,STAD_methylation,by.y ="Composite Element REF",by.x = "#id" )
###基因的平均甲基化水平
Drivergene_methylation<-aggregate(Drivergene_methylation[,3:ncol(Drivergene_methylation)],by=list(Drivergene_methylation$gene),mean)
rownames(Drivergene_methylation)<-Drivergene_methylation$Group.1
Drivergene_methylation<-Drivergene_methylation[,-1]
####基因突变谱
STAD_mutations<-fread("H:/免疫亚型分型/原始多组学数据/TCGA-STAD.mutect2_snv.tsv.gz",header = T, sep = '\t',data.table = F)

Immgene_mutation<-STAD_mutations[STAD_mutations$gene%in%Imm_signiture_result,]
###长数据变宽数据(获得0,1矩阵)
Immgene_mutation<-Immgene_mutation[,1:2]
Immgene_mutation$nMutation<-1
Immgene_mutation<-aggregate(Immgene_mutation$nMutation,by=list(Immgene_mutation$Sample_ID,Immgene_mutation$gene), FUN=sum)
Immgene_mutation<-spread(Immgene_mutation, key = "Group.1",
                         value = "x")
rownames(Immgene_mutation)<-Immgene_mutation[,1]
Immgene_mutation<-Immgene_mutation[,-1]
Immgene_mutation[!is.na(Immgene_mutation)]<-1
Immgene_mutation[is.na(Immgene_mutation)]<-0

driver_mutation<-STAD_mutations[STAD_mutations$gene%in%driver_gene,]
###长数据变宽数据(获得0,1矩阵)
driver_mutation<-driver_mutation[,1:2]
driver_mutation$nMutation<-1
driver_mutation<-aggregate(driver_mutation$nMutation,by=list(driver_mutation$Sample_ID,driver_mutation$gene), FUN=sum)
driver_mutation<-spread(driver_mutation, key = "Group.1",
                        value = "x")
rownames(driver_mutation)<-driver_mutation[,1]
driver_mutation<-driver_mutation[,-1]
driver_mutation[!is.na(driver_mutation)]<-1
driver_mutation[is.na(driver_mutation)]<-0

###突变和表达样本数量和顺序需要一致
intersect_sample <-Reduce(intersect, 
                          list(colnames(signiture.exp), colnames(driver.exp),
                               colnames(Immgene_methylation),colnames(Drivergene_methylation),
                               colnames(driver_mutation)))
#grep("-11",sample_intersect)
###和突变数据样本取交集后，都是癌症组织数据，没有癌旁
###不需要再筛选
###最终的输入数据
###表达
signiture.exp<-signiture.exp[,intersect_sample]
#删除值全为0的行
signiture.exp <- signiture.exp[rowSums(signiture.exp != 0) > 0, ]
#删除值全为0的列
signiture.exp <- signiture.exp[,colSums(signiture.exp != 0) > 0 ]

driver.exp<-driver.exp[,intersect_sample]
#删除值全为0的行
driver.exp <- driver.exp[rowSums(driver.exp != 0) > 0, ]
#删除值全为0的列
driver.exp <- driver.exp[,colSums(driver.exp != 0) > 0 ]

###甲基化
Immgene_methylation<-Immgene_methylation[,intersect_sample]
Immgene_methylation[is.na(Immgene_methylation)]<-0
#删除值全为0的行
Immgene_methylation <- Immgene_methylation[rowSums(Immgene_methylation != 0) > 0, ]
#删除值全为0的列
Immgene_methylation <- Immgene_methylation[,colSums(Immgene_methylation != 0) > 0 ]


Drivergene_methylation<-Drivergene_methylation[,intersect_sample]
Drivergene_methylation[is.na(Drivergene_methylation)]<-0
#删除值全为0的行
Drivergene_methylation <- Drivergene_methylation[rowSums(Drivergene_methylation != 0) > 0, ]
#删除值全为0的列
Drivergene_methylation <- Drivergene_methylation[,colSums(Drivergene_methylation != 0) > 0 ]

####突变
driver_mutation<-driver_mutation[,intersect_sample]
#删除值全为0的行
driver_mutation <- driver_mutation[rowSums(driver_mutation != 0) > 0, ]
#删除值全为0的列
driver_mutation <- driver_mutation[,colSums(driver_mutation != 0) > 0 ]

STAD_survival<-read.table("H:/免疫亚型分型/原始多组学数据/STAD_survival.txt",
                          header = T,sep = "\t",stringsAsFactors = F,quote = "")
sample_intersect1<-gsub("-01A","-01",intersect_sample)
sample_intersect1<-gsub("-01B","-01",sample_intersect1)
index<-match(sample_intersect1,STAD_survival$sample)
STAD_survival<-STAD_survival[index,]
#####样本临床信息
STAD_clinical<-fread("H:/免疫亚型分型/原始多组学数据/TCGA-STAD.GDC_phenotype.tsv.gz",header = T, sep = '\t',data.table = F)
STAD_clinical<-data.frame(sample=STAD_clinical$submitter_id.samples,
                          age=STAD_clinical$age_at_diagnosis.diagnoses,
                          grade=STAD_clinical$neoplasm_histologic_grade,
                          gender=STAD_clinical$gender.demographic,
                          stage=STAD_clinical$tumor_stage.diagnoses,
                          stageM=STAD_clinical$pathologic_M,
                          stageN=STAD_clinical$pathologic_N,
                          stageT=STAD_clinical$pathologic_T,
                          response=STAD_clinical$followup_treatment_success)
index<-match(intersect_sample,STAD_clinical$sample)
STAD_clinical<-STAD_clinical[index,]
surv.info<-data.frame(fustat=STAD_survival$OS,futime=STAD_survival$OS.time,
                      STAD_clinical[,2:ncol(STAD_clinical)])
rownames(surv.info)<-STAD_clinical$sample


STAD_list <- list(immgene.expr=signiture.exp,
                  driver.exp=driver.exp,
                  immgene.beta=Immgene_methylation,
                  driver.beta=Drivergene_methylation,
                  mut.status=driver_mutation,
                  surv.info=surv.info)   
save(STAD_list, file = "H:/免疫亚型分型/MOVICS/STAD_list.rda")


