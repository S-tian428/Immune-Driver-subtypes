rm(list = ls())
#install.packages("oncoPredict")
library(oncoPredict)
#BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")
load("G:/免疫亚型分型/MOVICS/STAD.multiomics.rda")
STAD.FPKM<-STAD_multiomics$STAD.FPKM
####GDSC中包含免疫治疗药物
###GDSC2_Expr基于RNA-Seq 数据，覆盖更多细胞系和基因表达数据
trainingExprData=readRDS(file="G:/免疫亚型分型/result5/Oncopredict/DataFiles/Training Data/GDSC2_Expr (RMA Normalized and Log Transformed).rds")
trainingExprData[1:4,1:4]
trainingPtype=readRDS(file="H:/免疫亚型分型/result5/Oncopredict/DataFiles/Training Data/GDSC2_Res.rds")
trainingPtype<-exp(trainingPtype)
trainingPtype[1:4,1:4]
dim(trainingPtype)
dim(trainingExprData)
setwd("G:/免疫亚型分型/result5/Oncopredict")
calcPhenotype(trainingExprData=trainingExprData,
              trainingPtype=trainingPtype,
              testExprData=as.matrix(STAD.FPKM),
              batchCorrect="eb",
              powerTransformPhenotype=TRUE,
              removeLowVaryingGenes=0.2,
              minNumSamples=10,
              selection=1,
              printOutput=TRUE,
              pcr=FALSE,
              removeLowVaringGenesFrom="homogenizeData",
              report_pc=FALSE,
              cc=FALSE,
              percent=80,
              rsq=FALSE)
resultPtype <- read.csv('./calcPhenotype_Output/DrugPredictions.csv', 
                        header = T ,stringsAsFactors = F ,check.names = F)
colnames(resultPtype)<-gsub("_.*","",colnames(resultPtype))
rownames(resultPtype)<-resultPtype[,1]
resultPtype<-resultPtype[,-1]
surv.info<-STAD_multiomics$surv.info
index<-match(rownames(resultPtype),rownames(surv.info))
resultPtype$Cluster<-surv.info$Cluster
i=1
wilcox_test_result<-NULL
for (i in 1:(ncol(resultPtype)-1)) {
  wt<- wilcox.test(resultPtype[,i] ~ Cluster, data = resultPtype)
  wt_p<-data.frame(Drug=colnames(resultPtype)[i],P=wt$p.value)
  wilcox_test_result<-rbind(wilcox_test_result,wt_p)
}

wilcox_test_result<-wilcox_test_result[wilcox_test_result$P<0.05,]
a<-aggregate(resultPtype[,1:198],by=list(resultPtype$Cluster),mean)
rownames(a)<-a[,1]
a<-a[,-1]
IC50_mean<-t(a)
IC50_mean<-data.frame(IC50_mean,stringsAsFactors = F)
index1<-match(wilcox_test_result$Drug,rownames(IC50_mean))
result<-cbind(wilcox_test_result,IC50_mean[index1,])
str(result)
result$Difference=result$CS2-result$CS1
result$response_cluster<-"CS1"
index_CS2<-which(result$Difference<0)
result$response_cluster[index_CS2]<-"CS2"
write.table(result,"Wilcox_test.txt",sep="\t",quote=FALSE,row.names = F)
#####绘制热图
result<-result[order(result$response_cluster),]
resultPtype<-resultPtype[order(resultPtype$Cluster),]
heatmap_onco<-resultPtype[,result$Drug]
#heatmap_onco<-t(heatmap_onco)
####注释行列
annotation_row<-data.frame(clust=resultPtype$Cluster)
rownames(annotation_row) <- rownames(resultPtype)
annotation_col<-data.frame(drug=paste(result$response_cluster,"Response",sep = " "))
rownames(annotation_col)<-result$Drug

colcolor <- c("#0072b5","#fdb462") 
names(colcolor) <- c("CS1 Response","CS2 Response") #类型颜色

rowcolor <- c("#2ec4b6","#e71d36") 
names(rowcolor) <- c("CS1","CS2") #类型颜色
ann_colors<-list(drug=colcolor,clust=rowcolor)
#mycol<-colorRampPalette(c("#17692c", "white", "#ff7f00"))(100)
#mycol<-colorRampPalette(c("#6a3d9a", "white", "#b15928"))(100)
library(pheatmap)
bk <- c(seq(-2,-0.1,by=0.01),seq(0,2,by=0.01))#e71d36 #0f7a3f#507aaf
p<-pheatmap(as.matrix(heatmap_onco)
            ,color = c(colorRampPalette(colors = c("#397FC7","white"))(length(bk)/2),
                       colorRampPalette(colors = c("white","#F1B656"))(length(bk)/2)),
            # legend_breaks=seq(-8,8,2),
            breaks=bk,
            cluster_rows = F,cluster_cols = F,fontsize_col =6,fontsize_row=6,
            scale = 'column',show_colnames = T,show_rownames = F,
            annotation_col = annotation_col,
            annotation_row =annotation_row ,
            annotation_colors = ann_colors)
library(ggpubr)
ggsave("Drug_heatmap.pdf",p,width = 8,height = 4)
