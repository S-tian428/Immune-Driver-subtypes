rm(list = ls())
library(MOVICS)
#install.packages("matrixStats")
load("G:/免疫亚型分型/MOVICS/STAD_list.rda")
movics.data<-list(immgene.expr=STAD_list$immgene.expr,
                  driver.exp=STAD_list$driver.exp,
                  immgene.beta=STAD_list$immgene.beta,
                  driver.beta=STAD_list$driver.beta,
                  mut.status=STAD_list$mut.status)
surv.info<-STAD_list$surv.info
#####
setwd("H:/免疫亚型分型/MOVICS/")
optk.brca <- getClustNum(data        = movics.data, # 4种组学数据
                         is.binary   = c(F,F,F,F,T), #前3个不是二分类的，最后一个是
                         try.N.clust = 2:8, # 尝试亚型数量，从2到8
                         fig.path="H:/免疫亚型分型/MOVICS/",
                         fig.name    = "CLUSTER NUMBER OF TCGA-STAD")#保存的文件名
set.seed(123)
######根据九种方法分别做(加上贝叶斯方法结果不好)
moic.res.list <- getMOIC(data        = movics.data,
                         methodslist = list("SNF", "CIMLR", "PINSPlus", "NEMO", "COCA", "MoCluster",
                                            "LRAcluster", "ConsensusClustering", "IntNMF"), # 9种算法
                         N.clust     = 2,
                         type        = c("gaussian", "gaussian", "gaussian","gaussian","binomial"))

save(moic.res.list, file = "H:/免疫亚型分型/MOVICS/moic.res.list.rda")
load("H:/免疫亚型分型/MOVICS/moic.res.list.rda")
####聚类热图
cmoic.brca <- getConsensusMOIC(moic.res.list = moic.res.list,
                               fig.name      = "CONSENSUS HEATMAP",
                               distance      = "euclidean",
                               linkage       = "average")
save(cmoic.brca,file = "H:/免疫亚型分型/MOVICS/cmoic.brca.rda")
###分类效果
getSilhouette(sil=cmoic.brca$sil,
              fig.path = getwd(),
              fig.name = "SILHOUETTE",
              height = 5.5,
              width = 5)
####生存
surv.brca <- compSurv(moic.res = cmoic.brca,
                      surv.info = surv.info,
                      convt.time = "m",
                      surv.median.line = "h",
                      xyrs.est = c(5,10),
                      fig.name = "KAPLAN-MEIER CURVE OF CONSENSUSMOIC")

result<-cmoic.brca$clust.res
write.table(result,"H:/免疫亚型分型/MOVICS/subtype_result.txt",
            sep = "\t",row.names = F,quote = F)

STAD_FPKM<-read.table("H:/免疫亚型分型/Expression/STAD_FPKM.txt",
                      header = T,sep = "\t",stringsAsFactors = F,quote = "")
colnames(STAD_FPKM)<-gsub("\\.","-",colnames(STAD_FPKM))
View(head(STAD_FPKM))
index<-which(colnames(STAD_FPKM)%in%result$samID)
subtype_expression<-STAD_FPKM[,index]
write.table(subtype_expression,"H:/免疫亚型分型/MOVICS/subtype_expression.txt",
            sep = "\t",row.names = T,quote = F)


#####展示多组学数据
####亚型间特征比较
#####样本生存（按样本顺序排序）
#surv.info<- brca.tcga$clin.info
####根据聚类结果获取多组学热力图
indata <- movics.data
# data normalization for heatmap
plotdata <- getStdiz(data       = indata,
                     halfwidth  = c(2,2,2,2,NA), # no truncation for mutation
                     centerFlag = c(T,T,T,T,F), # no center for mutation
                     scaleFlag  = c(T,T,T,T,F)) # no scale for mutation
####得到分类后每个特征的基因排序
feat   <- moic.res.list$CIMLR$feat.res
####分别提取各组学的前10个基因展示名字
feat1  <- feat[which(feat$dataset == "dat1"),][1:10,"feature"] 
feat2  <- feat[which(feat$dataset == "dat2"),][1:10,"feature"]
feat3  <- feat[which(feat$dataset == "dat3"),][1:10,"feature"]
feat4  <- feat[which(feat$dataset == "dat4"),][1:10,"feature"]
feat5  <- feat[which(feat$dataset == "dat5"),][1:10,"feature"]
annRow <- list(feat1, feat2, feat3, feat4,feat5)
####生成颜色搭配
# library(paletteer)
paletteer_c("scico::berlin", n = 15)
paletteer_d("RColorBrewer::Paired")
paletteer_dynamic("cartography::green.pal", 5)
library(RColorBrewer)
# 为每个组学的热图自定义颜色，不定义也可
ImmgeneExp.col   <- c("#00FF00", "#008000", "white", "#800000", "#FF0000")
DriverExp.col <- rev(brewer.pal(11, "PiYG"))
ImmgeneMeth.col <- c("#0074FE", "#96EBF9", "#FEE900", "#F00003")
DriverMeth.col<-c("navy", "white", "firebrick3")
DriverMut.col    <- c("grey90" , "black")
col.list   <- list(ImmgeneExp.col, DriverExp.col,
                   ImmgeneMeth.col, DriverMeth.col,
                   DriverMut.col)
surv.info$stage[grep("stage iv",surv.info$stage)]<-"IV"
surv.info$stage[grep("stage iii.*",surv.info$stage)]<-"III"
surv.info$stage[grep("stage ii.*",surv.info$stage)]<-"II"
surv.info$stage[grep("stage i.*",surv.info$stage)]<-"I"

annCol    <- surv.info[,c("gender", "grade", "stage","response"), drop = FALSE]
surv.info$response[grep("",surv.info$response)]
index_response<-which(surv.info$response%in%"Complete Remission/Response"|surv.info$response%in%"Partial Remission/Response"|
                        surv.info$response%in%"Progressive Disease"|surv.info$response%in%"Stable Disease" )
index_null<-setdiff(1:nrow(surv.info),index_response)
surv.info$response[index_null]<-"not reported"
# generate corresponding colors for sample annotation
annColors <- list(gender  = c("female" = "#E31A1CFF",
                              "male"   = "#1F78B4FF"
),
stage = c("I"    = "#17692CFF",
          "II"    = "#B2DF8AFF",
          "III"    = "#FDBF6FFF",
          "IV"    = "#FF7F00FF", 
          "not reported" = "black"),
grade=c("G1"    = "#6A3D9AFF",
        "G2"    = "#CAB2D6FF",
        "G3"    = "#FFFF99FF",
        "GX"    = "#B15928FF"),
response=c("Complete Remission/Response"= "#1F78B4FF",
           "Partial Remission/Response"    = "#FDBF6FFF",
           "Progressive Disease"    = "#FF7F00FF",
           "Stable Disease"    = "#E31A1CFF",
           "not reported" = "black"))

###绘制两种亚型样本各亚型热图
getMoHeatmap(data          = plotdata,
             row.title     = c("Immgene.expr","driver.exp","immgene.beta","driver.beta","driver.mut.status"),
             is.binary     = c(F,F,F,F,T), # the 4th data is mutation which is binary
             legend.name   = c("Immgene.FPKM","Driver.FPKM","Immgene Beta Value","Driver Beta Value","Driver Mutated"),
             clust.res     = cmoic.brca$clust.res, # cluster results
             clust.dend    = cmoic.brca$clust.dend, # no dendrogram
             show.rownames = c(F,F,F,F,F), # specify for each omics data
             show.colnames = FALSE, # show no sample names
             annRow        = annRow, # mark selected features
             color         = col.list,
             annCol        = annCol, # no annotation for samples
             annColors     = annColors, # no annotation color
             width         = 10, # width of each subheatmap
             height        = 5, # height of each subheatmap
             fig.name      = "COMPREHENSIVE HEATMAP OF ICLUSTERBAYES")
