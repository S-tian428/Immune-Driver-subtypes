rm(list = ls())
library(MOVICS)
load("H:/免疫亚型分型/MOVICS/STAD.multiomics.rda")
surv.info<-STAD_multiomics$surv.info
load("H:/免疫亚型分型/MOVICS/cmoic.brca.rda")
setwd("H:/免疫亚型分型/result2/")
####获得突变的0 1矩阵
STAD.maf<-STAD_multiomics$STAD.maf
unique(STAD.maf$Tumor_Sample_Barcode)

#> --all samples matched.
##比较总突变负荷
library(maftools)
#maf<-GDCquery_Maf(tumor ="STAD", pipelines = "mutect2")
tmb.brca <- compTMB(moic.res     = cmoic.brca,
                    maf          = STAD.maf,
                    rmDup        = TRUE, # remove duplicated variants per sample
                    rmFLAGS      = FALSE, # keep FLAGS mutations
                    exome.size   = 38, # estimated exome size
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "DISTRIBUTION OF TMB AND TITV")


library(tidyr)
###长数据变宽数据(获得0,1矩阵)
mutation<-STAD.maf[,c(1,16)]
mutation$nMutation<-1
mutation<-aggregate(mutation$nMutation,by=list(mutation$Tumor_Sample_Barcode,mutation$Hugo_Symbol), FUN=sum)
mutation<-spread(mutation, key = "Group.1",
                         value = "x")
rownames(mutation)<-mutation[,1]
mutation<-mutation[,-1]
mutation[!is.na(mutation)]<-1
mutation[is.na(mutation)]<-0
#####
clust.res<-cmoic.brca$clust.res
clust.res<-clust.res[clust.res$samID%in%colnames(mutation),]
cmoic.brca$clust.res<-clust.res
### 比较突变频率
# mutational frequency comparison
mut.brca <- compMut(moic.res     = cmoic.brca,
                    mut.matrix   = mutation, # binary somatic mutation matrix
                    doWord       = TRUE, # generate table in .docx format
                    doPlot       = TRUE, # draw OncoPrint
                    freq.cutoff  = 0.05, # keep those genes that mutated in at least 5% of samples
                    p.adj.cutoff = 0.05, # keep those genes with adjusted p value < 0.05 to draw OncoPrint
                    innerclust   = TRUE, # perform clustering within each subtype
                    clust.col = c( "#FF9F1C","#E71D36"),
                   # mut.col = "#0077b6",
                    #bg.col = "white",
                    #annCol       = annCol, # same annotation for heatmap
                    #annColors    = annColors, # same annotation color for heatmap
                    width        = 6, 
                    height       = 4,
                    fig.name     = "ONCOPRINT FOR SIGNIFICANT MUTATIONS",
                    tab.name     = "INDEPENDENT TEST BETWEEN SUBTYPE AND MUTATION")
####为了量化这些可能影响免疫治疗的基因组改变，
####MOVICS 提供了两个函数来计算总突变负荷 （TMB） 
######和基因组改变分数 （FGA）。
###具体来说，TMB 是指在肿瘤基因组中发现的突变数量
###FGA 是受拷贝数增加或丢失影响的基因组百分比。
##这两个属性提供了有关肿瘤基因组构成的更深入信息。
# compare FGA, FGG, and FGL
#不仅计算 FGA，还计算每个子类型中
#每个样本的比增益 （FGG） 或损耗 （FGL）。
STAD.segment<-STAD_multiomics$STAD.segment
fga.brca <- compFGA(moic.res     = cmoic.brca,
                    segment      = STAD.segment,
                    iscopynumber = FALSE, # this is a segmented copy number file
                    cnathreshold = 0.2, # threshold to determine CNA gain or loss
                    test.method  = "nonparametric", # statistical testing method
                    fig.name     = "BARPLOT OF FGA")
