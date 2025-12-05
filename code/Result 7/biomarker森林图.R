#install.packages("forestplot")
library(forestplot)
rm(list=ls())
OS_result<-read.table("G:/免疫亚型分型/result8/Biomarker_ICB/Sig_Cox_result.txt",
                      header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
OS_result$pvalue=round(OS_result$pvalue,3)
OS_result$HR=round(OS_result$HR,3)
OS_result$LowerCI=round(OS_result$LowerCI,3)
OS_result$upperCI=round(OS_result$upperCI,3)
OS_result$gene<-rownames(OS_result)
OS_result$HR..95..CI.<-paste0(OS_result$HR,"( ",OS_result$LowerCI," to ",OS_result$upperCI," )")
###森林图左侧
pdf("G:/免疫亚型分型/result8/Biomarker_ICB/forestplot.pdf",width = 7,height = 5)
n=nrow(OS_result)
nRow=n+1
ylim=c(1,nRow)
layout(matrix(c(1,2),nc=2),width=c(3,2))

#森林图左边的基因信息
xlim = c(0,3)
par(mar=c(4,2,1.5,1.5))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,xlab="",ylab="")#axes=F表示禁止生成坐标轴
text.cex=0.8 #放大0.8倍
text(0,n:1,OS_result$gene,adj=0,cex=text.cex);text(0,n+1,'gene',cex=text.cex,font=2,adj=0)#显示基因列信息
text(1.5-0.5*0.2,n:1,OS_result$pvalue,adj=1,cex=text.cex);text(1.5-0.5*0.2,n+1,'pvalue',cex=text.cex,font=2,adj=1)#显示pvalue列信息
text(3,n:1,OS_result$HR..95..CI.,adj=1,cex=text.cex);text(3,n+1,'Hazard ratio',cex=text.cex,font=2,adj=1,)#显示计算出的Hazard ratio列信息
####右侧森林图
par(mar=c(4,1,1.5,1),mgp=c(2,0.5,0))
xlim = c(0,max(as.numeric(OS_result$LowerCI),as.numeric(OS_result$upperCI)))
plot(1,xlim=xlim,ylim=ylim,type="n",axes=F,ylab="",xaxs="i",xlab="Hazard ratio")#设置x轴的标题
arrows(as.numeric(OS_result$LowerCI),n:1,as.numeric(OS_result$upperCI),n:1,angle=90,code=3,length=0.03,col="black",lwd=2)
abline(v=1,col="black",lty=2,lwd=2) #添加中线，设置中线的位置，颜色，类型，宽度
boxcolor = ifelse(as.numeric(OS_result$HR) > 1, '#B10026', '#1F77B4')#设置中线的取值
points(as.numeric(OS_result$HR), n:1, pch = 15, col = boxcolor, cex=2)#pch表示点的样式，设置点的大小，颜色
axis(1)
dev.off()
