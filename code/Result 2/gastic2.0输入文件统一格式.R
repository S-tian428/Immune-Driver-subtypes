rm(list = ls())
Segment<-read.table("G:/免疫亚型分型/result2/Gistic2.0/STAD_MaskedCopyNumberSegment.txt",
                                          header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
####只保留-01A
Segment$Sample <- substring(Segment$Sample,1,16)
Segment<-Segment[grep("-01A$",Segment$Sample),]
#Segment$Chromosome<-gsub("X","23",Segment$Chromosome)
#####亚型
load("G:/免疫亚型分型/MOVICS/cmoic.brca.rda")
Sample_cluster<-cmoic.brca$clust.res
CS1_sample<-Sample_cluster$samID[Sample_cluster$clust%in%"1"]
CS2_sample<-Sample_cluster$samID[Sample_cluster$clust%in%"2"]
CS1_Segment<-Segment[Segment$Sample%in%CS1_sample,]
CS2_Segment<-Segment[Segment$Sample%in%CS2_sample,]
write.table(CS1_Segment,"G:/免疫亚型分型/result2/Gistic2.0/Input/CS1_Segment.txt",
            row.names=F,col.names=F,sep="\t",quote=FALSE)
write.table(CS2_Segment,"G:/免疫亚型分型/result2/Gistic2.0/Input/CS2_Segment.txt",
            row.names=F,col.names=F,sep="\t",quote=FALSE)

SNP6<-read.table("G:/免疫亚型分型/result2/Gistic2.0/snp6.na35.remap.hg38.subset.txt",
                    header=T,sep = "\t", quote = "\"'",stringsAsFactors=FALSE)
SNP6<-SNP6[SNP6$freqcnv%in%"FALSE",]
SNP6<-SNP6[,1:3]
SNP6<-SNP6[order(SNP6$chr),]
#colnames(SNP6)<-c("Marker Name","Chromosome","Marker Position")
write.table(SNP6,"G:/免疫亚型分型/result2/Gistic2.0/Input/Marker_file.txt",
            row.names=F,col.names=F,sep="\t",quote=FALSE)
