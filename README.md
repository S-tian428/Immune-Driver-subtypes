Pipeline of this study
1.	Multi-omic data for gastric adenocarcinoma (STAD) of the TCGA project, including transcriptomic expression, DNA methylation, somatic mutations, copy number variations, as well as clinical information, were downloaded from the UCSC Xena database. 
2.	1,794 immune-related genes were obtained from the ImmPort database. Immune-related lncRNAs and miRNAs were downloaded from the ImmLnc database. 110 driver genes associated with gastric cancer were retrieved from the NCG database.
3.	We identified 123 differentially expressed immune-related lncRNAs, 137 differentially expressed immune-related miRNAs, and 196 DEGs using log2(FPKM+1) transcriptomic expression by limma.
4.	We established Immune&Driver molecular subtypes CS1 and CS2 using the MOVICS package by integrating multi-omics data of immune-related genes and driver genes.
5.	We performed univariable and multivariable regression Cox analyses to explore the influence of treatment exposure and clinicopathological factors on prognostic differences between CS1 and CS2.
6.	We compared the genomic heterogeneity and genome instability between CS1 and CS2.
7.	We conducted a comprehensive comparative analysis of the tumor immune microenvironment (TiME) between CS1 and CS2.
8.	We compared the potential response to immune checkpoint and chemokine therapy between CS1 and CS2.
9.	We incorporated three additional gastric cancer cohorts and successfully reproduced the CS1 and CS2 subtypes using subtype-specific markers by the Nearest Template Prediction (NTP) algorithm.
10.	We predicted and compared the potential drug sensitivity between CS1 and CS2 subtypes across the TCGA cohort and three additional gastric cancer cohorts.
11.	To explore the potential relevance of the Immune&Driver subtypes to immunotherapy response, we reproduced CS1- and CS2-like transcriptomic patterns in multiple melanoma datasets using subtype-specific markers by NTP.
12.	We identified ten CS biomarkers associated with both immunotherapy response and overall survival. Based on the transcriptional expression of these biomarkers, we developed a predictive model for immunotherapy response using seven machine-learning algorithms.
