rm(list = ls())  
library("GenomicFeatures")
library("ChIPseeker")
library("txdbmaker")

setwd("D:/Home/OsWRKY-NGS/ATAC-seq")

spombe <- makeTxDbFromGFF("./Oryza_sativa.IRGSP-1.0.57.gff3")

Mock <- readPeakFile('./peaks/Mock.bed')
RSV <- readPeakFile('./peaks/RSV.bed')

Mock_peakAnno <- annotatePeak(Mock, tssRegion =c(-3000, 3000), TxDb = spombe)
RSV_peakAnno <- annotatePeak(RSV, tssRegion =c(-3000, 3000), TxDb = spombe)

# plotAnnoPie(Mock_peakAnno)
plotAnnoPie(RSV_peakAnno)
