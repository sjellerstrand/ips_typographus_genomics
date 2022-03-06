### Calculate and plot outflank outliers

### Note: This analysis and script is a later addition upon request by the reviewers.
### Therefore, running this script might require more manual work when importing that data

rm(list=ls())

## Load libraries

library(tidyverse)
### Install required package devtools if not installed
if (!("devtools" %in% installed.packages())){install.packages("devtools")}
library(devtools)

if (!require("BiocManager", quietly = TRUE)){install.packages("BiocManager")}
library("BiocManager")

if (!("qvalue" %in% installed.packages())){BiocManager::install("qvalue")}
library("qvalue")

### Install OutFLANK if not installed
if (!("OutFLANK" %in% installed.packages())){install_github("whitlock/OutFLANK")}
library(OutFLANK) ### Used OutFLANK 0.2

if (!("vcfR" %in% installed.packages())){install.packages("vcfR")}
library(vcfR) ### Used version 1.12.0
library(ggplot2)
library(ggExtra)


## Upload Compressed VCF file and popmap
vcf <- read.vcfR("C:/Users/Simon/OneDrive - Lund University/Skrivbordet/outflank/bark_beetles_locations.vcf.gz", verbose=FALSE)
popmap <- read.delim("C:/Users/Simon/OneDrive - Lund University/Skrivbordet/outflank/popmap.txt", sep="\t", header=F)
pop <- matrix(NA, length(colnames(vcf@gt))-1, 2)
pop[,1] <- colnames(vcf@gt)[-1]

for(i in 1:nrow(popmap)) {
  pop[i,2] <- popmap[which(popmap[,1] == pop[i,1]),2]
}

ref <-  read.delim("C:/Users/Simon/OneDrive - Lund University/Skrivbordet/outflank/reference_customized.fasta.fai", sep="\t", head=F)
colnames(ref) <- c("contig", "length", "abs.pos","NA", "NA")

vcf2 <- extract.gt(vcf, element="GT")
vcf3 <- vcf2[which(is.biallelic(vcf)),]

for(i in 1:length(vcf3)) {
  if(!is.na(vcf3[i])) {
    vcf3[i] <- sum(as.numeric(unlist(strsplit(vcf3[i], split=c("[|///]")))))
  } else {
    vcf3[i] <- 9
  }
}

SNPdata <- t(vcf3)

## Calculate input data frame
FstDataFrame <- MakeDiploidFSTMat(SNPdata, colnames(SNPdata), pop[,2])

### Evaluate data
plot(FstDataFrame$FST, FstDataFrame$FSTNoCorr, xlim=c(-0.03,0.33), ylim=c(-0.03,0.33), pch=20)
abline(0,1)
plot(FstDataFrame$He, FstDataFrame$FSTNoCorr, pch=20, col="grey")
hist(FstDataFrame$FSTNoCorr, breaks=seq(0,0.33, by=0.001))
hist(FstDataFrame$FSTNoCorr[FstDataFrame$He>0.05], breaks=seq(0,0.33, by=0.001))
hist(FstDataFrame$FSTNoCorr[FstDataFrame$He>0.1], breaks=seq(0,0.3, by=0.001))




## Run outflank

outlier <- OutFLANK(FstDataFrame, NumberOfSamples=2)
OutFLANKResultsPlotter(outlier, withOutliers = TRUE,
                       NoCorr = TRUE, Hmin = 0.1, binwidth = 0.001, Zoom =
                         FALSE, RightZoomFraction = 0.05, titletext = NULL)

hist(out1$results$pvaluesRightTail)


sum(outlier$results$qvalues<0.01, na.rm=TRUE)
plot(outlier$results$He, outlier$results$FST, pch=20, col="grey")
points(outlier$results$He[outlier$results$qvalues<0.01], outlier$results$FST[outlier$results$qvalues<0.01], pch=21, col="blue")


top_candidates <- outlier$results$qvalues<0.01 & outlier$results$He>0.1

plot(outlier$results$He, outlier$results$FST, pch=20, col="grey")
points(outlier$results$He[top_candidates], outlier$results$FST[top_candidates], pch=21, col="blue")
topcan <- outlier$results[top_candidates,]
topcan[order(topcan$LocusName),]

results <- as.data.frame(matrix(NA, length(outlier$results$FST), 5))
colnames(results) <- c("CHROM", "pos", "val", "abs.pos", "outliers")

for(i in 1:nrow(results)) {
  results[i,1:3] <- c(unlist(strsplit(outlier$results[i,1], split="_")), outlier$results$FST[i])
  results[i,4] <-  ref$abs.pos[which(ref$contig == results[i,1])] + as.numeric(results[i,2])
}

results[top_candidates,5] <- "yes"


results[,3] <- as.numeric(results[,3])
results[,4] <- as.numeric(results[,4])
results[,5] <- as.factor(results[,5])
str(results)

scan <- ggplot(results, aes(abs.pos, val, color = outliers)) +
  geom_point(size = 7) + theme_light() +
  theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position = "none",
        axis.title.x=element_text(size=41), axis.text.x = element_blank(), axis.text.y = element_text(size=41, colour="black"),
        axis.title.y=element_text(size=41)) + ylab("OutFLANK Fst") + xlab("Reference genome") +
  scale_color_manual(values="red", na.value="#595959")
scan <- ggExtra::ggMarginal(scan, type = "histogram", margins="y", size=8)

jpeg("C:/Users/Simon/OneDrive - Lund University/Skrivbordet/outflank/Outlank_genome_scan.jpeg",
     width=3000, height=500, quality=100)
print(scan)
dev.off()
