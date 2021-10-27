#!/usr/bin/Rscript

### Perform multivariate analyses with principal components
## Export variables and load libraries
library(plyr)
library(tidyverse)
library(viridis)
library(heplots)
library(MASS)

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
  assign(args[[i]][1], args[[i]][2])
}

## Read in data
eigvec <- read.delim(paste(WORKDIR, "/", DATA, "/", VCF, "_", DATA, ".eigenvec", sep=""), sep=" ", head=F)
eigvec <- eigvec[,-1]

## Set names
names(eigvec)[1] <- "Individual"
names(eigvec)[2:ncol(eigvec)] <- paste0("PC", 1:(ncol(eigvec)-1))
PCs <- as.numeric(PCs)

## Sort out the individual sample locations and geographic regions
info <- read.delim(paste(METADATA, "/sample_info_", DATA, ".txt", sep=""), header=T)
order <- read.delim(paste(METADATA, "/order_", DATA, ".txt", sep=""), header=F)
translate <- read.delim(paste(METADATA, "/translate_population_names.txt", sep=""), header=F)
eigvec <-  cbind(eigvec[,1:(1+PCs)], rep(NA, nrow(eigvec)), rep(NA, nrow(eigvec)))
colnames(eigvec)[(ncol(eigvec)-1):ncol(eigvec)] <- c("Locality", "Region")
for(i in 1:nrow(eigvec)) {
eigvec[i,(ncol(eigvec)-1):ncol(eigvec)] <- as.matrix(info[which(eigvec[i,1] == info[,1]),3:4])
if(eigvec[i,ncol(eigvec)] != "Norrbotten") {
eigvec[i,ncol(eigvec)] <- "South"
}
}
eigvec$Locality <- factor(eigvec$Locality, levels=rev(order[,1]))
eigvec$Locality <- mapvalues(eigvec$Locality, from = translate[,1], to =  translate[,2])
eigvec$Region <- factor(eigvec$Region, levels=unique(eigvec$Region))

## Model variables
prin_comp <- as.matrix(eigvec[,2:(PCs+1)])
Locality <- eigvec$Locality
Region <- eigvec$Region

## MANOVA
#sampling locations
model1 <- lm(prin_comp ~ Locality)
summary.manova1 <- summary(manova(model1))
summary.aov.manova1 <- summary.aov(manova(model1))
effect.size1 <- etasq(model1, anova=T, partial=T)

# north south
model2 <- lm(prin_comp ~ Region)
summary.manova2 <- summary(manova(model2))
summary.aov.manova2 <- summary.aov(manova(model2))
effect.size2 <- etasq(model2, anova=T, partial=T)

# Write MANOVA statistics
sink(paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, "_summary_manova_sampling_locations.txt", sep=""))
print(summary.manova1)
print(summary.aov.manova1)
print(effect.size1)
sink()
sink(paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, "_summary_manova_regions.txt", sep=""))
print(summary.manova2)
print(summary.aov.manova2)
print(effect.size2)
sink()

## LDA
# sampling locations
model3 <- lda(Locality ~ prin_comp, CV=T)
ct3 <-table(Locality, model3$class)
correct_assign3 <- as.data.frame(cbind(diag(prop.table(ct3, 1)), colnames(prop.table(ct3, 1))))
accuracy3 <- c("model accuracy", sum(diag(prop.table(ct3))))
lda_summary3 <- cbind(eigvec, model3$posterior)

# north south
model4 <- lda(Region ~ prin_comp, CV=T)
ct4 <-table(Region, model4$class)
correct_assign4 <- as.data.frame(cbind(diag(prop.table(ct4, 1)), colnames(prop.table(ct4, 1))))
accuracy4<- c("model accuracy", sum(diag(prop.table(ct4))))
lda_summary4 <- cbind(eigvec, model4$posterior)

# Write LDA statistics
sink(paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, "_LDA_accuracy_sampling_locations.txt", sep=""))
print(correct_assign3)
print(accuracy3)
sink()
sink(paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, "_LDA_accuracy_sampling_regions.txt", sep=""))
print(correct_assign4)
print(accuracy4)
sink()
write.table(lda_summary3, paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, "_summary_LDA_sampling_locations.txt", sep=""),
row.names=F, col.names=T, quote=F)
write.table(lda_summary4, paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, "_summary_LDA_regions.txt", sep=""), row.names=F, col.names=T, quote=F)
q()