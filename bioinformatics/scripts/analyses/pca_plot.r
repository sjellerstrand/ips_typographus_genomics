#!/usr/bin/Rscript

### Plot Pricipal component 1 and 2
## Export variables and load libraries
rm(list=ls())
library(plyr)
library(tidyverse)
library(viridis)

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
  assign(args[[i]][1], args[[i]][2])
}

## read in data
pca <- read.delim(paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, ".eigenvec", sep=""), sep=" ", head=F)
eigenval <- scan(paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, ".eigenval", sep=""))
pca <- pca[,-1]

## set names
names(pca)[1] <- "ind"
names(pca)[2:ncol(pca)] <- paste0("PC", 1:(ncol(pca)-1))

## sort out the individual sample locations and geographic regions
info <- read.delim(paste(METADATA, "/sample_info_", DATA, ".txt", sep=""), header=T)
order <- read.delim(paste(METADATA, "/order_", DATA, ".txt", sep=""), header=F)
translate <- read.delim(paste(METADATA, "/translate_population_names.txt", sep=""), header=F)
geo <- matrix(NA, nrow(pca), 5)
for(i in 1:nrow(pca)) {
geo[i,1:5] <-as.character(c(pca[i,1], info[which(pca[i,1] == info[,1]),c(1, 3, 4, 5)]))
}
geo <- as.data.frame(geo)
geo[,3] <- factor(geo[,3], levels=rev(order[,1]))
colnames(geo) <- c("ID", "ID", "Locality", "Region", "Library")
pca <- cbind(pca, geo)
pca$Locality <- mapvalues(pca$Locality, from = translate[,1], to =  translate[,2])

## calculate the cumulative sum of the percentage variance explained
pve <- data.frame(PC = 1:length(eigenval), pve = eigenval/sum(eigenval)*100)
cumsum(pve$pve)

## Create plots
pca <- as_tibble(data.frame(pca, geo[,4]))
a <- ggplot(pve, aes(PC, pve)) + geom_bar(stat = "identity") + ylab("Percentage variance explained") + theme_light()
b <- ggplot(pca, aes(PC1, PC2), col)  + geom_point(size = 8, aes(color=Locality), ) + stat_ellipse(aes(color=Locality)) + 
coord_equal() + theme_light() + xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
theme(panel.grid = element_blank(), axis.title.x = element_text(size=25), axis.title.y = element_text(size=25),
axis.text.x = element_text(size=25, colour="black"), axis.text.y = element_text(size=25, colour="black"),
panel.border = element_rect(colour = "black", fill=NA, size=1), legend.key.size = unit(1.5, "cm")) +
scale_color_viridis(discrete = T, guide = guide_legend(override.aes = list(size = 5, linetype=rep("blank", length(unique(pca$Locality))))))

## Write plots as .jpeg
jpeg(paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, "_eigenval.jpeg", sep=""), width=500, height=500, quality=100)
print(a)
dev.off()
jpeg(paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, "_PCA.jpeg", sep=""), width=1000, height=1000, quality=100)
print(b)
dev.off()
q()