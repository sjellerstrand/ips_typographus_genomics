#!/usr/bin/Rscript

### Plot replicate comparisons
## Export variables and load libraries
rm(list=ls())
library(tidyverse)
library(ggplot2)
library(gridExtra)

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
  assign(args[[i]][1], args[[i]][2])
}

## read in data
err <- read.delim(paste(VCF_OUT, ".err", sep=""), header=T)
rep <- read.delim(replicates, header=F, sep=" ")
rep2 = rbind(rep, cbind(rep[,2], rep[,1]))
info <- read.delim(sample_info, header=T)
list <- list()
disc <- data.frame(matrix(NA, nrow(rep2), 3))

## Create plots
for(i in 1:nrow(rep2)) {
    data <- err[which(err[,4] == rep2[i,1] | err[,5] == rep2[i,1]) ,]
    data[6] <- rep("dataset", nrow(data))
    data[which(data[,4] == rep2[which(rep2[,1] == rep2[i,1]), 2] | data[,5] == rep2[which(rep2[,1] == rep2[i,1]), 2]), 6] <- "replicate"
    colnames(data)[c(2,6)] <- c("diff", "group")
    data[,2] <- as.numeric(data[,2])
    list[[i]] <- ggplot(data) + theme_light() + ggtitle(rep2[i,1]) +
    geom_point(data=subset(data, group == 'dataset'), aes(x=rep[i,1], y=diff),
    color="green", size=4, position=position_jitter(w=0.1, h=0)) +
    geom_point(data=subset(data, group == 'replicate'), aes(x=rep[i,1], y=diff),
    color="red", size=4, position=position_jitter(w=0.1, h=0))
    disc[i,] <- data[which(data[,6] == "replicate"), c(4, 5, 2)]
}
disc <- rbind(disc[1:(nrow(disc)/2),], c("Average", "Error-rate", mean(disc[1:(nrow(disc)/2),3])))

## print error rates and plot replicate comparisons
write.table(disc, paste(VCF_OUT, "_replicates.txt", sep=""), row.names=F, col.names=F, quote=F)
jpeg(paste(VCF_OUT, "_replicates.jpeg", sep=""), width=2000, height=1300, quality=100)
grid.arrange(grobs=list, nrow=2, ncol=13)
dev.off()
