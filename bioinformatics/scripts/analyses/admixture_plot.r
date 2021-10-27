#!/usr/bin/Rscript

### Plot admixture crossvalidation, and barplot for the best K.
## Export variables and load libraries
rm(list=ls())
library(plyr)
library(tidyverse)
library(reshape2)
library(viridis)

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
  assign(args[[i]][1], args[[i]][2])
}

## Find best K
cv <- read.delim(paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, ".CV", sep=""), sep=" ", head=F)
best1 <- cv[which(cv[,2] == min(cv[,2])),1]
if(best1 == 1) {
    cv2 <- cv[-which(cv[,1] == "1"),]
    best <- cv2[which(cv2[,2] == min(cv2[,2])),1]
    second <- paste(", second best K=", best, sep="")
} else {
    best <- best1
    second <- ""
}
## Plot cross-validation and report best K
jpeg(paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, "_CV.jpeg", sep=""), width=500, height=500, quality=100)
plot(cv[,1], cv[,2], main=paste("cross-validation for admixture analysis\nbest K=", best1, second, sep=""), ylab="CV error", xlab="K")
points(best, cv[which(best == cv[,1]),2], col="red", pch=18, cex=3)
dev.off()

## Connect data with correct sample and sort according to population order
tbl <- read.table(paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, ".", best, ".Q", sep=""), sep=" ", head=F)
pop <- read.table(paste(METADATA, "/popmap_", DATA, ".txt", sep=""), sep="\t", header=F)
samples <- read.table(paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, ".fam", sep=""), sep=" ", head=F)
order <- read.delim(paste(METADATA, "/order_", DATA, ".txt", sep=""), header=F)
translate <- read.delim(paste(METADATA, "/translate_population_names.txt", sep=""), header=F)
geo <- matrix(NA, nrow(samples), 2)
for(i in 1:nrow(samples)) {
geo[i,1:2] <- as.character(pop[which(samples[i,1] == pop[,1]), 1:2])
}
tbl2 <- as.matrix(cbind(geo, tbl))
tbl3 <- matrix(NA, 1, ncol(tbl2))
popsum <- matrix(NA, nrow(order), best+1)
for(i in 1:nrow(order)) {
    tbl3 <- rbind(tbl3, tbl2[which(tbl2[,2] == order[i,1]),])
     pop <- as.data.frame(lapply(as.data.frame(tbl2[which(tbl2[,2] == order[i,1]), 3:(2+best)]), as.numeric))
    popsum[i,] <- c(order[i,1], colMeans(pop))
}
tbl4 <- as.data.frame(tbl3[-1,])
tbl5 <- cbind(tbl4[1:2], lapply(tbl4[,3:(2+best)], as.numeric))

## write population assignment summary
write.table(tbl5, paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, ".indsum", sep=""), sep="\t", row.names=F, col.names=F, quote=F)
popsum <- as.data.frame(popsum)
write.table(popsum, paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, ".popsum", sep=""), sep="\t", row.names=F, col.names=F, quote=F)

## Barplots for best K and assigment summary per individual and population
tbl6 <- cbind(lapply(as.data.frame(tbl5[,1:2]), factor), tbl5[,3:(2+best)])
colnames(tbl6)[1:2] <- c("id", "location")
tbl7 <- melt(tbl6, id=c("id", "location"))
tbl7[,2] <- factor(tbl7[,2], levels=order[,1])
tbl7$location <- mapvalues(tbl7$location, from = translate[,1], to = translate[,2])
tbl8 <- as_tibble(tbl7)
colnames(tbl8)[3:4] <- c("popK", "prob")

plot1 <-
  ggplot(tbl8, aes(factor(id), prob, fill = factor(popK))) +
  geom_col(color = "gray", size = 0.1) +
  facet_grid(~fct_inorder(location), switch = "x", scales = "free", space = "free") +
  theme_minimal() + labs(x = "Sampling locations", title = "K=2", y = "Ancestry") +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(expand = expansion(add = 1)) +
  theme(panel.spacing.x = unit(0.1, "lines"), axis.text.x = element_blank(), panel.grid = element_blank(),
  legend.position="none" ,strip.text.x = element_text(angle = 90, vjust = 0.5, hjust=1, size=40), 
  axis.title.x = element_text(size=50), axis.title.y = element_text(size=40), title = element_text(size=50), 
  axis.text.y = element_text(size=40)) + scale_fill_viridis_d()

translate2 <- levels(tbl7$location)
translate2 <- cbind(translate2, seq(length(translate2), 1))
translate2 <- cbind(translate2, paste(translate2[,2], ". ", translate2[,1], sep=""))
tbl7$location <- mapvalues(tbl7$location, from = translate2[,1], to = translate2[,3])
tbl8 <- as_tibble(tbl7)
colnames(tbl8)[3:4] <- c("popK", "prob")

plot2 <-
ggplot(tbl8, aes(factor(id), prob, fill = factor(popK))) +
geom_col(color = "gray", size = 0.1) +
facet_grid(~fct_inorder(location), switch = "x", scales = "free", space = "free") +
theme_minimal() + labs(x = "Sampling locations", title = "K=2", y = "Ancestry") +
scale_y_continuous(expand = c(0, 0)) +
scale_x_discrete(expand = expansion(add = 1)) +
theme(panel.spacing.x = unit(0.1, "lines"), axis.text.x = element_blank(), panel.grid = element_blank(),
legend.position="none" ,strip.text.x = element_text(angle = 270, vjust = 0.5, hjust=0, size=55), 
axis.title.x = element_blank(), axis.title.y = element_text(angle = 270, vjust = 0, hjust=0.5, size=60),
title = element_blank(), axis.text.y = element_text(angle = 270, vjust = 0.5, hjust=0.5, size=45, colour="black"), 
plot.margin = unit(c(3, 2, 2, 2), "cm")) +
scale_y_continuous(trans = "reverse") + scale_fill_viridis_d()

## Write plots as .jpeg
jpeg(paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, "_K", best, ".jpeg", sep=""), width=3000, height=1000, quality=100)
plot1
dev.off()
jpeg(paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, "_K", best, "_modified.jpeg", sep=""), width=3000, height=1000, quality=100)
plot2
dev.off()
q()
