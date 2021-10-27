#!/usr/bin/Rscript

### Plot statistics
## Export variables and load libraries
rm(list=ls())
library(tidyverse)
library (ggplot2)
library(dplyr)
library(gridExtra)

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
assign(args[[i]][1], args[[i]][2])
}

## Read in data
data <- read.delim(paste(VCF_OUT, ".gatk.annotations", sep=""), sep="\t", header=T)
a1<- ggplot(data, aes(QD)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + ggtitle("QD") + theme_light()
b1 <- summary(data$QD[1:5])
a2<- ggplot(data, aes(FS)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + ggtitle("FS") + theme_light()
b2 <- summary(data$FS[1:5])
a3<- ggplot(data, aes(MQ)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + ggtitle("MQ") + theme_light()
b3 <- summary(data$MQ[1:5])
a4<- ggplot(data, aes(MQRankSum)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + ggtitle("MQRankSum") + theme_light()
b4 <- summary(data$MQRankSum[1:5])
a5<- ggplot(data, aes(ReadPosRankSum)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + ggtitle("ReadPosRankSum") + theme_light()
b5 <- summary(data$ReadPosRankSum[1:5])
a6<- ggplot(data, aes(ExcessHet)) + geom_density(fill = "dodgerblue1", colour = "black", alpha = 0.3) + ggtitle("ExcessHet") + theme_light()
b6 <- summary(data$ExcessHet[1:5])

## Write plots as .jpeg
sumstat <-rbind(b1, b2, b3, b4, b5, b6)
sumstat <- cbind(c("QD", "FS", "MQ", "MQRankSum", "ReadPosRankSum", "ExcessHet"), sumstat)
tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))
tbl <- tableGrob(sumstat, rows=NULL, theme=tt)
jpeg(paste(VCF_OUT,"_gatk_filters.jpeg", sep=""), width=1300, height=1300, quality=100)
grid.arrange(a1, a2, a3, a4, a5, a6, tbl, nrow=4, ncol=2)
dev.off()
q()
