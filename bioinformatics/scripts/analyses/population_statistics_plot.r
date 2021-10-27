#!/usr/bin/Rscript

### Plot population statistics and Tajimas D probability density plots
## Export variables and load libraries
rm(list=ls())
library(plyr)
library(tidyverse)
library(reshape2)
library(viridis)
library(gridExtra)

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
  assign(args[[i]][1], args[[i]][2])
}

## Read in data
pop_stats<- read.delim(paste(OUTDIR, "/", DATA, "/population_summary_", DATA, ".txt", sep=""), sep=" ", header=T)
order <- as.matrix(read.delim(paste(METADATA, "/order_", DATA, ".txt", sep=""), header=F))
translate <- read.delim(paste(METADATA, "/translate_population_names.txt", sep=""), header=F)
pop_stats <- cbind(pop_stats, rep(NA, nrow(pop_stats)))
Tajima <- list()
for(i in 1:length(order)) {
     Tajima[[i]] <- read.delim(paste(OUTDIR, "/", DATA, "/TajimasD/", order[i], ".Tajima.D", sep=""), sep="\t", header=T)$Tajima
     names(Tajima)[i] <- order[i]
     pop_stats[which(pop_stats[,1] == order[i]),ncol(pop_stats)] <- median(as.numeric(Tajima[[i]]), na.rm=T)
}
colnames(pop_stats)[ncol(pop_stats)] <- "TajimasD"
pop_stats[,1] <- factor(pop_stats[,1], levels=order[,1])
pop_stats$location <- mapvalues(pop_stats$location, from = translate[,1], to = translate[,2])
translate2 <- levels(pop_stats$location)
translate2 <- cbind(translate2, seq(length(translate2), 1))
translate2 <- cbind(translate2, paste(translate2[,2], ". ", translate2[,1], sep=""))
pop_stats$location <- mapvalues(pop_stats$location, from = translate2[,1], to = translate2[,3])
pop_stats[,2:ncol(pop_stats)] <- lapply(pop_stats[2:ncol(pop_stats)], as.numeric)
output_stats <- list()

## Create global stats heatmaps
for(i in 2:ncol(pop_stats)) {
     # Create heatmap
     min <-  min(pop_stats[,i], na.rm=T)
     max <- max(pop_stats[,i], na.rm=T)
     diff <- (max-min)/3
     if(diff !=0) {
     heatmap <- ggplot(data = pop_stats[,c(1,i)], aes_string(1, "location", fill = colnames(pop_stats)[i])) +
     geom_tile(color = "white") + theme_minimal() + theme(panel.grid = element_blank(), axis.title.x = element_blank(),
     axis.title.y=element_blank(), axis.text.y = element_text(size=25, hjust=0, colour="black"), axis.text.x = element_blank(),
     legend.title =element_text(size=30), legend.text =element_text(size=25), legend.key.size =unit(1.5, "cm")) +
     coord_fixed() + scale_fill_viridis(name=colnames(pop_stats)[i], discrete=F, na.value="white", option="plasma",
     direction = 1, labels = c(min, signif(min+diff, digits=8), signif(max-diff, digits=8), max),
     breaks = c(min, signif(min+diff, digits=8), signif(max-diff, digits=8), max))

     # Write plot as .jpeg
     jpeg(paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, "_", colnames(pop_stats)[i], "_heatmap.jpeg", sep=""), width=1000, height=1000, quality=100)
     print(heatmap)
     dev.off()
}
}

## Plot Tajimas'D probability distributions
Tajima2 <- melt(Tajima)
colnames(Tajima2) <- c("data", "Locality")
Tajima2[,2] <- factor(Tajima2[,2], levels=rev(order[,1]))
Tajima2$Locality <- mapvalues(Tajima2$Locality, from = translate[,1], to = translate[,2])
Tajima2$Locality  <- mapvalues(Tajima2$Locality, from = translate2[,1], to = translate2[,3])
Tajima_plot <-
ggplot(Tajima2, aes(data), col) + geom_vline(xintercept = 0, size = 1, colour = "grey", linetype = "dashed") +
geom_density(aes(color=Locality), size=1) + theme_light() + theme(panel.grid = element_blank(),
axis.title.x=element_text(size=35), axis.text.x = element_text(size=30, colour="black"),
axis.text.y = element_text(size=30, colour="black"), panel.border = element_rect(colour = "black", fill=NA, size=1),
axis.title.y=element_text(size=35), legend.title =element_blank(), legend.text =element_text(size=25),
plot.margin = unit(c(2, 1, 1, 1), "cm"), legend.key.size = unit(1.1, "cm")) + xlab("Tajimas' D") + ylab("Density") +
scale_color_viridis(discrete = T)

## Write plots as .jpeg
jpeg(paste(OUTDIR, "/", DATA, "/Tajimas_D_", DATA, ".jpeg", sep=""), width=1000, height=800, quality=100)
print(Tajima_plot)
dev.off()
q()