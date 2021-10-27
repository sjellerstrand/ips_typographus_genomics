#!/usr/bin/Rscript

### Plot pairwise heatmaps and perform mantel tests
## Export variables and load libraries
rm(list=ls())
library(plyr)
library(tidyverse)
library(reshape2)
library(viridis)
library(raster)
library(ade4)

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
assign(args[[i]][1], args[[i]][2])
}

## Read in data
pairwise <- read.delim(paste(OUTDIR, "/", DATA, "/pairwise_summary_", DATA, ".txt", sep=""), sep=" ", header=T)
order <- as.matrix(read.delim(paste(METADATA, "/order_", DATA, ".txt", sep=""), header=F))
pop <- read.table(paste(METADATA, "/popmap_", DATA, ".txt", sep=""), sep="\t", header=F)
translate <- read.delim(paste(METADATA, "/translate_population_names.txt", sep=""), header=F)
prin_comp <- read.delim(paste(PCA, "/", DATA, "/bark_beetles_", DATA, ".eigenvec", sep=""), sep=" ", header=F)
prin_comp <- prin_comp[,-1]

## Calculate and add Eucledian distances from Principal components
prin_comp <- cbind(prin_comp, rep(NA, nrow(prin_comp)))
for(i in 1:nrow(prin_comp)) {
prin_comp[i,ncol(prin_comp)] <- pop[which(pop[,1] == prin_comp[i,1]), 2]
}
colnames(prin_comp)[ncol(prin_comp)] <- "Locality"
euc_mean <- cbind(order, as.data.frame(matrix(NA, length(order), as.numeric(PCs))))
for(i in 1:length(order)) {
for(j in 1:PCs) {
euc_mean[i,j+1] <- mean(prin_comp[which(prin_comp[,ncol(prin_comp)] == order[i]),j+1])
}
}
pairwise <- cbind(pairwise[,1:(ncol(pairwise)-4)], rep(NA, nrow(pairwise)), pairwise[,(ncol(pairwise)-3):ncol(pairwise)])
colnames(pairwise)[ncol(pairwise)-4] <- "Eucledian"
for(i in 1:nrow(pairwise)) {
pairwise[i,(ncol(pairwise)-4)] <- dist(rbind(euc_mean[which(euc_mean[,1] == pairwise[i,1]),2:ncol(euc_mean)], euc_mean[which(euc_mean[,1] == pairwise[i,2]),2:ncol(euc_mean)]))
}

## Add empty "diagonal" and mirror triangle
for(i in 1:length(order[,1])) {
pairwise <- rbind(pairwise, c(order[i,1], order[i,1], rep(NA, ncol(pairwise)-2)))
}

pairwise[,1] <- factor(pairwise[,1], levels=order[,1])
pairwise[,2] <- factor(pairwise[,2], levels=order[,1])
pairwise[,3:ncol(pairwise)] <- lapply(pairwise[,3:ncol(pairwise)], as.numeric)

## Create and write heatmaps for each pairwise variable
pairwise2 <- as.data.frame(rbind(as.matrix(pairwise), as.matrix(cbind(pairwise[,2], pairwise[,1], pairwise[,3:ncol(pairwise)]))))
pairwise2[,1] <- factor(pairwise2[,1], levels=order[,1])
pairwise2[,2] <- factor(pairwise2[,2], levels=order[,1])
pairwise2$location1 <- mapvalues(pairwise2$location1, from = translate[,1], to = translate[,2])
pairwise2$location2 <- mapvalues(pairwise2$location2, from = translate[,1], to = translate[,2])
pairwise2[,3:ncol(pairwise2)] <- lapply(pairwise2[,3:ncol(pairwise2)], as.numeric)
translate2 <- levels(pairwise2$location1)
translate2 <- cbind(translate2, seq(length(translate2), 1))
translate2 <- cbind(translate2, paste(translate2[,2], ". ", translate2[,1], sep=""))
pairwise2$location1 <- mapvalues(pairwise2$location1, from = translate2[,1], to = translate2[,3])
pairwise2$location2 <- mapvalues(pairwise2$location2, from = translate2[,1], to = translate2[,3])
for(i in 3:(ncol(pairwise2)-4)) {
# Create heatmap
min <-  min(pairwise2[,i], na.rm=T)
max <- max(pairwise2[,i], na.rm=T)

diff <- (max-min)/3
heatmap <- ggplot(data = pairwise2[,c(1,2,i)], aes_string("location1", "location2", fill = colnames(pairwise2)[i])) +
geom_tile(color = "white") + theme_minimal() + theme(panel.grid = element_blank(), axis.title.x = element_blank(),
axis.title.y=element_blank(), axis.text.y = element_text(hjust = 0, size=25, colour="black"),
axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, size=25, colour="black"), legend.title =element_text(size=30),
legend.text =element_text(size=25), legend.key.size =unit(1.5, "cm")) + coord_fixed() +
scale_fill_viridis(name=colnames(pairwise2)[i], discrete=F, na.value="white", option="plasma", direction = 1,
labels = c(min, signif(min+diff, digits=8), signif(max-diff, digits=8), max), breaks = c(min, signif(min+diff, digits=8), signif(max-diff, digits=8), max))

# Write plot as .jpeg
jpeg(paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, "_", colnames(pairwise2)[i], "_heatmap.jpeg", sep=""), width=1000, height=1000, quality=100)
print(heatmap)
dev.off()
}

## Remove Kinda
if(length(which(levels(pairwise[,1]) == "Kinda")) != 0)  {
pairwise <- pairwise[-which(pairwise[,1] == "Kinda"),]
}
if(length(which(levels(pairwise[,2]) == "Kinda")) != 0)  {
pairwise <- pairwise[-which(pairwise[,2] == "Kinda"),]
}

## Calculate geographic distances
pairwise <- cbind(pairwise, rep(NA, nrow(pairwise)))
colnames(pairwise)[ncol(pairwise)] <- "distance"
for(i in 1:nrow(pairwise)) {
pairwise[i, ncol(pairwise)] <- pointDistance(c(pairwise$loc1_longitude[i], pairwise$loc1_latitude[i]),
c(pairwise$loc2_longitude[i], pairwise$loc2_latitude[i]), lonlat=T)/1000
}
pairwise$distance <- log10(pairwise$distance)

# Transform data
geo_dist <- t(acast(pairwise, value.var="distance" ,formula=location1 ~ location2))
geo_dist<- as.dist(geo_dist)
pairwise$mean_Fst <- pairwise$mean_Fst/(1-pairwise$mean_Fst)
pairwise$weighted_Fst <- pairwise$weighted_Fst/(1-pairwise$weighted_Fst)
pairwise[,5:(ncol(pairwise)-5)] <- log10(pairwise[,5:(ncol(pairwise)-5)])

## Perform Mantel tests
for(i in 3:(ncol(pairwise)-5)) {
# Create distance matrix
assign(paste(colnames(pairwise)[i], "_dist", sep=""), t(acast(pairwise, value.var=colnames(pairwise)[i] ,formula=location1 ~ location2)))
assign(paste(colnames(pairwise)[i], "_dist", sep=""), as.dist(eval(parse(text=paste(colnames(pairwise)[i], "_dist", sep="")))))

# Perform Mantel test
mantel_test <- mantel.randtest(geo_dist, eval(parse(text=paste(colnames(pairwise)[i], "_dist", sep=""))), nrepet=1000000)

# Plot Mantel test
IBD <- ggplot(data = pairwise, aes_string("distance", colnames(pairwise)[i])) + geom_point(size = 5, colour="#595959") + theme_light() +
theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), axis.title.x=element_text(size=30), axis.title.y=element_text(size=30), axis.text.y=element_text(size=30, colour="black"),
axis.text.x=element_text(size=30, colour="black")) + xlab("Geographic distance [log10 km]") + ylab(colnames(pairwise)[i]) +
geom_smooth(method = "lm", se=F, col="black") + scale_color_manual(values="#595959")
IBD <- ggExtra::ggMarginal(IBD, type = "histogram")

# Write plots
jpeg(paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, "_", colnames(pairwise)[i], "_IBD.jpeg", sep=""), width=1000, height=1000, quality=100)
print(IBD)
dev.off()

# Write Mantel test statistics
sink(paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, "_", colnames(pairwise)[i], "_manteltest.txt", sep=""))
print(mantel_test)
sink()
}
q()
