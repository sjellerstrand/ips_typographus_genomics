#!/usr/bin/Rscript

### Plot genome scans
## Export variables and load libraries
rm(list=ls())
library(tidyverse)
library(gridExtra)

args <- strsplit(commandArgs(trailingOnly=T), "=", fixed=T)
for(i in 1:length(args)) {
  assign(args[[i]][1], args[[i]][2])
}

## Import and modify data by adding window mid point and absolute position in reference genome

# List variables
data_list <- rbind(c("Fst", "WEIGHTED_FST"), c("Dxy", "dxy_Norrbotten_South"), c("pi", "PI"),
c("obs_het", "Observed_heterozygosity"), c("TajimasD", "TajimaD"))

# Reference genome
ref <-  read.delim(paste(REF, ".fai", sep=""), sep="\t", head=F)
colnames(ref) <- c("contig", "length", "abs.pos","NA", "NA")

# Fst
Fst <- read.delim(paste(OUTDIR, "/", DATA, "/Fst/", VCF, "_", DATA, ".Norrbotten_South.windowed.weir.fst", sep=""), sep="\t", head=T)
Fst <- cbind(Fst, matrix(NA, nrow(Fst), 3))
colnames(Fst)[(ncol(Fst)-2):ncol(Fst)] <- c("mid.pos", "abs.pos", "outlier")
for(i in 1:nrow(Fst)) {
      if(Fst$CHROM[i] != Fst$CHROM[i+1] || i+1 > nrow(Fst)) {
          if(ref$length[which(ref$contig == Fst$CHROM[i])] < Fst$BIN_END[i]-Fst$BIN_START[i]+1) {
               Fst$mid.pos[i] <- Fst$BIN_START[i] + round((ref$length[which(ref$contig == Fst$CHROM[i])] - Fst$BIN_START[i])/2)
          } else {
               Fst$mid.pos[i] <- Fst$BIN_START[i] + round((Fst$BIN_END[i] - Fst$BIN_START[i]+1)/2)
          }
    } else {
          Fst$mid.pos[i] <- Fst$BIN_START[i] + round((Fst$BIN_END[i] - Fst$BIN_START[i]+1)/2)
    }
    Fst$abs.pos[i] <- ref$abs.pos[which(ref$contig == Fst$CHROM[i])] + Fst$mid.pos[i]
}

# Absolute divergence Dxy
Dxy <- read.delim(paste(OUTDIR, "/", DATA, "/Dxy/Norrbotten_South", sep=""), sep=",", head=T)
Dxy<- cbind(Dxy, matrix(NA, nrow(Dxy), 3))
colnames(Dxy)[(ncol(Dxy)-2):ncol(Dxy)] <- c("mid.pos", "abs.pos", "outlier")
for(i in 1:nrow(Dxy)) {
     if(Dxy$scaffold[i] != Dxy$scaffold[i+1] || i+1 > nrow(Dxy)) {
          if(ref$length[which(ref$contig == Dxy$scaffold[i])] < Dxy$end[i]-Dxy$start[i]+1) {
               Dxy$mid.pos[i] <- Dxy$start[i] + round((ref$length[which(ref$contig == Dxy$scaffold[i])] - Dxy$start[i])/2)
          } else {
               Dxy$mid.pos[i] <- Dxy$start[i] + round((Dxy$end[i] - Dxy$start[i]+1)/2)
          }
    } else {
        Dxy$mid.pos[i] <- Dxy$start[i] + round((Dxy$end[i] - Dxy$start[i]+1)/2)
    }
    Dxy$abs.pos[i] <- ref$abs.pos[which(ref$contig == Dxy$scaffold[i])] + Dxy$mid.pos[i]
}

# Nucleotide diversity pi
pi <- read.delim(paste(OUTDIR, "/", DATA, "/Nucleotide_diversity/global.windowed.pi", sep=""), sep="\t", head=T)
pi <- cbind(pi, matrix(NA, nrow(pi), 3))
colnames(pi)[(ncol(pi)-2):ncol(pi)] <- c("mid.pos", "abs.pos", "outlier")
for(i in 1:nrow(pi)) {
    if(pi$CHROM[i] != pi$CHROM[i+1] || i+1 > nrow(pi)) {
          if(ref$length[which(ref$contig == pi$CHROM[i])] < pi$BIN_END[i]-pi$BIN_START[i]+1) {
               pi$mid.pos[i] <- pi$BIN_START[i] + round((ref$length[which(ref$contig == pi$CHROM[i])] - pi$BIN_START[i])/2)
          } else {
               pi$mid.pos[i] <- pi$BIN_START[i] + round((pi$BIN_END[i] - pi$BIN_START[i]+1)/2)
          }
    } else {
        pi$mid.pos[i] <- pi$BIN_START[i] + round((pi$BIN_END[i] - pi$BIN_START[i]+1)/2)
    }
    pi$abs.pos[i] <- ref$abs.pos[which(ref$contig == pi$CHROM[i])] + pi$mid.pos[i]
}

# Observed heterozygosity
obs_het <- read.delim(paste(OUTDIR, "/", DATA, "/Observed_heterozygosity/global2.hwe", sep=""), sep="\t", head=F)
obs_het <- cbind(obs_het, matrix(NA, nrow(obs_het), 2))
colnames(obs_het) <- c("CHROM", "mid.pos", "Observed_heterozygosity", "abs.pos", "outlier")
for(i in 1:nrow(obs_het)) {
    obs_het$abs.pos[i] <- ref$abs.pos[which(ref$contig ==  obs_het$CHROM[i])] +  obs_het$mid.pos[i]
}

# Tajimas D
TajimasD <- read.delim(paste(OUTDIR, "/", DATA, "/TajimasD/global.Tajima.D", sep=""), sep="\t", head=T)
TajimasD  <- cbind(TajimasD , matrix(NA, nrow(TajimasD), 3))
colnames(TajimasD )[(ncol(TajimasD)-2):ncol(TajimasD )] <- c("mid.pos", "abs.pos", "outlier")
for(i in 1:nrow(TajimasD )) {
    if(TajimasD$CHROM[i] != TajimasD$CHROM[i+1] || i+1 > nrow(TajimasD)) {
        if(ref$length[which(ref$contig == TajimasD$CHROM[i])] < TajimasD$BIN_START[2] - TajimasD$BIN_START[1]) {
             TajimasD$mid.pos[i] <- TajimasD$BIN_START[i] + round((ref$length[which(ref$contig == TajimasD $CHROM[i])] - TajimasD$BIN_START[i])/2)
         } else {
             TajimasD$mid.pos[i] <- TajimasD$BIN_START[i] + round((TajimasD$BIN_START[2] - TajimasD$BIN_START[1])/2+1)
         }
    } else {
        TajimasD$mid.pos[i] <- TajimasD$BIN_START[i] + round((TajimasD$BIN_START[2] - TajimasD$BIN_START[1])/2+1)
    }
    TajimasD$abs.pos[i] <- ref$abs.pos[which(ref$contig == TajimasD $CHROM[i])] + TajimasD$mid.pos[i]
}


## Find outliers and make plots
for(i in 1:nrow(data_list)) {
    # Find outliers
    quan <- quantile(eval(parse(text=paste(data_list[i,1], "$", data_list[i,2], sep=""))) , c(0.975,0.995), na.rm = T)
    q95 <- as.data.frame(matrix(NA, nrow=1, ncol(eval(parse(text=paste(data_list[i,1]))))))
    colnames(q95) <- colnames(eval(parse(text=paste(data_list[i,1]))))
    q99 <- as.data.frame(matrix(NA, nrow=1, ncol(eval(parse(text=paste(data_list[i,1]))))))
    colnames(q99) <- colnames(eval(parse(text=paste(data_list[i,1]))))
    for(j in 1:nrow(eval(parse(text=paste(data_list[i,1]))))) {
          if(eval(parse(text=paste(data_list[i,1], "$", data_list[i,2], "[",j ,"]", sep=""))) < quan[1] ||
          is.na(eval(parse(text=paste(data_list[i,1], "$", data_list[i,2], "[",j ,"]", sep=""))))) {
               eval(parse(text=paste(data_list[i,1], "$outlier[",j ,'] <- "out_no"', sep="")))
          } else if(eval(parse(text=paste(data_list[i,1], "$", data_list[i,2], "[",j ,"]", sep=""))) < quan[2]) {
               eval(parse(text=paste(data_list[i,1], "$outlier[",j ,'] <- "out_95"', sep="")))
               q95[nrow(q95)+1,] <- eval(parse(text=paste(data_list[i,1], "[",j ,",]", sep="")))
          } else {
               eval(parse(text=paste(data_list[i,1], "$outlier[",j ,'] <- "out_99"', sep=""))) 
               q99[nrow(q99)+1,] <- eval(parse(text=paste(data_list[i,1], "[",j ,",]", sep="")))
          }
    }
    q95 <- q95[-1,]
    q99 <- q99[-1,]

    # Write outlier data
     write.table(q95, paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, "_", eval(parse(text=paste('"', data_list[i,2], '"'))),"_outliers_q95.txt", sep=""),
     row.names=F, col.names=T, quote=F)
     write.table(q99, paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, "_", eval(parse(text=paste('"', data_list[i,2], '"'))),"_outliers_q99.txt", sep=""),
     row.names=F, col.names=T, quote=F)


    # Create plots
    scan <- ggplot(eval(parse(text=paste(data_list[i,1]))), aes(abs.pos, eval(parse(text=paste(data_list[i,2]))), colour = outlier)) +
    geom_point(size = 7) + theme_light() + geom_hline(yintercept = quan[1], size = 1, colour = "grey", linetype = "dashed") +
    geom_hline(yintercept = quan[2], size = 1, colour = "grey", linetype = "dashed") +
    theme(panel.grid = element_blank(), panel.border = element_rect(colour = "black", fill=NA, size=1), legend.position = "none",
    axis.title.x=element_text(size=41), axis.text.x = element_blank(), axis.text.y = element_text(size=41, colour="black"),
    axis.title.y=element_text(size=41)) + ylab(eval(parse(text=paste('"', data_list[i,2], '"')))) + xlab("Reference genome") +
    scale_color_manual(values=c("orange", "red", "#595959"))
    scan <- ggExtra::ggMarginal(scan, type = "histogram", margins="y", size=8)

    # Write plot as .jpeg
    jpeg(paste(OUTDIR, "/", DATA, "/", VCF, "_", DATA, "_", eval(parse(text=paste('"', data_list[i,2], '"'))), "_genome_scan.jpeg", sep=""),
    width=3000, height=1000, quality=100)
    print(scan)
    dev.off()
}
q()
