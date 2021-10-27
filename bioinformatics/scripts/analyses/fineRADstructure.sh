#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 40:00:00
#SBATCH -J fineRADstructure
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters

## Paired or SE reads
reads=SE; # "paired" or "SE"

## General
MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
WORKDIR=$MAINDIR/data/radseq/radseq2020/vcfs_final\
/$reads\_vcfs/non_pruned;
OUTDIR=$MAINDIR/working/Simon;
PLOT=$MAINDIR/scripts/radseq2020/plotting;
fineRAD=$MAINDIR/bin/radseq2020/fineRADstructure;

### Load modules
module load bioinfo-tools fineSTRUCTURE/4.0.1 R_packages/4.0.0;

### Create folders
mkdir $OUTDIR/fineRADstructure.$reads \
$OUTDIR/fineRADstructure.$reads/all \
$OUTDIR/fineRADstructure.$reads/transect \
$OUTDIR/fineRADstructure.$reads/locations;

### Perform finestructure analyses
OUTDIR=$OUTDIR/fineRADstructure.$reads;
DATASETS=$(echo all locations transect);
VCF=bark_beetles;

for DATA in ${DATASETS[@]}; do
## Convert input from VCF file
$fineRAD/RADpainter hapsFromVCF $WORKDIR/$VCF\_$DATA.vcf.gz \
> $OUTDIR/$DATA/$VCF\_$DATA.txt;

## Calculate the co-ancestry matrix
$fineRAD/RADpainter paint $OUTDIR/$DATA/$VCF\_$DATA.txt;

## Assign individuals to populations
fs fs -X -Y -x 10000000 -y 10000000 -z 10000 \
$OUTDIR/$DATA/$VCF\_$DATA\_chunks.out \
$OUTDIR/$DATA/$VCF\_$DATA\_chunks.mcmc.xml;

## Tree building
fs fs -X -Y -m T -x 1000000 \
$OUTDIR/$DATA/$VCF\_$DATA\_chunks.out \
$OUTDIR/$DATA/$VCF\_$DATA\_chunks.mcmc.xml \
$OUTDIR/$DATA/$VCF\_$DATA\_chunks.mcmcTree.xml;

## Run R-script to plot heatmap
Rscript $PLOT/fineRADstructurePlot.R --args  OUTDIR=$OUTDIR \
VCF=$VCF DATA=$DATA fineRAD=$fineRAD;
done;
