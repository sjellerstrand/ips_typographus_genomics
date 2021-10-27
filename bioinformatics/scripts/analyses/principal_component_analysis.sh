#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -J principal_component_analysis
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters

## Paired or SE reads
reads=SE; # "paired" or "SE"

## General
MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
WORKDIR=$MAINDIR/data/radseq/radseq2020/vcfs_final\
/$reads\_vcfs/pruned;
OUTDIR=$MAINDIR/working/Simon;
METADATA=$MAINDIR/data/radseq/radseq2020/metadata;
PLOT=$MAINDIR/scripts/radseq2020/plotting;

### Load modules
module load bioinfo-tools R_packages/4.0.0 plink/1.90b4.9;

### Create folders
mkdir $OUTDIR/PCA.$reads \
$OUTDIR/PCA.$reads/all \
$OUTDIR/PCA.$reads/transect \
$OUTDIR/PCA.$reads/locations;

### Perform principal component analyses
OUTDIR=$OUTDIR/PCA.$reads;
DATASETS=$(echo all locations transect);
VCF=bark_beetles;

for DATA in ${DATASETS[@]}; do
## Calculate Pricipal components
plink --vcf $WORKDIR/$VCF\_$DATA.pruned.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# --pca 20 \
--out $OUTDIR/$DATA/$VCF\_$DATA;

## Run R-script to plot pricipal component 1 and 2
Rscript $PLOT/pca_plot.r --args OUTDIR=$OUTDIR VCF=$VCF \
DATA=$DATA METADATA=$METADATA;
done;
