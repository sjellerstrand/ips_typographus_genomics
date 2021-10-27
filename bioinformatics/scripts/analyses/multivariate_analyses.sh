#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J multivariate_analyses
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters

## Number of top principal components to use
PCs=4;

## Paired or SE reads
reads=SE; # "paired" or "SE"

## General
MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
WORKDIR=$MAINDIR/results/radseq2020/PCA.$reads;
OUTDIR=$MAINDIR/working/Simon;
METADATA=$MAINDIR/data/radseq/radseq2020/metadata;
PLOT=$MAINDIR/scripts/radseq2020/plotting;
VCF=bark_beetles;

### Load modules
module load R_packages/4.0.0;

### Create folders
mkdir $OUTDIR/multivariate_analyses.$reads \
$OUTDIR/multivariate_analyses.$reads/all \
$OUTDIR/multivariate_analyses.$reads/transect \
$OUTDIR/multivariate_analyses.$reads/locations;

### Perform multivariate analyses with principal components
OUTDIR=$OUTDIR/multivariate_analyses.$reads;
DATASETS=$(echo all locations transect);

for DATA in ${DATASETS[@]}; do
Rscript $PLOT/multivariate_analyses_plot.r --args PCs=$PCs WORKDIR=$WORKDIR \
OUTDIR=$OUTDIR VCF=$VCF DATA=$DATA METADATA=$METADATA;
done;
