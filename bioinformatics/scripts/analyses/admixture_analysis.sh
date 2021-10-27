#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 3:00:00
#SBATCH -J admixture_analysis
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
module load bioinfo-tools ADMIXTURE/1.3.0 plink/1.90b4.9 \
gnuparallel/20180822 R_packages/4.0.0;

### Create folders
mkdir $OUTDIR/admixture.$reads \
$OUTDIR/admixture.$reads/all \
$OUTDIR/admixture.$reads/transect \
$OUTDIR/admixture.$reads/locations;

### Perform admixture analyses
OUTDIR=$OUTDIR/admixture.$reads;
DATASETS=$(echo all locations transect);
VCF=bark_beetles;

for DATA in ${DATASETS[@]}; do
## Create bed files
plink --vcf $WORKDIR/$VCF\_$DATA.pruned.vcf.gz --double-id --allow-extra-chr \
--set-missing-var-ids @:# --make-bed --out $OUTDIR/$DATA/$VCF\_$DATA;
awk '{$1=0;print $0}' $OUTDIR/$DATA/$VCF\_$DATA.bim \
> $OUTDIR/$DATA/$VCF\_$DATA.bim.tmp;
mv $OUTDIR/$DATA/$VCF\_$DATA.bim.tmp $OUTDIR/$DATA/$VCF\_$DATA.bim;

## Run admixture
K=$(seq 1 1 $(cat $METADATA/popmap_$DATA.txt | cut -f 2 | sort | uniq | wc -l) | sort -n -r);
cd $OUTDIR/$DATA;
export OUTDIR DATA VCF;
parallel -k 'admixture $OUTDIR/$DATA/$VCF\_$DATA.bed {} --cv=100 -s $RANDOM \
> $OUTDIR/$DATA/$VCF\_$DATA.{}' ::: $K;
cat $OUTDIR/$DATA/$VCF\_$DATA* | grep "CV error" | cut -d' ' -f 3,4 | \
sed -e 's/(K=\(.*\)):/\1/' | sort -n > $OUTDIR/$DATA/$VCF\_$DATA.CV;

## Run R-script to plot admixture crossvalidation, and barplot for the best K
Rscript $PLOT/admixture_plot.r --args OUTDIR=$OUTDIR VCF=$VCF \
DATA=$DATA METADATA=$METADATA;
done;
