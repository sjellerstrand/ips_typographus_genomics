#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J population_statistics
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters

## Paired or SE reads
reads=SE; # "paired" or "SE"

## General
MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
WORKDIR1=$MAINDIR/data/radseq/radseq2020/vcfs_final/$reads\_vcfs/non_pruned;
WORKDIR2=$MAINDIR/data/radseq/radseq2020/vcfs_final/$reads\_vcfs/pruned;
OUTDIR=$MAINDIR/working/Simon;
METADATA=$MAINDIR/data/radseq/radseq2020/metadata;
PLOT=$MAINDIR/scripts/radseq2020/plotting;

### Load modules
module load bioinfo-tools vcftools/0.1.16 Stacks/2.53 \
gnuparallel/20180822 R_packages/4.0.0;

### Create folders
mkdir $OUTDIR/population_statistics.$reads \
$OUTDIR/population_statistics.$reads/all \
$OUTDIR/population_statistics.$reads/transect \
$OUTDIR/population_statistics.$reads/locations;

### Calculate population statistics
OUTDIR=$OUTDIR/population_statistics.$reads;
DATASETS=$(echo all locations transect);
VCF=bark_beetles;
POPS=$(cat $METADATA/popmap_all.txt | cut -f 2 | sort | uniq);
for pop in ${POPS[@]}; do
cat $METADATA/sample_info_all.txt | cut -f 1,3 | awk -v var="$pop" '$2~var{print}' > $OUTDIR/$pop;
done;
export WORKDIR1 WORKDIR2 OUTDIR VCF;

for DATA in ${DATASETS[@]}; do
export DATA;
POPS=$(cat $METADATA/order_$DATA.txt);

## Tajima's D
mkdir $OUTDIR/$DATA/TajimasD;
tajimasd() {
pop=$1;
vcftools --gzvcf $WORKDIR1/$VCF\_$DATA.vcf.gz --keep $OUTDIR/$pop \
--TajimaD 100000 --out $OUTDIR/$DATA/TajimasD/$pop;
};
export -f tajimasd;
parallel 'tajimasd {}' ::: $POPS 

## Nucleotide diversity pi
mkdir $OUTDIR/$DATA/Nucleotide_diversity;
pi() {
pop=$1;
vcftools --gzvcf $WORKDIR1/$VCF\_$DATA.vcf.gz --keep $OUTDIR/$pop \
--window-pi 100000 --window-pi-step 25000 \
--out $OUTDIR/$DATA/Nucleotide_diversity/$pop;
echo $pop $(cat $OUTDIR/$DATA/Nucleotide_diversity/$pop.windowed.pi | tail -n +2 | \
awk '{ total += $5; count++ } END { print total/count }');
};
export -f pi;
parallel 'pi {}' ::: $POPS > $OUTDIR/$DATA/pi_summary_$DATA.txt;

## Observed heterozygosity
mkdir $OUTDIR/$DATA/Observed_heterozygosity;
obs_het() {
pop=$1;
vcftools --gzvcf $WORKDIR2/$VCF\_$DATA.pruned.vcf.gz --keep $OUTDIR/$pop \
--hardy --out $OUTDIR/$DATA/Observed_heterozygosity/$pop;
echo $pop $(cat $OUTDIR/$DATA/Observed_heterozygosity/$pop.hwe | tail -n +2 | \
cut -f1,2,3 | awk '{ print $1":"$2"\t"$3}' | awk '{gsub("/", "\t"); print}' | \
awk '{sumALL[$1] += $2+$3+$4; sumHET[$1] += $3}; \
END{for (id in sumALL) {print id, "\t", sumHET[id]/sumALL[id]}}' | \
awk '{ total += $2; count++ } END { print total/count }');
};
export -f obs_het;
parallel 'obs_het {}' ::: $POPS > $OUTDIR/$DATA/obs_het_summary_$DATA.txt;

## Private alleles
mkdir $OUTDIR/$DATA/Private_alleles;
populations -V $WORKDIR1/$VCF\_$DATA.vcf.gz -O $OUTDIR/$DATA/Private_alleles \
--popmap $METADATA/popmap_$DATA.txt;
cat $OUTDIR/$DATA/Private_alleles/$VCF\_$DATA.p.sumstats_summary.tsv | cut -f1,2 | \
tail -n $(echo $(cat $OUTDIR/$DATA/Private_alleles/$VCF\_$DATA.p.sumstats_summary.tsv | wc -l) / 2 -2 | bc) \
> $OUTDIR/$DATA/priv_allele_summary_$DATA.txt;

## Summary of population data
for pop in  ${POPS[@]}; do
echo $pop \
$(cat $OUTDIR/$DATA/pi_summary_$DATA.txt | grep $pop | cut -d' ' -f2) \
$(cat $OUTDIR/$DATA/obs_het_summary_$DATA.txt | grep $pop | cut -d' ' -f2) \
$(cat $OUTDIR/$DATA/priv_allele_summary_$DATA.txt | grep $pop | cut -f2);
done > $OUTDIR/$DATA/population_summary.temp.txt;
echo -e "location Nucleotide_diversity Observed_heterozygosity Private_alleles" | \
cat - $OUTDIR/$DATA/population_summary.temp.txt > $OUTDIR/$DATA/population_summary_$DATA.txt;
rm $OUTDIR/$DATA/population_summary.temp.txt;

## Plot population statistics and Tajimas D probability density plots
Rscript $PLOT/population_statistics_plot.r --args  OUTDIR=$OUTDIR \
VCF=$VCF DATA=$DATA METADATA=$METADATA;
done;
