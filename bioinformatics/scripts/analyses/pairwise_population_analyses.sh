#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 02:00:00
#SBATCH -J pairwise_population_analyses
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
## Number of top principal components to use for eucledian distances
PCs=4;

## Paired or SE reads
reads=SE; # "paired" or "SE"

## General
MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
WORKDIR=$MAINDIR/data/radseq/radseq2020/vcfs_final\
/$reads\_vcfs/non_pruned;
OUTDIR=$MAINDIR/working/Simon;
METADATA=$MAINDIR/data/radseq/radseq2020/metadata;
PLOT=$MAINDIR/scripts/radseq2020/plotting;
genomics=$MAINDIR/bin/radseq2020/genomics_general;
PCA=$MAINDIR/results/radseq2020/PCA.$reads;

### Load modules
module load bioinfo-tools vcftools/0.1.16 python/3.7.2 \
bcftools/1.10 gnuparallel/20180822 R_packages/4.0.0;

### Create folders
mkdir $OUTDIR/pairwise_populations.$reads \
$OUTDIR/pairwise_populations.$reads/all \
$OUTDIR/pairwise_populations.$reads/transect \
$OUTDIR/pairwise_populations.$reads/locations;

### Perform pairwise population analyses
OUTDIR=$OUTDIR/pairwise_populations.$reads;
DATASETS=$(echo all locations transect);
VCF=bark_beetles;
POPS=$(cat $METADATA/popmap_all.txt | cut -f 2 | sort | uniq);
for pop in ${POPS[@]}; do
cat $METADATA/sample_info_all.txt | cut -f 1,3 | awk -v var="$pop" '$2~var{print}' > $OUTDIR/$pop;
done;
export WORKDIR OUTDIR VCF genomics;

for DATA in ${DATASETS[@]}; do
popmap=$METADATA/popmap_$DATA.txt;
export DATA popmap;
pop=$(cat $METADATA/order_$DATA.txt)
pop1=$(seq 1 1 $(echo $(cat $METADATA/popmap_$DATA.txt | cut -f 2 | sort | uniq | wc -l) -1 | bc));

## Make a list of pairwise comparisons
for i in ${pop1[@]}; do
pop2=$(seq  $(echo $(echo $i) +1 | bc) 1 $(echo $(echo "$pop1" | wc -l) +1 | bc));
for j in ${pop2[@]}; do
echo $(echo $pop | cut -d' ' -f$i) $(echo $pop | cut -d' ' -f$j) 
done;
done > $OUTDIR/$DATA/pairwise_list_$DATA.txt;
pairs=$OUTDIR/$DATA/pairwise_list_$DATA.txt;
pairnr=$(seq 1 1 $(cat $pairs | wc -l));

## Pairwise Fst
mkdir $OUTDIR/$DATA/Fst;
fst() {
p1=$1;
p2=$2;
vcftools --gzvcf $WORKDIR/$VCF\_$DATA.vcf.gz  --weir-fst-pop $OUTDIR/$p1 \
--weir-fst-pop $OUTDIR/$p2 --fst-window-size 100000 --fst-window-step 25000 \
--out $OUTDIR/$DATA/Fst/$VCF\_$DATA.$p1\_$p2 &> $OUTDIR/$DATA/Fst/$p1\_$p2;
echo $p1 $p2 $(cat $OUTDIR/$DATA/Fst/$p1\_$p2 | grep "Cockerham" | tail -n +2 | cut -d' ' -f7);
};
export -f fst;
parallel --colsep ' ' 'fst {1} {2}' :::: $pairs > $OUTDIR/$DATA/Fst/Fst_summary_$DATA.txt;

## Pairwise Dxy
mkdir $OUTDIR/$DATA/Dxy;
python3 $genomics/VCF_processing/parseVCF.py -i $WORKDIR/$VCF\_$DATA.vcf.gz | \
bgzip > $OUTDIR/$DATA/Dxy/parsed.$DATA.gz;
dxy() {
p1=$1;
p2=$2;
python3 $genomics/popgenWindows.py -g $OUTDIR/$DATA/Dxy/parsed.$DATA.gz \
-o $OUTDIR/$DATA/Dxy/$p1\_$p2 -p $p1 -p $p2 -f phased -w 100000 -s 25000 --popsFile $popmap;
echo $p1 $p2 $(cat $OUTDIR/$DATA/Dxy/$p1\_$p2 | tail -n +2 | cut -d',' -f8 | \
awk '{ total += $1; count++ } END { print total/count }');
};
export -f dxy;
parallel --colsep ' ' 'dxy {1} {2}' :::: $pairs > $OUTDIR/$DATA/Dxy/Dxy_summary_$DATA.txt;

## Summary of pairwise data
for i in  ${pairnr[@]}; do
echo $(cat $pairs | head -n $i | tail -n 1) \
$(cat $OUTDIR/$DATA/Fst/Fst_summary_$DATA.txt | grep "$(cat $pairs | head -n $i | tail -n 1)" | cut -d' ' -f3,4) \
$(cat $OUTDIR/$DATA/Dxy/Dxy_summary_$DATA.txt | grep "$(cat $pairs | head -n $i | tail -n 1)" | cut -d' ' -f3) \
$(cat $METADATA/coordinates.txt | grep "$(cat $pairs | head -n $i | tail -n 1 | cut -d' ' -f1)" | cut -f2,3) \
$(cat $METADATA/coordinates.txt | grep "$(cat $pairs | head -n $i | tail -n 1 | cut -d' ' -f2)" | cut -f2,3);
done > $OUTDIR/$DATA/pairwise_summary.temp.txt;
echo -e "location1 location2 mean_Fst weighted_Fst Dxy \
loc1_latitude loc1_longitude loc2_latitude loc2_longitude" | \
cat - $OUTDIR/$DATA/pairwise_summary.temp.txt > $OUTDIR/$DATA/pairwise_summary_$DATA.txt;
rm $OUTDIR/$DATA/pairwise_summary.temp.txt;

## Plot pairwise heatmaps and perform mantel tests
Rscript $PLOT/pairwise_population_analyses_plot.r --args  PCs=$PCs OUTDIR=$OUTDIR \
VCF=$VCF DATA=$DATA METADATA=$METADATA PCA=$PCA;
done;
