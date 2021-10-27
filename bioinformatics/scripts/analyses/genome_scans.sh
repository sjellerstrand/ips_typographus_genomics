#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 02:00:00
#SBATCH -J genome_scans
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
REF=$MAINDIR/data/reference/Itypographus_ref_customized\
/reference_customized.fasta;
METADATA=$MAINDIR/data/radseq/radseq2020/metadata;
PLOT=$MAINDIR/scripts/radseq2020/plotting;
genomics=$MAINDIR/bin/radseq2020/genomics_general;

### Load modules
module load bioinfo-tools vcftools/0.1.16 python/3.7.2 \
bcftools/1.10 R_packages/4.0.0;

### Create folders
mkdir $OUTDIR/genome_scans.$reads \
$OUTDIR/genome_scans.$reads/locations;

### Perform genome scans
OUTDIR=$OUTDIR/genome_scans.$reads;
DATASETS=$(echo locations);
VCF=bark_beetles;

for DATA in ${DATASETS[@]}; do

cat $METADATA/sample_info_$DATA.txt | tail -n +2 | cut -f 1,4 | \
awk '$2 == "Norrbotten" {print $0}' > $OUTDIR/$DATA/Norrbotten;
cat $METADATA/sample_info_$DATA.txt | tail -n +2 | cut -f 1,4 | \
awk '$2 != "Norrbotten" {print $1,"\t","South"}' > $OUTDIR/$DATA/South;
cat  $OUTDIR/$DATA/Norrbotten  $OUTDIR/$DATA/South > $OUTDIR/$DATA/popmap.txt;
popmap=$OUTDIR/$DATA/popmap.txt;

## Global population statistics

# Global Tajima's D
mkdir $OUTDIR/$DATA/TajimasD;
vcftools --gzvcf $WORKDIR1/$VCF\_$DATA.vcf.gz \
--TajimaD 100000 --out $OUTDIR/$DATA/TajimasD/global;

# Global nucleotide diversity pi
mkdir $OUTDIR/$DATA/Nucleotide_diversity;
vcftools --gzvcf $WORKDIR1/$VCF\_$DATA.vcf.gz \
--window-pi 100000 --window-pi-step 25000 \
--out $OUTDIR/$DATA/Nucleotide_diversity/global;

# Global observed heterozygosity
mkdir $OUTDIR/$DATA/Observed_heterozygosity;
vcftools --gzvcf $WORKDIR2/$VCF\_$DATA.pruned.vcf.gz \
--hardy --out $OUTDIR/$DATA/Observed_heterozygosity/global;
cat $OUTDIR/$DATA/Observed_heterozygosity/global.hwe | tail -n +2 | \
cut -f1,2,3 | awk '{ print $1":"$2"\t"$3}' | awk '{gsub("/", "\t"); print}' | \
awk '{sumALL[$1] += $2+$3+$4; sumHET[$1] += $3}; \
END{for (id in sumALL) {print id, "\t", sumHET[id]/sumALL[id]}}' | \
awk '{gsub(":", "\t"); print}' | sort -V -k1,1 -k2,2 \
> $OUTDIR/$DATA/Observed_heterozygosity/global2.hwe

## North and South pairwise statistics

# Pairwise Fst
mkdir $OUTDIR/$DATA/Fst;
vcftools --gzvcf $WORKDIR1/$VCF\_$DATA.vcf.gz --weir-fst-pop $OUTDIR/$DATA/Norrbotten \
--weir-fst-pop $OUTDIR/$DATA/South --fst-window-size 100000 --fst-window-step 25000 \
--out $OUTDIR/$DATA/Fst/$VCF\_$DATA.Norrbotten_South &> $OUTDIR/$DATA/Fst/Norrbotten_South;
echo Norrbotten South $(cat $OUTDIR/$DATA/Fst/Norrbotten_South | grep "Cockerham" | tail -n +2 | cut -d' ' -f7) \
> $OUTDIR/$DATA/Fst/Fst_summary_$DATA.txt;

# Pairwise Dxy
mkdir $OUTDIR/$DATA/Dxy;
python3 $genomics/VCF_processing/parseVCF.py -i $WORKDIR1/$VCF\_$DATA.vcf.gz | \
bgzip > $OUTDIR/$DATA/Dxy/parsed.$DATA.gz;
python3 $genomics/popgenWindows.py -g $OUTDIR/$DATA/Dxy/parsed.$DATA.gz \
-o $OUTDIR/$DATA/Dxy/Norrbotten_South -p Norrbotten -p South -f phased \
-w 100000 -s 25000 --popsFile $popmap;
echo Norrbotten South $(cat $OUTDIR/$DATA/Dxy/Norrbotten_South | tail -n +2 | cut -d',' -f8 | \
awk '{ total += $1; count++ } END { print total/count }') > $OUTDIR/$DATA/Dxy/Dxy_summary_$DATA.txt;

## Plot genome scans
Rscript $PLOT/genome_scans_plot.r --args  OUTDIR=$OUTDIR \
REF=$REF VCF=$VCF DATA=$DATA METADATA=$METADATA;
done;
