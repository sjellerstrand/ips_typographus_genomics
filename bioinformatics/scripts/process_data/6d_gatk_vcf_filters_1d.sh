#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 1:00:00
#SBATCH -J gatk_vcf_filters_1d
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters

## Paired or SE reads
reads=SE; # "paired" or "SE"

## General

MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
OUTDIR0=$MAINDIR/working/Simon/gatk_vcf_filters_1.$reads;
VCF_IN=$OUTDIR0/filter6/bark_beetles_filter6;
REF=$MAINDIR/data/reference/Itypographus_ref_customized\
/reference_customized.fasta;
METADATA=$MAINDIR/data/radseq/radseq2020/metadata;
functions=$MAINDIR/scripts/radseq2020/quality_control;
allelic_balance=$MAINDIR/bin/radseq2020/vcf_filters;

### Load modules
module load bioinfo-tools vcftools/0.1.16 vcflib/1.0.1 python/3.7.2 \
bcftools/1.10 R_packages/4.0.0 plink/1.90b4.9;

### Define functions
filter_stats=$functions/filter_stats.sh;
pca=$functions/pca.sh;
replicate_comp=$functions/replicate_comp.sh;
contam=$functions/contam.sh;

### Apply filters
DATA=2;
replicates=$METADATA/replicates_$DATA.txt;
sample_info=$METADATA/sample_info_$DATA.txt;
popmap=$METADATA/popmap_$DATA.txt;
remove=$METADATA/remove_inds.1b.txt;

## Filter 7: Remove troublesome individuals and filter low and high depth variants
filter=7;
mkdir $OUTDIR0/filter$filter;
OUTDIR=$OUTDIR0/filter$filter;
VCF_OUT=$OUTDIR/bark_beetles_filter$filter;
vcftools --gzvcf $VCF_IN.vcf.gz --remove $remove --max-missing 0.8 \
--min-meanDP 48 --minDP 30 --max-meanDP 52 --mac 3 \
--recode --recode-INFO-all --stdout | bgzip -c > $VCF_OUT.vcf.gz;
tabix $VCF_OUT.vcf.gz;
source $filter_stats;
source $pca;
source $replicate_comp;

## Filter 8: Allelic balance per heterozygous site and individual
VCF_IN=$VCF_OUT;
filter=8;
mkdir $OUTDIR0/filter$filter;
OUTDIR=$OUTDIR0/filter$filter;
VCF_OUT=$OUTDIR/bark_beetles_filter$filter;
python3 $allelic_balance/allelicBalance_SJE.py -i $VCF_IN.vcf \
-o $VCF_OUT.temp.vcf -r1 0.14 -r2 0.28;
vcftools --vcf $VCF_OUT.temp.vcf --mac 3 --recode \
--recode-INFO-all --stdout | bgzip -c > $VCF_OUT.vcf.gz;
rm $VCF_OUT.temp.vcf;
source $filter_stats;
source $pca;
source $replicate_comp;
source $contam;
