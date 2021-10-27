#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 1:00:00
#SBATCH -J gatk_vcf_filter_1a
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters

## Paired or SE reads
reads=SE; # "paired" or "SE"

## General
MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
VCF_IN=$MAINDIR/data/radseq/radseq2020/genotypegvcfs1.$reads\
/bark_beetles_genotyped.vcf.gz;
OUTDIR0=$MAINDIR/working/Simon;
REF=$MAINDIR/data/reference/Itypographus_ref_customized\
/reference_customized.fasta;
METADATA=$MAINDIR/data/radseq/radseq2020/metadata;
functions=$MAINDIR/scripts/radseq2020/quality_control;

### Load modules
module load bioinfo-tools  GATK/4.1.4.1 vcftools/0.1.16 bcftools/1.10 \
R_packages/4.0.0 plink/1.90b4.9;

### Define functions
filter_stats=$functions/filter_stats.sh;
pca=$functions/pca.sh;
replicate_comp=$functions/replicate_comp.sh;
gatk_filters=$functions/gatk_filters.sh;

### Create folders
mkdir $OUTDIR0/gatk_vcf_filters_1.$reads;

### Apply filters
OUTDIR0=$OUTDIR0/gatk_vcf_filters_1.$reads;
DATA=1;
replicates=$METADATA/replicates_$DATA.txt;
sample_info=$METADATA/sample_info_$DATA.txt;
popmap=$METADATA/popmap_$DATA.txt;

## Filter 1: Extract SNPs
filter=1;
mkdir $OUTDIR0/filter$filter;
OUTDIR=$OUTDIR0/filter$filter;
VCF_OUT=$OUTDIR/bark_beetles_filter$filter;
gatk SelectVariants -R $REF -V $VCF_IN --select-type-to-include SNP \
-O $VCF_OUT.vcf.gz;
source $filter_stats;
source $pca;
source $replicate_comp;

## Filter 2: Minimum filters
VCF_IN=$VCF_OUT;
filter=2;
mkdir $OUTDIR0/filter$filter;
OUTDIR=$OUTDIR0/filter$filter;
VCF_OUT=$OUTDIR/bark_beetles_filter$filter;
vcftools --gzvcf $VCF_IN.vcf.gz --remove-indels --max-missing 0.5 \
--min-meanDP 10 --minDP 5 --minQ 30 --mac 3 \
--recode --recode-INFO-all --stdout | bgzip -c > $VCF_OUT.vcf.gz;
source $filter_stats;
source $pca;
source $replicate_comp;

## Filter 3: Remove repeat regions (soft masked lower case bases in the reference)
VCF_IN=$VCF_OUT;
filter=3;
mkdir $OUTDIR0/filter$filter;
OUTDIR=$OUTDIR0/filter$filter;
VCF_OUT=$OUTDIR/bark_beetles_filter$filter;
vcftools --gzvcf $VCF_IN.vcf.gz --mask $REF.repeats_masked \
--recode --recode-INFO-all --stdout | bgzip -c > $VCF_OUT.vcf.gz;
tabix $VCF_OUT.vcf.gz;
source $filter_stats;
source $pca;
source $replicate_comp;
source $gatk_filters;
