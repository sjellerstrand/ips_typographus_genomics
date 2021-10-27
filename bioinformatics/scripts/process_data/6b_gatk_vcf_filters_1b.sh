#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 2:00:00
#SBATCH -J gatk_vcf_filters_1b
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters

## Paired or SE reads
reads=SE; # "paired" or "SE"

## For GATK hard filters
QD=5.0; # Recommended=2
MQ=45.0; # Recommended=40
MQRankSum=-5.0; # Recommended=-12.5
ReadPosRankSum=-3.0; # Recommended=-8
ExcessHet=54.69; # 54.69? Equals z-score of 4.5

## General
MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
OUTDIR0=$MAINDIR/working/Simon/gatk_vcf_filters_1.$reads;
VCF_IN=$OUTDIR0/gatk_vcf_filters_1.SE/filter3/bark_beetles_filter3;
REF=$MAINDIR/data/reference/Itypographus_ref_customized\
/reference_customized.fasta;
METADATA=$MAINDIR/data/radseq/radseq2020/metadata;
functions=$MAINDIR/scripts/radseq2020/quality_control;
VCF_comb=$MAINDIR/data/radseq/radseq2020/genotypegvcfs1.$reads\
/bark_beetles_combined.vcf.gz;

### Load modules
module load bioinfo-tools GATK/4.1.4.1 vcftools/0.1.16 \
bcftools/1.10 R_packages/4.0.0 plink/1.90b4.9;

### Define functions
filter_stats=$functions/filter_stats.sh;
pca=$functions/pca.sh;
replicate_comp=$functions/replicate_comp.sh;
contam=$functions/contam.sh;

### Apply filters
DATA=1;
replicates=$METADATA/replicates_$DATA.txt;
sample_info=$METADATA/sample_info_$DATA.txt;
popmap=$METADATA/popmap_$DATA.txt;

## Filter 4: GATK Hard filters
filter=4;
mkdir $OUTDIR0/filter$filter;
OUTDIR=$OUTDIR0/filter$filter;
VCF_OUT=$OUTDIR/bark_beetles_filter$filter;
gatk VariantFiltration -R $REF -V  $VCF_IN.vcf.gz --filter-expression \
"QD < $QD || MQ < $MQ || MQRankSum < $MQRankSum \
|| ReadPosRankSum < $ReadPosRankSum || ExcessHet > $ExcessHet" \
--filter-name "GATK_hard_filters" -O $VCF_OUT.temp.vcf.gz;
gatk SelectVariants -R $REF -V $VCF_OUT.temp.vcf.gz \
--exclude-filtered true  -O $VCF_OUT.vcf.gz;
rm $VCF_OUT.temp.vcf.gz;
tabix $VCF_OUT.vcf.gz;
source $filter_stats;
source $pca;
source $replicate_comp;

## Filter 5: Allelic balance per heterozygous site
VCF_IN=$VCF_OUT;
filter=5;
mkdir $OUTDIR0/filter$filter;
OUTDIR=$OUTDIR0/filter$filter;
VCF_OUT=$OUTDIR/bark_beetles_filter$filter;
bcftools view -g het $VCF_IN.vcf.gz -O z | \
bcftools query -f '[%CHROM:%POS\t%DP\t%AD\n]' | \
awk '{gsub(",", "\t"); print}' | awk '$3 > 0 && $4 > 0' | \
awk '{sumAD[$1] += $2; sumalt[$1] += $4}; \
END{for (id in sumAD) {print id, "\t", sumAD[id], "\t", sumalt[id]/sumAD[id]}}' | \
awk '{gsub(":", "\t"); print}' > $VCF_OUT.allelic_balance;
cat $VCF_OUT.allelic_balance | \
awk '{ if($4 > 0.28 && $4 < 0.72 || $4 < 0.01 || $4 > 0.99) print}' \
> $VCF_OUT.balanced_loci;
vcftools --gzvcf $VCF_IN.vcf.gz --positions  $VCF_OUT.balanced_loci \
--recode --recode-INFO-all --stdout | bgzip -c > $VCF_OUT.vcf.gz;
source $filter_stats;
source $pca;
source $replicate_comp;

## Filter 6: Strand bias
VCF_IN=$VCF_OUT;
filter=6;
mkdir $OUTDIR0/filter$filter;
OUTDIR=$OUTDIR0/filter$filter;
VCF_OUT=$OUTDIR/bark_beetles_filter$filter;
bcftools view -g het $VCF_IN.vcf.gz -O z | \
bcftools query -f '[%CHROM\t%POS\n]' | uniq > $VCF_OUT.sites;
bcftools query $VCF_comb -R $VCF_OUT.sites -f '[%CHROM:%POS\t%SB\n]' | \
awk '{gsub(",", "\t"); print}' | awk '{SRF[$1] += $2; SRR[$1] += $3; \
SAF[$1] += $4; SAR[$1] += $5}; END{for (id in SRF) {print id, "\t", \
SRF[id], "\t", SRR[id], "\t", SAF[id], "\t", SAR[id]}}' | \
awk '{gsub(":", "\t"); print}' > $VCF_OUT.strandbias;
cat $VCF_OUT.strandbias | \
awk '{ if($5 / ($6 + 0.01) > 100 && $3 / ($4 + 0.01) > 100 || \
$6 / ($5 + 0.01) > 100 && $4 / ($3 + 0.01) > 100) print}' \
> $VCF_OUT.good_loci;
vcftools --gzvcf $VCF_IN.vcf.gz --positions $VCF_OUT.good_loci \
--recode --recode-INFO-all --stdout | bgzip -c > $VCF_OUT.vcf.gz;
source $filter_stats;
source $pca;
source $replicate_comp;
source $contam;
