#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t12:00:00
#SBATCH -J gatk_vcf_filters_1f
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters

## Paired or SE reads
reads=SE; # "paired" or "SE"

## General
MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
OUTDIR0=$MAINDIR/working/Simon/gatk_vcf_filters_1.$reads;
VCF_IN=$OUTDIR0/filter8/bark_beetles_filter8;
REF=$MAINDIR/data/reference/Itypographus_ref_customized\
/reference_customized.fasta;
METADATA=$MAINDIR/data/radseq/radseq2020/metadata;
functions=$MAINDIR/scripts/radseq2020/quality_control;

### Load modules
module load bioinfo-tools vcftools/0.1.16 vcflib/1.0.1 \
bcftools/1.10 R_packages/4.0.0 plink/1.90b4.9;

### Define functions
filter_stats=$functions/filter_stats.sh;
pca=$functions/pca.sh;
contam=$functions/contam.sh;

### Create folders
mkdir $MAINDIR/working/Simon/vcfs_final \
$MAINDIR/working/Simon/vcfs_final/$reads \
$MAINDIR/working/Simon/vcfs_final/$reads/non_pruned \
$MAINDIR/working/Simon/vcfs_final/$reads/pruned;

### Apply filters
DATA=all;
replicates=$METADATA/replicates_$DATA.txt;
sample_info=$METADATA/sample_info_$DATA.txt;
popmap=$METADATA/popmap_$DATA.txt;
order=$METADATA/order_2.txt;
remove=$METADATA/remove_inds.2b.txt;

## Filter 9: Remove duplicate samples and down sample dataset
filter=9;
mkdir $OUTDIR0/filter$filter;
OUTDIR=$OUTDIR0/filter$filter;
VCF_OUT=$OUTDIR/bark_beetles_filter$filter;
vcftools --gzvcf $VCF_IN.vcf.gz --remove $remove --max-missing 0.8 --mac 3 \
--recode --recode-INFO-all --stdout | bgzip -c > $VCF_OUT.vcf.gz;
tabix $VCF_OUT.vcf.gz;
source $filter_stats;
source $pca;

## Filter 10: Missingness per population
VCF_IN=$VCF_OUT;
filter=10;
mkdir $OUTDIR0/filter$filter;
OUTDIR=$OUTDIR0/filter$filter;
VCF_OUT=$OUTDIR/bark_beetles_filter$filter;
cat $popmap | cut -f 2 | sort | uniq > $VCF_OUT.populations;
POPS=$(cat $VCF_OUT.populations);
for pop in ${POPS[@]}; do
cat $sample_info | cut -f 1,3 | awk -v var="$pop" '$2~var{print}' > $VCF_OUT.$pop.keep;
vcftools --gzvcf $VCF_IN.vcf.gz --keep $VCF_OUT.$pop.keep --missing-site --out $VCF_OUT.$pop;
done;
cat $VCF_OUT*.lmiss | awk '!/CHR/' | awk '$6 > 0.4' | cut -f 1,2 | sort -k1,1 -k2,2 | uniq > $VCF_OUT.missing_loci;
rm $VCF_OUT.*.keep $VCF_OUT.*.lmiss $VCF_OUT.*.log;
vcftools --gzvcf $VCF_IN.vcf.gz --exclude-positions  $VCF_OUT.missing_loci \
--max-missing 0.9 --recode --recode-INFO-all --stdout | bgzip -c > $VCF_OUT.vcf.gz;
source $filter_stats;
source $pca;

## Create three datasets out of filtered vcf
VCF_IN=$VCF_OUT;
cat $sample_info | awk '$4 == "Transect"' | cut -f 1 >$METADATA/keep_transect.txt;
### Manually create a file with locations samples so that some transect samples can be included.
### Save the file to $METADATA/keep_locations.txt

## All samples
cp $VCF_IN.vcf.gz $MAINDIR/working/Simon/vcfs_final/$reads/non_pruned/bark_beetles_all.vcf.gz;
mkdir $OUTDIR0/all;
OUTDIR=$OUTDIR0/all;
VCF_OUT=$OUTDIR/bark_beetles_all;
plink --vcf $VCF_IN.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out $VCF_OUT;
cat $VCF_OUT.prune.in | awk '{gsub(":", "\t"); print}' > $VCF_OUT.keep_sites;
vcftools --gzvcf $VCF_IN.vcf.gz --positions $VCF_OUT.keep_sites \
--recode --recode-INFO-all --stdout | bgzip -c > $VCF_OUT.pruned.vcf.gz;
VCF_OUT=$VCF_OUT.pruned;
source $filter_stats;
mv $VCF_OUT.vcf.gz $MAINDIR/working/Simon/vcfs_final/$reads/pruned;

## Transect samples
mkdir $OUTDIR0/transect;
OUTDIR=$OUTDIR0/transect;
VCF_OUT=$OUTDIR/bark_beetles_transect;
vcftools --gzvcf $VCF_IN.vcf.gz --keep $METADATA/keep_transect.txt --mac 3 \
--recode --recode-INFO-all --stdout | bgzip -c > $VCF_OUT.vcf.gz;
source $filter_stats;
cp $VCF_OUT.vcf.gz $MAINDIR/working/Simon/vcfs_final/$reads/non_pruned;
plink --vcf $VCF_OUT.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out $VCF_OUT;
cat $VCF_OUT.prune.in | awk '{gsub(":", "\t"); print}' > $VCF_OUT.keep_sites;
vcftools --gzvcf $VCF_OUT.vcf.gz --positions $VCF_OUT.keep_sites \
--recode --recode-INFO-all --stdout | bgzip -c > $VCF_OUT.pruned.vcf.gz;
mv $VCF_OUT.pruned.vcf.gz $MAINDIR/working/Simon/vcfs_final/$reads/pruned;

## Location samples
mkdir $OUTDIR0/locations;
OUTDIR=$OUTDIR0/locations;
VCF_OUT=$OUTDIR/bark_beetles_locations;
vcftools --gzvcf $VCF_IN.vcf.gz --keep $METADATA/keep_locations.txt --mac 3 \
--recode --recode-INFO-all --stdout | bgzip -c > $VCF_OUT.vcf.gz;
source $filter_stats;
cp $VCF_OUT.vcf.gz $MAINDIR/working/Simon/vcfs_final/$reads/non_pruned;
plink --vcf $VCF_OUT.vcf.gz --double-id --allow-extra-chr --set-missing-var-ids @:# \
--indep-pairwise 50 10 0.1 --out $VCF_OUT;
cat $VCF_OUT.prune.in | awk '{gsub(":", "\t"); print}' > $VCF_OUT.keep_sites;
vcftools --gzvcf $VCF_OUT.vcf.gz --positions $VCF_OUT.keep_sites \
--recode --recode-INFO-all --stdout | bgzip -c > $VCF_OUT.pruned.vcf.gz;
mv $VCF_OUT.pruned.vcf.gz $MAINDIR/working/Simon/vcfs_final/$reads/pruned;

## Creates new info files
cp $order $METADATA/order_all.txt;
cat $sample_info | fgrep "$(cat $METADATA/keep_transect.txt)" \
> $METADATA/sample_info_temp.txt;
cat $sample_info | head -n 1 | cat - $METADATA/sample_info_temp.txt \
> $METADATA/sample_info_transect.txt;
rm $METADATA/sample_info_temp.txt;
cat $popmap | fgrep "$(cat $METADATA/keep_transect.txt)" > $METADATA/popmap_transect.txt;
patterns=$(cat $METADATA/sample_info_transect.txt | cut -f 3 | sort | uniq | tr '\n' '|');
patterns="${patterns::-1}";
cat $order | egrep -i $patterns > $METADATA/order_transect.txt;
cat $sample_info | fgrep "$(cat $METADATA/keep_locations.txt)" \
> $METADATA/sample_info_temp.txt;
cat $sample_info | head -n 1 | cat - $METADATA/sample_info_temp.txt \
> $METADATA/sample_info_locations.txt;
rm $METADATA/sample_info_temp.txt;
cat $popmap | fgrep "$(cat $METADATA/keep_locations.txt)" > $METADATA/popmap_locations.txt;
patterns=$(cat $METADATA/sample_info_locations.txt | cut -f 3 | sort | uniq | tr '\n' '|');
patterns="${patterns::-1}";
cat $order | egrep -i $patterns > $METADATA/order_locations.txt;
