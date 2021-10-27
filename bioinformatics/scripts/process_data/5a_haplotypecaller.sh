#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 40:00:00
#SBATCH -J haplotypecaller1
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters

## Paired or SE reads
reads=SE; # "paired" or "SE"

## General
MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
WORKDIR=$MAINDIR/data/radseq/radseq2020/bamfiles/$reads\_bams;
OUTDIR=$MAINDIR/working/Simon;
REF=$MAINDIR/data/reference/Itypographus_ref_customized\
/reference_customized.fasta;
METADATA=$MAINDIR/data/radseq/radseq2020/metadata;

### Load modules
module load bioinfo-tools GATK/4.1.4.1 samtools/1.10 picard/2.23.4 \
gnuparallel/20180822;

### Create folders
mkdir $OUTDIR/haplotypecaller1.$reads;

### Index reference and create sequence dictionary
samtools faidx $REF;
java -jar /sw/apps/bioinfo/picard/2.23.4/rackham/picard.jar \
CreateSequenceDictionary -R $REF -O $(echo "${REF%.*}".dict);

### Call variants
OUTDIR=$OUTDIR/haplotypecaller1.$reads;
INDS=$METADATA/sample_info_1.txt;
mitochondrion=$(echo KY952782.1);
export REF WORKDIR OUTDIR mitochondrion reads;

## Nuclear genome
parallel --colsep '\t' --header : 'gatk --java-options "-Xmx4g" HaplotypeCaller \
-R $REF -XL $mitochondrion -I $WORKDIR/{1}.$reads.bam -O $OUTDIR/{1}.g.vcf.gz \
-contamination 0.10 -ERC GVCF --dont-use-soft-clipped-bases -mbq 20' :::: $INDS;

## Mitochondrial genome
parallel --colsep '\t' --header : 'gatk --java-options "-Xmx4g" HaplotypeCaller \
-R $REF -L $mitochondrion -I $WORKDIR/{1}.$reads.bam -O $OUTDIR/{1}.mitochondrion.g.vcf.gz \
-contamination 0.10 -ERC GVCF -ploidy 1 --dont-use-soft-clipped-bases -mbq 20' :::: $INDS;
