#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p core
#SBATCH -n 6
#SBATCH -t 40:00:00
#SBATCH -J genotypegvcfs1
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters

## Paired or SE reads
reads=SE; # "paired" or "SE"

## General
MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
WORKDIR=$MAINDIR/data/radseq/radseq2020/haplotypecaller1.$reads;
OUTDIR=$MAINDIR/working/Simon;
REF=$MAINDIR/data/reference/Itypographus_ref_customized\
/reference_customized.fasta;
METADATA=$MAINDIR/data/radseq/radseq2020/metadata;

### Load modules
module load bioinfo-tools GATK/4.1.4.1;

### Create folders
mkdir $OUTDIR/genotypegvcfs1.$reads;

### Call variants
OUTDIR=$OUTDIR/genotypegvcfs1.$reads;
INDS=$(cat $METADATA/sample_info_1.txt | tail -n +2 | cut -f 1);
combine1=$(for IND in ${INDS[@]}; do echo -V $WORKDIR/$IND.g.vcf.gz; done);
combine2=$(for IND in ${INDS[@]}; do echo -V $WORKDIR/$IND.mitochondrion.g.vcf.gz; done);

## Nuclear genome
gatk --java-options "-Xmx96g -Xms96g" CombineGVCFs -R $REF $combine1 \
-O $OUTDIR/bark_beetles_combined.vcf.gz;
gatk --java-options "-Xmx96g -Xms96g" GenotypeGVCFs -stand-call-conf 20 \
-RF ProperlyPairedReadFilter -R $REF -V $OUTDIR/bark_beetles_combined.vcf.gz \
-O $OUTDIR/bark_beetles_genotyped.vcf.gz;

## Mitochondrial genome
gatk --java-options "-Xmx96g -Xms96g" CombineGVCFs -R $REF $combine2 \
-O $OUTDIR/bark_beetles_combined.mitochondrion.vcf.gz;
gatk --java-options "-Xmx96g -Xms96g" GenotypeGVCFs -R $REF \
-V $OUTDIR/bark_beetles_combined.mitochondrion.vcf.gz \
-O $OUTDIR/bark_beetles_genotyped.mitochondrion.vcf.gz -stand-call-conf 20;
