#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 00:30:00
#SBATCH -J customize_reference
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
WORKDIR=$MAINDIR/data/reference;
OUTDIR=$MAINDIR/working/Simon;

### Load modules
module load bioinfo-tools vcftools/0.1.16;

### Create folders
mkdir $OUTDIR/Itypographus_ref_customized;

## Customise reference
OUTDIR=$OUTDIR/Itypographus_ref_customized;

## Concatenate nuclear and mitochondrial references
cat $WORKDIR/intermediate/Itypographus_ref_intermediate/SortedAndRenamed.fasta \
$WORKDIR/Itypographus_mitochondrion_ref/ips_mitochondrion.fasta \
> $OUTDIR/reference_customized.temp1.fasta;
cat $WORKDIR/intermediate/Itypographus_ref_repeats_masked/SortedAndRenamed.masked.fa \
$WORKDIR/Itypographus_mitochondrion_ref/ips_mitochondrion.fasta \
> $OUTDIR/reference_customized.temp1.fasta.repeats_masked;

## Remove contig 177 and 187 from reference genome
## (these contigs are redundant and were removed from the final assembly)
cat $OUTDIR/reference_customized.temp1.fasta | \
sed -n '1,/>IpsContig177/p;/IpsContig178/,$p' | \
sed -n '1,/>IpsContig187/p;/IpsContig188/,$p' | \
sed '/>IpsContig177/d' | sed '/>IpsContig187/d' \
> $OUTDIR/reference_customized.fasta;
cat $OUTDIR/reference_customized.temp1.fasta.repeats_masked | \
sed -n '1,/>IpsContig177/p;/IpsContig178/,$p' | \
sed -n '1,/>IpsContig187/p;/IpsContig188/,$p' | \
sed '/>IpsContig177/d' | sed '/>IpsContig187/d' \
> $OUTDIR/reference_customized.temp2.fasta.repeats_masked;

## Mask repeat regions
cat $OUTDIR/reference_customized.temp2.fasta.repeats_masked | \
sed -e '/^>/!s/[N]/1/g' | sed -e '/^>/!s/[A-Z]/0/g' \
> $OUTDIR/reference_customized.fasta.repeats_masked;

## Remove temporary files
rm $OUTDIR/reference_customized.temp*