#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p node
#SBATCH -n 1
#SBATCH -t 10:00:00
#SBATCH -J align_reads
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
WORKDIR=$MAINDIR/data/radseq/radseq2020/quality_trim;
OUTDIR=$MAINDIR/working/Simon;
REF=$MAINDIR/data/reference/Itypographus_ref_customized\
/reference_customized.fasta;
METADATA=$MAINDIR/data/radseq/radseq2020/metadata;

### Load modules
module load bioinfo-tools bwa/0.7.17 samtools/1.10 gnuparallel/20180822;

### Create folders
mkdir $OUTDIR/bamfiles $OUTDIR/bamfiles/paired_bams \
$OUTDIR/bamfiles/SE_bams;

### Index reference
bwa index $REF;

### Align reads to reference, add @RG header, sort bam files,
### and clip overlapping pairs on the PE read
OUTDIR=$OUTDIR/bamfiles;
OUTDIR1=$OUTDIR/paired_bams;
OUTDIR2=$OUTDIR/SE_bams;
INDS=$METADATA/sample_info_1.txt;
export REF WORKDIR OUTDIR OUTDIR1 OUTDIR2;

## Define function
align_pairs() {

# Input parameters
IND=$1;
LIB=$2;

# Align paired reads
bwa mem -M $REF $WORKDIR/$IND.1.paired.fq.gz \
$WORKDIR/$IND.2.paired.fq.gz -R "@RG\tID:$IND\tSM:$IND\tLB:$LIB\tPL:illumina" | \
samtools view -q 1 -b | samtools sort -T $REF > $OUTDIR1/$IND.paired.bam;

# Create SE bam files
samtools view $OUTDIR1/$IND.paired.bam -T $REF -f 64 -b \
> $OUTDIR2/$IND.SE.bam;

# Index bam files
samtools index $OUTDIR1/$IND.paired.bam;
samtools index $OUTDIR2/$IND.SE.bam;
};

## Excecute function in parallell
export -f align_pairs;
parallel --colsep '\t' --header : 'align_pairs {1} {5}' :::: $INDS;

mv $OUTDIR $MAINDIR/data/radseq/radseq2020/;
sbatch freebayes1.SE.sh
sbatch haplotypecaller1.SE.sh 
