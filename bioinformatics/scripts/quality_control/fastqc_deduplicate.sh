#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 2:00:00
#SBATCH -J fastqc_deduplicate
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
WORKDIR=$MAINDIR/data/radseq/radseq2020/deduplicate;
OUTDIR=$MAINDIR/working/Simon;

### Load modules
module load bioinfo-tools FastQC/0.11.8 gnuparallel/20180822;

### Create folders
mkdir $OUTDIR/deduplicate;

### Concatenate reads
OUTDIR=$OUTDIR/deduplicate;
zcat $WORKDIR/*.1.fq.gz > $OUTDIR/deduplicate_SE.fq;
zcat $WORKDIR/*.2.fq.gz > $OUTDIR/deduplicate_PE.fq;

### Create quality report for deduplicated data
files=$(echo deduplicate_SE.fq deduplicate_PE.fq);

export WORKDIR OUTDIR;
parallel 'fastqc $OUTDIR/{} -o $OUTDIR' ::: $files;

### Remove temporary files
rm $OUTDIR/deduplicate_SE.fq $OUTDIR/deduplicate_PE.fq;
