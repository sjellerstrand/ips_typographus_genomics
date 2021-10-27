#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 4:00:00
#SBATCH -J fastqc_trim
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
WORKDIR=$MAINDIR/data/radseq/radseq2020/quality_trim;
OUTDIR=$MAINDIR/working/Simon;

### Load modules
module load bioinfo-tools FastQC/0.11.8 gnuparallel/20180822;

### Create folders
mkdir $OUTDIR/quality_trim;

### Concatenate reads
OUTDIR=$OUTDIR/quality_trim;
zcat $WORKDIR/*.1.paired.fq.gz $WORKDIR/*.1.unpaired.fq.gz > $OUTDIR/trim_SE.fq;
zcat $WORKDIR/*.2.paired.fq.gz $WORKDIR/*.2.unpaired.fq.gz > $OUTDIR/trim_PE.fq;

### Create quality report for quality trimmed reads
files=$(echo trim_SE.fq trim_PE.fq);

export WORKDIR OUTDIR;
parallel 'fastqc $WORKDIR/{} -o $OUTDIR' ::: $files;

### Remove temporary files
rm $OUTDIR/trim_SE.fq $OUTDIR/trim_PE.fq;
