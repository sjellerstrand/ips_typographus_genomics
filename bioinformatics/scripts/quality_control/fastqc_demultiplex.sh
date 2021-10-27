#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p core
#SBATCH -n 2
#SBATCH -t 5:00:00
#SBATCH -J fastqc_demultiplex
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
WORKDIR=$MAINDIR/data/radseq/radseq2020/demultiplex;
OUTDIR=$MAINDIR/working/Simon;

### Load modules
module load bioinfo-tools FastQC/0.11.8 gnuparallel/20180822;

### Create folders
mkdir $OUTDIR/demultiplex;

### Concatenate reads
OUTDIR=$OUTDIR/demultiplex;
zcat $WORKDIR/*.1.fq.gz > $OUTDIR/demultiplex_SE.fq;
zcat $WORKDIR/*.2.fq.gz > $OUTDIR/demultiplex_PE.fq;

### Create quality report for demultiplexed data in parallel
files=$(echo demultiplex_SE.fq demultiplex_PE.fq);

export MAINDIR OUTDIR;
parallel 'fastqc $OUTDIR/{} -o $OUTDIR' ::: $files;

### Remove temporary files
rm $OUTDIR/demultiplex_SE.fq $OUTDIR/demultiplex_PE.fq;
