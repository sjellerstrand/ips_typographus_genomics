#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 1:30:00
#SBATCH -J fastqc_raw
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
WORKDIR=$MAINDIR/data/radseq/radseq2020/rawdata/bark_beetle_Jun_16;
OUTDIR=$MAINDIR/working/Simon;

### Load modules
module load bioinfo-tools FastQC/0.11.8 gnuparallel/20180822;

### Create folders
mkdir $OUTDIR/fastqc $OUTDIR/fastqc/raw16jun;

### Create quality report for raw data 16 June
OUTDIR=$OUTDIR/fastqc/raw16jun;
files=$(find $WORKDIR -name "*.fastq.gz");

export OUTDIR;
parallel 'fastqc -o $OUTDIR' ::: $files;
