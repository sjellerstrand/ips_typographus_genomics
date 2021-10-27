#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 2:00:00
#SBATCH -J quality_trimming
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
WORKDIR=$MAINDIR/data/radseq/radseq2020/deduplicate;
OUTDIR=$MAINDIR/working/Simon;
METADATA=$MAINDIR/data/radseq/radseq2020/metadata;

### Load modules
module load bioinfo-tools trimmomatic/0.36 gnuparallel/20180822;

### Create folders
mkdir $OUTDIR/quality_trim;

### Remove adapters and quality trim reads
OUTDIR=$OUTDIR/quality_trim;
ADAPTERDIR=$METADATA/adapters;
INDS=$METADATA/sample_info_1.txt;

## Trim reads
export WORKDIR OUTDIR ADAPTERDIR;
parallel --colsep '\t' --header : 'trimmomatic PE $WORKDIR/{1}.1.fq.gz $WORKDIR/{1}.2.fq.gz \
$OUTDIR/{1}.1.paired.fq.gz $OUTDIR/{1}.1.unpaired.fq.gz $OUTDIR/{1}.2.paired.fq.gz \
$OUTDIR/{1}.2.unpaired.fq.gz ILLUMINACLIP:$ADAPTERDIR/adapters.{6}.fa:3:30:10:1:TRUE \
TRAILING:20 SLIDINGWINDOW:4:20 MINLEN:120' :::: $INDS;
