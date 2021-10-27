#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 1:30:00
#SBATCH -J remove_PCR_duplicates
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
WORKDIR=$MAINDIR/data/radseq/radseq2020/demultiplex;
OUTDIR=$MAINDIR/working/Simon;
METADATA=$MAINDIR/data/radseq/radseq2020/metadata;

### Load modules
module load bioinfo-tools Stacks/2.53 gnuparallel/20180822;

### Create folders
mkdir $OUTDIR/deduplicate;

### Remove PCR duplicates
OUTDIR=$MAINDIR/working/deduplicate;
INDS=$(cat $METADATA/sample_info_1.txt | tail -n +2 | cut -f 1);

export WORKDIR OUTDIR;
parallel 'clone_filter -1 $WORKDIR/{}.1.fq.gz -2 $WORKDIR/{}.2.fq.gz -o $OUTDIR -i gzfastq' ::: $INDS;

### Rename files
for IND in ${INDS[@]}; do
mv $OUTDIR/$IND.1.1.fq.gz $OUTDIR/$IND.1.fq.gz;
mv $OUTDIR/$IND.2.2.fq.gz $OUTDIR/$IND.2.fq.gz;
done;
