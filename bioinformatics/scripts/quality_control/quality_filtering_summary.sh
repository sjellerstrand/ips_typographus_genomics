#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p core
#SBATCH -n 8
#SBATCH -t 1:00:00
#SBATCH -J quality_filtering_summary
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
WORKDIR=$MAINDIR/data/radseq/radseq2020;
OUTDIR=$MAINDIR/working/Simon;
METADATA=$MAINDIR/data/radseq/radseq2020/metadata;

### Load modules
module load bioinfo-tools samtools/1.10 gnuparallel/20180822;

### Create folders
mkdir $OUTDIR/quality_filter_summary;

### Report quality filtering summary by counting remaing read pairs at each step
OUTDIR=$OUTDIR/quality_filter_summary;
INDS=$(cat $METADATA/sample_info_1.txt | tail -n +2 | cut -f 1);
export WORKDIR OUTDIR:

## Read pairs after demultiplexing
demulti() { IND=$1; echo $IND $(expr $(zcat $WORKDIR/demultiplex/$IND.1.fq.gz | wc -l) / 4);};
export -f demulti;
parallel 'demulti {}' ::: $INDS | sort -n -k1 > $OUTDIR/reads_demultiplexing.txt;

## Read pairs after PCR duplicate removal
dedup() { IND=$1; echo $IND $(expr $(zcat $WORKDIR/deduplicate/$IND.1.fq.gz | wc -l) / 4);};
export -f dedup;
parallel 'dedup {}' ::: $INDS | sort -n -k1 > $OUTDIR/reads_deduplicate.txt;

## Read pairs after quality trimming
trim() { IND=$1; echo $IND $(expr $(zcat $WORKDIR/quality_trim/$IND.1.paired.fq.gz | wc -l) / 4);};
export -f trim;
parallel 'trim {}' ::: $INDS | sort -n -k1 > $OUTDIR/reads_trim.txt;

## SE read average length after quality trimming
se_length() { IND=$1; echo $IND $(zcat $WORKDIR/quality_trim/$IND.1.paired.fq.gz | \
awk '{if(NR%4==2) {count++; bases += length} }END{print bases/count}');};
export -f se_length;
parallel 'se_length {}' ::: $INDS | sort -n -k1 > $OUTDIR/SE_length.txt;

## PE read average length after quality trimming
pe_length() { IND=$1; echo $IND $(zcat $WORKDIR/quality_trim/$IND.2.paired.fq.gz | \
awk '{if(NR%4==2) {count++; bases += length} }END{print bases/count}');};
export -f pe_length;
parallel 'pe_length {}' ::: $INDS | sort -n -k1 > $OUTDIR/PE_length.txt;