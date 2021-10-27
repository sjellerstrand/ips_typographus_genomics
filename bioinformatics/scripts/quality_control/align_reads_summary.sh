### Check for errors in slurm and other reports
grep -iE "\b(err|e:|warn|w:|fail|abort)" <file>;

Quality filtering summary

quality_filtering_summary.sh

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

## Reads after demultiplexing
demulti() { IND=$1; echo $IND $(expr $(zcat $WORKDIR/demultiplex/$IND.1.fq.gz | wc -l) / 2);};
export -f demulti;
parallel 'demulti {}' ::: $INDS | sort -n -k1 > $OUTDIR/reads_demultiplexing.txt;

## Read pairs after PCR duplicate removal
dedup() { IND=$1; echo $IND $(expr $(zcat $WORKDIR/dedeuplicate/$IND.1.fq.gz | wc -l) / 2);};
export -f dedup;
parallel 'dedup {}' ::: $INDS | sort -n -k1 > $OUTDIR/reads_deduplicate.txt;

## Read pairs after quality trimming
trim() { IND=$1; echo $IND $(expr $(zcat $WORKDIR/quality_trim/$IND.* | wc -l) / 4);};
export -f trim;
parallel 'trim {}' ::: $INDS | sort -n -k1 > $OUTDIR/reads_trim.txt;

## SE read average length after quality trimming
se_length() { IND=$1; echo $IND $(zcat $WORKDIR/quality_trim/$IND.1.* | \
awk '{if(NR%4==2) {count++; bases += length} }END{print bases/count}');};
export -f se_length;
parallel 'se_length {}' ::: $INDS | sort -n -k1 > $OUTDIR/SE_length.txt;

## PE read average length after quality trimming
pe_length() { IND=$1; echo $IND $(zcat $WORKDIR/quality_trim/$IND.2.* | \
awk '{if(NR%4==2) {count++; bases += length} }END{print bases/count}');};
export -f pe_length;
parallel 'pe_length {}' ::: $INDS | sort -n -k1 > $OUTDIR/PE_length.txt;

Mapping summary

 align_reads_summary.sh

#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p core
#SBATCH -n 10
#SBATCH -t 3:00:00
#SBATCH -J align_reads_summary
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters

## Paired or SE reads
reads=SE; # "paired" or "SE"

## General
MAINDIR=/crex/proj/snic2020-6-222/Projects/Itypographus;
WORKDIR=$MAINDIR/data/radseq/radseq2020/bamfiles/$reads\_bams;
OUTDIR=$MAINDIR/working/Simon;
METADATA=$MAINDIR/data/radseq/radseq2020/metadata;

### Load modules
module load bioinfo-tools samtools/1.10 gnuparallel/20180822;

### Create folders
mkdir $OUTDIR/align_reads_summary;

### Evaluate performance by counting remaing read pairs at each step
OUTDIR=$OUTDIR/align_reads_summary;
INDS=$(cat $METADATA/sample_info_1.txt | tail -n +2 | cut -f 1);
export WORKDIR OUTDIR;

## Reads mapped
reads_mapped() { IND=$1; echo $IND $(samtools flagstat $WORKDIR/$IND.$reads.bam | \
awk -F "[+]" 'NR == 5 {print $1}');};
export -f reads_mapped;
parallel 'reads_mapped {}' ::: $INDS | sort -n -k1 > $OUTDIR/reads_mapped.txt;

## Total average read depth
tot_depth() { IND=$1; echo $IND $(samtools depth -a $WORKDIR/$IND.$readsbam | \
awk '{c++;s+=$3}END{print s/c}');};
export -f tot_depth;
parallel 'tot_depth {}' ::: $INDS | sort -n -k1 > $OUTDIR/total_depth.txt;

## Mapped average read depth
mapped_depth() { IND=$1; echo $IND $(samtools depth $WORKDIR/$IND.$reads.bam | \
awk '{c++;s+=$3}END{print s/c}');};
export -f mapped_depth;
parallel 'mapped_depth {}' ::: $INDS | sort -n -k1 > $OUTDIR/mapped_depth.txt;

## Breadth of coverage
coverage_breadth() { IND=$1; echo $IND $(samtools depth -a $WORKDIR/$IND.$reads.bam | \
awk '{c++; if($3>0) total+=1}END{print (total/c)*100}');};
export -f coverage_breadth;
parallel 'coverage_breadth {}' ::: $INDS | sort -n -k1 > $OUTDIR/coverage_breadth.txt;

## Merge files and add column headers
paste $OUTDIR/reads_mapped.txt $OUTDIR/total_depth.txt $OUTDIR/mapped_depth.txt \
$OUTDIR/coverage_breadth.txt -d ' ' | cut -f 1,2,4,6,8 -d ' ' > $OUTDIR/mapping_summary.temp;
echo -e "#Reads mapped\t#Total read depth\t#Mapped read depth\t#Breadth of coverage" | \
cat - $OUTDIR/mapping_summary.temp > $OUTDIR/mapping_summary.$reads.txt;

### Remove temporary files
rm $OUTDIR/reads_mapped.txt $OUTDIR/total_depth.txt $OUTDIR/mapped_depth.txt $OUTDIR/coverage_breadth.txt \
$OUTDIR/mapping_summary.temp;
