#!/bin/bash -l

#SBATCH -A snic2020-5-582
#SBATCH -p core
#SBATCH -n 1
#SBATCH -t 5:00:00
#SBATCH -J demultiplex
#SBATCH --mail-user=simon.jacobsen_ellerstrand@biol.lu.se
#SBATCH --mail-type=FAIL

### Set parameters
INPUTDIR=$1;
barcodes=$2;
OUTDIR=$3;
report=$4;

### Load modules
module load bioinfo-tools Stacks/2.53;

### Demultiplex data
process_radtags -P -p $INPUTDIR -i gzfastq -b $barcodes -o $OUTDIR -y gzfastq \
--barcode_dist_1 1 -e sbfI -r -c -q --inline_null --adapter_1 AGATCGGAAGAGCACACGTCTGAACTCCAGTCA \
--adapter_2 AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT --adapter_mm 2 &> $report;
