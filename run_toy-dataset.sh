#!/bin/bash
set -e
set -u
set -o pipefail

# conda activate netcis
# make sure cutadapt, bowtie2, and samtools can be found on PATH
cd /project/cs-myers/MathewF/projects/Laura-SB-Analysis/NetCIS || exit

ntasks=4
npara=2
fastq="toy-data/2020_SB-fastq"
output_prefix="toy-data/2020_SB"
bowtieIndex="/project/cs-myers/MathewF/software/bowtie2-2.4.5/indexes/GRCm39/GRCm39"
input="toy-data/input.tsv"
python netcis/preprocess_reads.py $fastq $output_prefix $bowtieIndex $ntasks $npara $input
python netcis/preprocess_insertions.py $output_prefix $npara $input
