#!/bin/bash
set -e
set -u
set -o pipefail

# conda activate netcis
# make sure cutadapt, bowtie2, and samtools can be found on PATH
cd /research/labs/immunology/rogerslm/m277102/projects/NetCIS || exit

# ntask multiplied by npara should not exceed the number of available cores
ntasks=4
npara=2
fastq="toy-data/2020_SB-fastq"
output_prefix="toy-data/2020_SB"
bowtieIndex="/research/labs/immunology/rogerslm/tools/bowtie2_indexes/GRCm39/GRCm39"
input="toy-data/input.tsv"
python netcis/preprocess_reads.py $fastq $output_prefix $bowtieIndex $ntasks $npara $input

# njobs=$((ntasks * npara))
# mapq_thresh=0.05
# python netcis/preprocess_insertions.py $output_prefix $njobs $input $mapq_thresh

# edge_thresh=50000
# python netcis/cis_networks.py --output_prefix $output_prefix --threshold $edge_thresh --jobs $njobs
