#!/bin/bash
set -e
set -u
set -o pipefail

# contains fasta files. ex.) chr1.fa,  chr12.fa, etc.
file_dir=/research/labs/immunology/rogerslm/m277102/projects/2023_SB/ref_data/seq_kit_ta_files
output_dir=/research/labs/immunology/rogerslm/m277102/projects/2023_SB/ref_data/ta_files

files=$(ls $file_dir)
for file in $files
do
    echo "$file"
    seqkit locate -P -p TA --bed < "$file_dir/$file" > "$output_dir/$file.bed"
done
