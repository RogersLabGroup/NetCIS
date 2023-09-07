#!/bin/bash
set -e
set -u
set -o pipefail

# contains fasta files. ex.) chr1.fa,  chr12.fa, etc.
file_dir=$1
output_dir=$2

files=(chr1.fa chr2.fa chr3.fa chr4.fa chr5.fa chr6.fa chr7.fa chr8.fa chr9.fa chr10.fa chr11.fa chr12.fa chr13.fa chr14.fa chr15.fa chr16.fa chr17.fa chr18.fa chr19.fa chrX.fa chrY.fa chrM.fa)
for file in "${files[@]}"
do
    seqkit locate -P -p TA --bed < "$file_dir/$file" > "$output_dir/$file.bed"
done
