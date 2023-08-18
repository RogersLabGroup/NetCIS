#!/bin/bash
set -e
set -u
set -o pipefail

# contains fasta files. ex.) chr1.fa,  chr12.fa, etc.
file_dir=GRCm39_major/

files=$(ls $file_dir)
for file in $files
do
    echo "$file"
    seqkit locate -P -p TA --bed < "$file_dir/$file" > "GRCm39_TAs/$file.bed"
done
