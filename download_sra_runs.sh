#!/bin/bash
set -e
set -u
set -o pipefail

# get SRA info 
# esearch -db sra -query PRJNA641272 | efetch -format runinfo | cut -d ',' -f 1 > sra_runs.txt

# download SRA runs in parallel
sra_runs=$1
jobs=$2
parallel -j "$jobs" "prefetch {} && fasterq-dump {}" ::: "$sra_runs"
