#!/bin/bash
set -e
set -u
set -o pipefail

# conda activate netcis
# make sure cutadapt, bowtie2, and samtools can be found on PATH
cd /research/labs/immunology/rogerslm/m277102/projects/NetCIS || exit


fastq="/research/labs/immunology/rogerslm/data/2020_SBscreen"
arg2="/research/labs/immunology/rogerslm/m277102/projects/IAS-revised/test_data/supp-table1.xlsx"
arg3="/research/labs/immunology/rogerslm/m277102/projects/IAS-revised/test_data/Run5-Lane1-mapping_input.xlsx"
arg4="/research/labs/immunology/rogerslm/m277102/projects/IAS-revised/test_data/Run5-Lane2-mapping_input.xlsx"
arg5="/research/labs/immunology/rogerslm/m277102/projects/IAS-revised/test_data/Run5-Lane3-mapping_input.xlsx"
arg6="2020_SB"
python /research/labs/immunology/rogerslm/m277102/projects/IAS-revised/prepare_IAS_input.py $fastq $arg2 $arg3 $arg4 $arg5 $arg6

# ntask multiplied by npara should not exceed the number of available cores
ntasks=1
npara=16
output_prefix="2020_SB/all_files"
bowtieIndex="/research/labs/immunology/rogerslm/tools/bowtie2_indexes/GRCm39/GRCm39"
input="2020_SB/input.tsv"
python netcis/preprocess_reads.py $fastq $output_prefix $bowtieIndex $ntasks $npara $input

njobs=$((ntasks * npara))
python netcis/preprocess_insertions.py $output_prefix $njobs $input

thres=50000
python netcis/cis_networks.py --output_prefix $output_prefix --threshold $thres --jobs $njobs

# TODO:
python netcis/analysis.py
