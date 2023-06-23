#!/bin/bash
set -e
set -u
set -o pipefail

# conda activate netcis
# make sure cutadapt, bowtie2, and samtools can be found on PATH
cd /research/labs/immunology/rogerslm/m277102/projects/NetCIS || exit


fastq="/research/labs/immunology/rogerslm/data/2020_SBscreen"
# arg2="/research/labs/immunology/rogerslm/m277102/projects/IAS-revised/test_data/supp-table1.xlsx"
# arg3="/research/labs/immunology/rogerslm/m277102/projects/IAS-revised/test_data/Run5-Lane1-mapping_input.xlsx"
# arg4="/research/labs/immunology/rogerslm/m277102/projects/IAS-revised/test_data/Run5-Lane2-mapping_input.xlsx"
# arg5="/research/labs/immunology/rogerslm/m277102/projects/IAS-revised/test_data/Run5-Lane3-mapping_input.xlsx"
# arg6="2020_SB"
# python /research/labs/immunology/rogerslm/m277102/projects/IAS-revised/prepare_IAS_input.py $fastq $arg2 $arg3 $arg4 $arg5 $arg6

# ntask multiplied by npara should not exceed the number of available cores
ntasks=4
npara=12
output_prefix="2020_SB/all_files"
bowtieIndex="/research/labs/immunology/rogerslm/tools/bowtie2_indexes/GRCm39/GRCm39"
input="2020_SB/input.tsv"
chrom_bed="2020_SB/keep-chroms.bed"
irl_tpn="AAATTTGTGGAGTAGTTGAAAAACGAGTTTTAATGACTCCAACTTAAGTGTATGTAAACTTCCGACTTCAACTG"
irr_tpn="GGATTAAATGTCAGGAATTGTGAAAAAGTGAGTTTAAATGTATTTGGCTAAGGTGTATGTAAACTTCCGACTTCAACTG"
primer="GTAATACGACTCACTATAGGGCTCCGCTTAAGGGAC"
python netcis/preprocess_reads.py -d $fastq -o $output_prefix -b $bowtieIndex -i $input -l $irl_tpn -r $irr_tpn -p $primer -t $ntasks -n $npara -t $ntasks -n $npara -c $chrom_bed

njobs=$((ntasks * npara))
mapP_thresh=0.05
chrom_mapper="2020_SB/chrom_mapper.tsv"
python netcis/preprocess_insertions.py -o $output_prefix -i $input -c $chrom_mapper -t $mapP_thresh -j $njobs

edge_thresh=50000
python netcis/cis_networks.py -o $output_prefix -t $edge_thresh -j 4

# TODO:
# python netcis/analysis.py
