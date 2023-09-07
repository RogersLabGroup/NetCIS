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
# chrom_bed="toy-data/keep-chroms.bed"  # don't need chrom bed anymore, since I've personalized the mapper indexes
irl_tpn="AAATTTGTGGAGTAGTTGAAAAACGAGTTTTAATGACTCCAACTTAAGTGTATGTAAACTTCCGACTTCAACTG"
irr_tpn="GGATTAAATGTCAGGAATTGTGAAAAAGTGAGTTTAAATGTATTTGGCTAAGGTGTATGTAAACTTCCGACTTCAACTG"
primer="GTAATACGACTCACTATAGGGCTCCGCTTAAGGGAC"
python netcis/preprocess_reads.py -d $fastq -o $output_prefix -b $bowtieIndex -i $input -l $irl_tpn -r $irr_tpn -p $primer -t $ntasks -n $npara -t $ntasks -n $npara -m 0

njobs=$((ntasks * npara))
mapP_thresh=0.05
# chrom_mapper="toy-data/chrom_mapper.tsv"
python netcis/preprocess_insertions.py -o $output_prefix -i $input -t $mapP_thresh -j $njobs  # -c $chrom_mapper

edge_thresh=50000
python netcis/cis_networks.py -o $output_prefix -t $edge_thresh -j 4

ta_dir="/project/cs-myers/MathewF/software/bowtie2-2.4.5/indexes/GRCm39_TAs/"
gene_annot="/project/cs-myers/MathewF/projects/Laura-SB-Analysis/NetCIS/toy-data/MRK_List2.rpt"
ta_error=5
verbose=1
python netcis/analysis.py  -o $output_prefix -b $ta_dir -g $gene_annot -t $ta_error -v $verbose
