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
chrom_bed="toy-data/keep-chroms.bed"
irl_tpn="AAATTTGTGGAGTAGTTGAAAAACGA"
irr_tpn="GGATTAAATGTCAGGAATTGTGAAAA"
adaptor="TACCCATACGACGTCCCAGA"
python netcis/preprocess_reads.py -d $fastq -o $output_prefix -b $bowtieIndex -i $input -l $irl_tpn -r $irr_tpn -a $adaptor -t $ntasks -p $npara -c $chrom_bed

njobs=$((ntasks * npara))
mapP_thresh=0.05
chrom_mapper="2020_SB/chrom_mapper.tsv"
python netcis/preprocess_insertions.py -o $output_prefix -i $input -c $chrom_mapper -t $mapP_thresh -j $njobs

edge_thresh=50000
python netcis/cis_networks.py -o $output_prefix -t $edge_thresh -j $njobs
