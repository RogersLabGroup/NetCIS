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
# chrom_bed="toy-data/keep-chroms.bed"  # don't need chrom bed anymore, since I've personalized the mapper indexes
irl_tpn="AAATTTGTGGAGTAGTTGAAAAACGAGTTTTAATGACTCCAACTTAAGTGTATGTAAACTTCCGACTTCAACTG"
irr_tpn="GGATTAAATGTCAGGAATTGTGAAAAAGTGAGTTTAAATGTATTTGGCTAAGGTGTATGTAAACTTCCGACTTCAACTG"
primer="GTAATACGACTCACTATAGGGCTCCGCTTAAGGGAC"
python netcis/preprocess_reads.py -d $fastq -o $output_prefix -b $bowtieIndex -i $input -l $irl_tpn -r $irr_tpn -p $primer -t $ntasks -n $npara  # -c $chrom_bed

njobs=$((ntasks * npara))
mapP_thresh=0.05
# chrom_mapper="toy-data/chrom_mapper.tsv"
python netcis/preprocess_insertions.py -o $output_prefix -i $input -t $mapP_thresh -j $njobs  # -c $chrom_mapper

edge_thresh=50000
python netcis/cis_networks.py -o $output_prefix -t $edge_thresh -j $njobs
