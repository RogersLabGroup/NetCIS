#!/bin/bash
set -e
set -u
set -o pipefail

# conda activate netcis
cd /project/cs-myers/MathewF/projects/Laura-SB-Analysis || exit


# ntask multiplied by npara should not exceed the number of available cores
ntasks=4
npara=7

mapping=GRCm39
name="2020_SB"
output_prefix=publication/output/$name/results
verbose=1

# Step 1: Preprocessing reads
fastq=publication/fastq/$name
input=publication/input/$name.tsv
mapper_index=/project/cs-myers/MathewF/software/alignment_indexes/bowtie2/$mapping  # add in README how to generate this
irl_tpn=AAATTTGTGGAGTAGTTGAAAAACGAGTTTTAATGACTCCAACTTAAGTGTATGTAAACTTCCGACTTCAACTG
irr_tpn=GGATTAAATGTCAGGAATTGTGAAAAAGTGAGTTTAAATGTATTTGGCTAAGGTGTATGTAAACTTCCGACTTCAACTG
primer=GTAATACGACTCACTATAGGGCTCCGCTTAAGGGAC
# python NetCIS/netcis/preprocess_reads.py -d $fastq -o $output_prefix -b $mapper_index -i $input -l $irl_tpn -r $irr_tpn -p $primer -t $ntasks -n $npara -m $mapQ_threshold


# Step 2: Preprocessing insertions
njobs=$((ntasks * npara))
mapQ_threshold=13  # default
# python NetCIS/netcis/preprocess_insertions.py -o $output_prefix -i $input -j $njobs -m $mapQ_threshold
# python NetCIS/netcis/insertions_to_bed.py -o $output_prefix --treatment RT
# python NetCIS/netcis/insertions_to_bed.py -o $output_prefix --treatment LT
# python NetCIS/netcis/insertions_to_bed.py -o $output_prefix --treatment S


# Step 3: Generating pCIS networks
edge_threshold=5000
# python NetCIS/netcis/pcis_networks.py -o $output_prefix -t $edge_threshold -j $njobs


# Step 4: pCIS -> CIS
python NetCIS/netcis/cis_networks.py -o $output_prefix -a LT -b S -j $njobs -v $verbose
python NetCIS/netcis/cis_networks.py -o $output_prefix -a RT -b S -j $njobs -v $verbose
python NetCIS/netcis/cis_networks.py -o $output_prefix -a LT -b RT -j $njobs -v $verbose


# Step 5: Analysis
gene_annot=/project/cs-myers/MathewF/projects/Laura-SB-Analysis/NetCIS/toy-data/MRK_List2.rpt
pval_threshold=0.05
marker_expander=5000
marker_type=Gene  # Gene, Other Genome Feature, QTL, Pseudogene, Cytogenetic Marker, DNA Segment, Transgene, Complex/Cluster/Region
feature_type=""  # Gene: protein coding gene, many other possibilites for 'Gene' or other marker types

# python NetCIS/netcis/analysis.py -o $output_prefix -a LT -b S -g $gene_annot -p $pval_threshold -x $marker_expander -m $marker_type -f "$feature_type" -v $verbose
# python NetCIS/netcis/analysis.py -o $output_prefix -a RT -b S -g $gene_annot -p $pval_threshold -x $marker_expander -m $marker_type -f "$feature_type" -v $verbose
# python NetCIS/netcis/analysis.py -o $output_prefix -a LT -b RT -g $gene_annot -p $pval_threshold -x $marker_expander -m $marker_type -f "$feature_type" -v $verbose
