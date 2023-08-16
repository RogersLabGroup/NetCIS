
# NetCIS

Network based Common Insertional Site analysis

## Clone from github

[Github Repo](https://github.com/FischyM/NetCIS)

```bash
git clone git@github.com:FischyM/NetCIS.git
```

## Create python environment

```bash
conda create -n netcis python=3.10 numpy pandas scipy networkx seaborn docopt jupyterlab tqdm Biopython
conda activate netcis
pip install pysam
```

Jupyterlab is optional and can be installed as needed.

You will also need cutadapt, bowtie2, and samtools to be accesible from the command line.

## Prepare data

Data should be fastq files directly from sequencing. Preprocessing of the data happens in preprocess_reads.py where the following is done:

- remove adaptor tags from the reads (cutadapt)
- map reads to reference genome (bowtie2)
- process insertions
- save individual insertions per fastq file

## File Structure

TODO:

- input file for preprocess reads
- chrom mapper file for preprocess insertions
  - GRCm39 https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27/
  - https://genome.ucsc.edu/cgi-bin/hgTracks?chromInfoPage=&hgsid=1560703641_1YwiSDzyFEZ8nuDrTobTnwtYvReT
- gene annotation https://www.informatics.jax.org/genes.shtml 7/10/2023
  - https://www.informatics.jax.org/downloads/reports/MRK_List1.rpt
  - https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt

## Limitations

It is assummed that the SB screen data has IRL/IRR libraries and paired reads.

## Other Preprocessing Things

/project/cs-myers/MathewF/software/bowtie2-2.4.5/indexes

- download reference genome
- keep only major chromosomes
- concatenate individual chromosomes into one file
- create mapper index from ref genome (bowtie2)
- find all TA dimer sites in ref genome chroms into .bed format (seqkit)



# build index for bowtie2, hisat2, minimap2 for all the contigs and just the main contigs
# UCSC
# https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/ 
# mm39.fa.gz, mm39.chromFa.tar.gz
# https://www.biostars.org/p/113126/
# echo "$(ls chr*.fa | sort -V | grep -vP 'chr[^X|Y|\d]'; ls chr*.fa | sort -V | grep -vP 'chr[\d|X|Y]')" | xargs cat > GRCm39.fa
# cd /project/cs-myers/MathewF/software/bowtie2-2.4.5/indexes/GRCm39_major
# bowtie2-build chr10.fa,chr11.fa,chr12.fa,chr13.fa,chr14.fa,chr15.fa,chr16.fa,chr17.fa,chr18.fa,chr19.fa,chr1.fa,chr2.fa,chr3.fa,chr4.fa,chr5.fa,chr6.fa,chr7.fa,chr8.fa,chr9.fa,chrM.fa,chrX.fa,chrY.fa /project/cs-myers/MathewF/software/alignment_indexes/bowtie2/GRCm39-major
# hisat2-build chr10.fa,chr11.fa,chr12.fa,chr13.fa,chr14.fa,chr15.fa,chr16.fa,chr17.fa,chr18.fa,chr19.fa,chr1.fa,chr2.fa,chr3.fa,chr4.fa,chr5.fa,chr6.fa,chr7.fa,chr8.fa,chr9.fa,chrM.fa,chrX.fa,chrY.fa /project/cs-myers/MathewF/software/alignment_indexes/hisat2/GRCm39-major
# minimap2 -x sr -d /project/cs-myers/MathewF/software/alignment_indexes/minimap2/GRCm39-major.mmi /project/cs-myers/MathewF/software/bowtie2-2.4.5/indexes/GRCm39-major.fa.gz
# clear; minimap2 -ax sr /project/cs-myers/MathewF/software/bowtie2-2.4.5/indexes/GRCm39-major.fa.gz toy-data/2020_SB-fastq/trim2-1-1_RT_20170405001_20170405000_S465_L008_R1_001.fastq.gz toy-data/2020_SB-fastq/trim2-1-1_RT_20170405001_20170405000_S465_L008_R2_001.fastq.gz 
# clear; minimap2 -ax sr /project/cs-myers/MathewF/software/alignment_indexes/minimap2/GRCm39-major.mmi toy-data/2020_SB-fastq/trim2-1-1_RT_20170405001_20170405000_S465_L008_R1_001.fastq.gz toy-data/2020_SB-fastq/trim2-1-1_RT_20170405001_20170405000_S465_L008_R2_001.fastq.gz 
# 