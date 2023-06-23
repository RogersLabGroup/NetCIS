
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

## Limitations

It is assummed that the SB screen data has IRL/IRR libraries and paired reads.

## Other Preprocessing Things

/project/cs-myers/MathewF/software/bowtie2-2.4.5/indexes

- download reference genome
- keep only major chromosomes
- concatenate individual chromosomes into one file
- create mapper index from ref genome (bowtie2)
- find all TA dimer sites in ref genome chroms into .bed format (seqkit)
