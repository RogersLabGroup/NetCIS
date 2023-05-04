
# NetCIS

Network based Common Insertional Site analysis

## Clone from github

[Github Repo](https://github.com/FischyM/NetCIS)

```bash
git clone git@github.com:FischyM/NetCIS.git
```

## Create python environment

```bash
conda create -n netcis python=3.10
conda activate netcis
conda install numpy pandas scipy networkx seaborn docopt pysam
```

You will also need cutadapt, bowtie2, and samtools to be accesible from the command line.

## Prepare data

Data should be fastq files directly from sequencing. Preprocessing of the data happens in preprocess_reads.py where the following is done:

- remove adaptor tags from the reads (cutadapt)
- map reads to reference genome (bowtie2)
- process insertions
- save individual insertions per fastq file

## Limitations

It is assummed that the SB screen data has IRL/IRR libraries and paired reads.
