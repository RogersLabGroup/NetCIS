
# NetCIS

Network based Common Insertional Site analysis

## Requirements

- It is assummed that the SB screen has been prepared as per LM-PCR protocol [TODO: link] and the fastq files are have IRL/IRR libraries and are sequenced with Illumina's short-read paired-end technology.
- Compute environment
  - NetCIS has only been tested on Linux TODO:
- Computer with enough storage space for fastq files and generated intermediate and result files

## Clone from github

[Github Repo](https://github.com/FischyM/NetCIS)

```bash
git clone git@github.com:FischyM/NetCIS.git
```

## Create python environment

```bash
conda create -n netcis -c bioconda -c conda-forge python=3.11 numpy pandas scipy networkx seaborn docopt jupyterlab tqdm Biopython gseapy pygenometracks pysam cutadapt, bowtie2, and samtools
conda activate netcis
```

Jupyterlab is optional and can be installed as needed.

## Prepare data

Data should be fastq files directly from sequencing. Preprocessing of the data happens in preprocess_reads.py where the following is done:

- remove adaptor tags from the reads (cutadapt)
- map reads to reference genome (bowtie2)
- process insertions
- save individual insertions per fastq file

## Input files

- fastq files
  - can be obtained from TODO:
- input file of fastq files
  - the following columns are sample_id, IRL_read1, IRL_read2, IRR_read1, IRR_read2, treatment, optional_columns, etc
    - first 6 columns are required
    - sample_id MUST be unique and does not need to be meaningful
    - optional_columns allows you to add any extra columns and they will propagate throughout NetCIS. For example, we have added tumor_model and pd1_treated
- reference mapping file for bowtie2
  - we have created our own GRCm39 reference mapping file with conventional contig names (chr1, chr2, etc.) and it can be obtained from our Zenodo repository TODO:
  - GRCm39 <https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27/>
  - <https://genome.ucsc.edu/cgi-bin/hgTracks?chromInfoPage=&hgsid=1560703641_1YwiSDzyFEZ8nuDrTobTnwtYvReT>
- gene annotation <https://www.informatics.jax.org/genes.shtml> 7/10/2023
  - <https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt>

## Limitations



## Other Preprocessing Things

/project/cs-myers/MathewF/software/bowtie2-2.4.5/indexes

- download reference genome
- keep only major chromosomes
- concatenate individual chromosomes into one file
- create mapper index from ref genome (bowtie2)
- find all TA dimer sites in ref genome chroms into .bed format (seqkit)

# build index for bowtie2, hisat2, minimap2 for all the contigs and just the main contigs

UCSC - <https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/>

```bash
mm39.fa.gz, mm39.chromFa.tar.gz
```

How to process keep only chromosomes 1-22,X,Y - <https://www.biostars.org/p/113126/>

```bash
echo "$(ls chr*.fa | sort -V | grep -vP 'chr[^X|Y|\d]'; ls chr*.fa | sort -V | grep -vP 'chr[\d|X|Y]')" | xargs cat > GRCm39.fa

cd /project/cs-myers/MathewF/software/bowtie2-2.4.5/indexes/GRCm39_major
```

# Disclosures

ChatGPT 3.5 used to assist in generating per-function documentation to be used with Sphinx to generate full html documentation.

Input: Provide documentation for this code similar to how scikit-learn does - (python file for documentation)

Output: Sure, here's a structured documentation for the provided code, similar to how scikit-learn documents its functions and classes:

Input: I need documentation for every function

Output: Certainly, here's the documentation for each function in your code:
