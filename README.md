
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
  - columns are sample name. irl r1, irl r2, irr r1, irr r2
  - sample name should be formatted as <treatment>-<id> and only uses one dash ("-")
- chrom mapper file for preprocess insertions
  - GRCm39 <https://www.ncbi.nlm.nih.gov/assembly/GCF_000001635.27/>
  - <https://genome.ucsc.edu/cgi-bin/hgTracks?chromInfoPage=&hgsid=1560703641_1YwiSDzyFEZ8nuDrTobTnwtYvReT>
- gene annotation <https://www.informatics.jax.org/genes.shtml> 7/10/2023
  - <https://www.informatics.jax.org/downloads/reports/MRK_List1.rpt>
  - <https://www.informatics.jax.org/downloads/reports/MRK_List2.rpt>

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
