# NetCIS: Network-based Common Insertion Site analysis

Forward genetic screens have become powerful tools for biological discoveries as they do not require any a priori knowledge of genes controlling a phenotype of interest. Sleeping Beauty (SB) has become a tool that can cause gain and loss of gene function with its ability to insert a DNA transposon at any TA dinucleotide site and has become a tool used to recapitulate hematopoietic and epithelial malignancies that mimic human tumors. Many tools have been created to identify which mutational insertions are driving these generated malignancies, deemed Common Insertion Sites (CIS), however, none currently implement case-control comparisons in a direct data-driven way. We have developed NetCIS: a Network-based Common Insertion Site analysis tool to robustly identify CIS in a case-control phenotype selection screen.

NetCIS uses mathematical graphs to identify pseudo-CIS (pCIS) in case and control samples. These pCIS are used to determine a CIS, where a case pCIS overlaps with a control pCIS, or in a situation where there is no overlap in a pCIS. Overlapping pCISes are tested for significance in their read counts per insertion sites with the Wilcoxon rank sums test as well as the number of sample cases and controls that are found in the pCIS using the Fisher exact test. We prioritize a data-driven approach and do not require that insertions are annotated to a reference TA site nor that a CIS is annotated to a gene or other genomic feature.

We tested NetCIS on a previously published SB screen that used mutagenized T Cells to identify CIS mutations that can improve intratumoral trafficking. We were able to recapitulate in our top ranked CIS the previously validated genes Ehhadh, Sprr1b, and Aak1, in addition to discovering novel annotated CIS to microRNA 101c and the protein tyrosine phosphatase receptor type D (Ptprd), along with many unannotated CIS that would have been missed by other methods.

NetCIS code and accompanying tutorial can be found at <https://github.com/FischyM/NetCIS>.

## Requirements

Users should have a basic knowledge of Linux command line and will benefit with knowledge of Conda and Python.

Users should have access to a Linux x86_64 computer with a minimum of 8 cores and 32 GB RAM, and at least 300 GB of hard drive space. More cores will lead to faster run times, especially for step 1.

NetCIS has been tested on Ubuntu 22.04.5, with 32 cores of Intel(R) Xeon(R) CPU E5-2630 v3 @ 2.40GHz and 192 GB RAM. There is no gurantee for NetCIS to run on other operating systems.

## Setup

If not already installed, download and install Conda following this documentation <https://docs.anaconda.com/miniconda/install/#quick-command-line-install>. Then run the following commands to add the following to the conda channels list.

```bash
conda config --add channels bioconda
conda config --add channels conda-forge
```

Make a directory to download NetCIS and related tutorial files to.  

```bash
mkdir netcis_tutorial
cd netcis_tutorial
```

NetCIS can be obtained with the following command.

```bash
git clone git@github.com:FischyM/NetCIS.git
```

Once downloaded, you can recreate the netcis environment with conda and pip.

```bash
conda env create --name netcis --file NetCIS/requirements.yml
conda activate netcis
pip install ranky
```

Download input files from Zenodo

```bash
mkdir input_files
cd input_files
TODO:
```

Extract archive (this might take some time). All extracted files will require 100 GB of space. This will temporarily use about 200 GB of space, but once extracted the archive can be deleted to save space.

```bash
tar -xf 2020_SB_fastq.tar.gz
rm 2020_SB_fastq.tar.gz
```

Move shell script to run the full NetCIS analysis

```bash
mv run_2020_SB.sh ..
cd ..
```

## Run NetCIS

Run the full analysis using the script provided.

```bash
./run_2020_SB.sh
```

Alternatively you can run each step at a time, however, the commands will be quite long, but are provided below. All arguments are detailed further in their respective python files which can be displayed by running the script with just the "-h" or "--help" argument.

```bash
python NetCIS/netcis/preprocess_reads.py --help
```

### Step 1: Preprocessing reads

- Ran with 16 parallel processes (-n) with 2 cores (-t) each for a total of 32 cores used
- Time to run: ~70 minutes
- Max memory used: ~20 GB

```bash
python NetCIS/netcis/preprocess_reads.py -d input_files/2020_SB/ -o output/results -b input_files/GRCm39 -i input_files/2020_SB.tsv -l AAATTTGTGGAGTAGTTGAAAAACGAGTTTTAATGACTCCAACTTAAGTGTATGTAAACTTCCGACTTCAACTG -r GGATTAAATGTCAGGAATTGTGAAAAAGTGAGTTTAAATGTATTTGGCTAAGGTGTATGTAAACTTCCGACTTCAACTG -p GTAATACGACTCACTATAGGGCTCCGCTTAAGGGAC -t 2 -n 16 -q 13 -v 1

python NetCIS/netcis/preprocess_insertions.py -o output/results -i input_files/2020_SB.tsv -j 32 -v 1

python NetCIS/netcis/insertions_to_bed.py -o output/results -v 1
```

### Step 2: Generating pCIS networks

- Time to run: <1 minute
- Max memory used: <2 GB

```bash
python NetCIS/netcis/pcis_networks.py -o output/results -t 20000 -j 32 -v 1
```

### Step 3: pCIS to CIS

- Time to run: <1 minute
- Max memory used: <2 GB

```bash
python NetCIS/netcis/cis_networks.py -o output/results -a LT -b S -j 32 -t 20000 -v 1
```

### Step 4: Analysis

For the first run. For gene set enrichment, we remove redundant gene sets. This uses scikit-learn's pairwise_distance function which by default will use as many cores and memory as available. Once this is done, a copy of these pathways are saved in order to speed up subsequent runs of this step.

- Time to run: ~7.5 minutes
- Max memory used: ~33 GB, about 1 GB per core

For subsequent runs

- Time to run: ~2.5 minutes
- Max memory used: ~3 GB

```bash
python NetCIS/netcis/analysis.py -o output/results -a LT -b S -g input_files/MRK_List2.rpt -s input_files/m5.all.v2023.2.Mm.symbols.gmt -p 0.05 -x 5000 -m Gene -f "" -t 20000 -v 1
```

Results are separated by the combination of case (LT) and control (S) and edge threshold (20000) and in this instance can be found in output/results-analysis/LT-S/20000
