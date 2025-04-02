# NetCIS: Network-based Common Insertion Site analysis

Forward genetic screens are powerful tools for biological discoveries as they do not require any a priori knowledge of genes controlling a phenotype of interest. Sleeping Beauty (SB) is a DNA transposon mutagenesis system that induces gain- and/or loss-of-gene-function when it inserts into or near a gene. SB is traditionally used to model heterogeneic and spontaneous malignancies and discover tumor driver genes. Mutations caused by SB insertions are typically observed using sequencing-based methods, then analyzed using computational tools to determine which insertions are driving the phenotype, termed Common Insertion Sites (CIS). Many computational analysis tools have been created to call CIS; however, none implement case-control comparisons in a direct and data-driven way. We have developed NetCIS: a Network-based Common Insertion Site analysis tool to robustly identify CIS in forward genetic screens with a case-control study design.

NetCIS uses mathematical graphs to identify pseudo-CIS (pCIS) in case and control samples. These pCIS are used to define genomic regions where pCIS from the case group can be compared to pCIS from the control group, reducing the problem space from the whole genome to observed insertion sites. Overlapping pCIS are tested for significance in their read counts per insertion sites with the Wilcoxon rank sums test and the Fisher’s exact test. We prioritize a data-driven approach and do not require that insertions are annotated to a reference TA site nor that a CIS is annotated to a gene or other genomic feature.

NetCIS was tested using data from a published SB screen studying intratumoral T cell trafficking(1). Briefly, T cells were mutagenized and allowed to infiltrate into tumors. Mutations in T cells that got into the tumor (case) were compared to mutations in T cells that were in the spleen but not in the tumor (control) to identify T cell mutations that support tumor infiltration. Three of these genes—Ehhadh, Sprr1b, and Aak1—were subsequently validated using in vivo adoptive transfer experiments(2). Using NetCIS, we successfully replicated identification of these experimentally validated genes. Importantly, we also discovered novel CIS, including miRNA Mir101c, Protein tyrosine phosphatase receptor type D (Ptprd), and several unannotated CIS that are missed by other methods.

1: Rogers LM, Wang Z, Mott SL, et al. A Genetic Screen to Identify Gain- and Loss-of-Function Modifications that Enhance T-cell Infiltration into Tumors. Cancer Immunol. Res. 2020; 8:1206–1214
2: Vianzon VV, Hanson RM, Garg I, et al. Rank aggregation of independent genetic screen results highlights new strategies for adoptive cellular transfer therapy of cancer. Front. Immunol. 2023; 14:1235131

NetCIS code and accompanying tutorial can be found at <https://github.com/RogersLabGroup/NetCIS>.  Tutorial data can be found at <https://zenodo.org/records/14733375>

## Requirements

Users should have a basic knowledge of Linux command line and will benefit with knowledge of Conda and Python.

Users should have access to a Linux x86_64 computer with a minimum of 8 cores and 32 GB RAM, and at least 300 GB of hard drive space. More cores will lead to faster run times, especially for step 1.

NetCIS has been tested on Ubuntu 22.04.5, with 32 cores of Intel(R) Xeon(R) CPU E5-2630 v3 @ 2.40GHz and 192 GB RAM. There is no guarantee for NetCIS to run on other operating systems.

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
git clone https://github.com/RogersLabGroup/NetCIS.git
```

Once downloaded, you can recreate the netcis environment with conda and pip.

```bash
conda env create --name netcis --file NetCIS/requirements.yml
conda activate netcis
pip install ranky
```

Download input files from Zenodo (zenodo_get has been installed with the netcis environment)

```bash
mkdir input_files
cd input_files
zenodo_get 14733375
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

Alternatively, you can run each step at a time, with the commands provided in Steps 1-4 below. All arguments are further detailed in their respective python files, which can be displayed by running the script with just the "-h" or "--help" argument.

```bash
python NetCIS/netcis/preprocess_reads.py --help
```

### Step 1: Preprocessing reads

- Ran with 16 parallel processes (-n) with 2 cores (-t) each for a total of 32 cores used
- Time to run: ~70 minutes
- Max memory used: ~20 GB

```bash
python NetCIS/netcis/preprocess_reads.py -d input_files/2020_SB/ -o output/results -b input_files/GRCm39 -i input_files/2020_SB.tsv -l AAATTTGTGGAGTAGTTGAAAAACGAGTTTTAATGACTCCAACTTAAGTGTATGTAAACTTCCGACTTCAACTG -r GGATTAAATGTCAGGAATTGTGAAAAAGTGAGTTTAAATGTATTTGGCTAAGGTGTATGTAAACTTCCGACTTCAACTG -p GTAATACGACTCACTATAGGGCTCCGCTTAAGGGAC -t 2 -n 16 -q 13 -v 1
```

```bash
python NetCIS/netcis/preprocess_insertions.py -o output/results -i input_files/2020_SB.tsv -j 32 -v 1
```

```bash
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

(For the first run) For gene set enrichment, redundant gene sets are removed. Removal uses scikit-learn's pairwise_distance function, which by default will use as many cores and memory as available. Once this is done, a copy of these pathways are saved in order to speed up subsequent runs of this step.

- Time to run: ~7.5 minutes
- Max memory used: ~33 GB, about 1 GB per core

(For subsequent runs)

- Time to run: ~2.5 minutes
- Max memory used: ~3 GB

```bash
python NetCIS/netcis/analysis.py -o output/results -a LT -b S -g input_files/MRK_List2.rpt -s input_files/m5.all.v2023.2.Mm.symbols.gmt -p 0.05 -x 5000 -m Gene -f "" -t 20000 -v 1
```

Results are separated by the combination of case (LT) and control (S) and edge threshold (20000) and in this example can be found in output/results-analysis/LT-S/20000.

Each subdirectory of results has the following structure:

- gene_set_enrichment_{case}/ (plots for case specific gene set enrichment)
- gene_set_enrichment_{control}/ (plots for control specific gene set enrichment)
- pyGT_helper/ (contains files to help build plots in pyGT_top_CIS)
- pyGT_top_CIS/ (genome track plots of the top CIS)
- volcano_plots/ (.png and .svg of the resulting volcano plots for rank sums and fisher's exact tests)
- CIS.tsv (Results for all common insertion sites with the following structure)
  - case_index: index for the case pCIS graph
  - case_pos_min: minimum genomic position for case pCIS graph
  - case_pos_max: maximum genomic position for case pCIS graph
  - control_index: index for the control pCIS graph
  - control_pos_min: minimum genomic position for control pCIS graph
  - control_pos_max: maximum genomic position for control pCIS graph
  - LFC: Log2 fold change of case/control read count (count per million)
  - ranksums: Rank sum statistical test, 2 sided
  - binomial: binomial statistical test, 2 sided
  - fishers_exact: fisher's exact statistical test, 2 sided
  - total_num_samples: total number of samples in CIS
  - case_num_samples: total number of case samples in CIS
  - control_num_samples: total number of control samples in CIS
  - total_IS: total number of insertion sites (unique genomic positions)
  - case_IS: total number of case insertion sites
  - control_IS: total number of control insertion sites
  - total_read_count: total read count as count per million
  - case_total_read_count: total read count from case samples
  - control_total_read_count: total read count from control samples
  - case: name of case group
  - control: name of control group
  - chrom: chromosome of CIS
  - genes: genes annotated to this CIS
  - ranksums-neglog: negative log of rank sums
  - fishers_exact-neglog: negative log of fisher's exact
  - binomial-neglog: negative log of binomial
  - ranksums-BY: rank sums corrected for multiple hypothesis testing using Benjaminini-Yekutieli method
  - fishers_exact-BY: fisher's exact corrected for multiple hypothesis testing using Benjaminini-Yekutieli method
  - binomial-BY: binomial corrected for multiple hypothesis testing using Benjaminini-Yekutieli method
  - ranksums-BY-neglog: negative log of corrected rank sums
  - fishers_exact-BY-neglog: negative log of corrected fisher's exact
  - binomial-BY-neglog: negative log of corrected binomial
  - rank: calculated Borda count rank using the median method
  - enriched: which group has more reads, case or control
  - CIS_start: genomic location of the beginning of the CIS
  - CIS_end: genomic location of the ending of the CIS
  - genome_viewer: string that is used to find the specific CIS in a genome viewer application
