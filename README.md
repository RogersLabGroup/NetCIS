# NetCIS: Network-based Common Insertion Site analysis

Forward genetic screens have become powerful tools for biological discoveries as they do not require any a priori knowledge of genes controlling a phenotype of interest. Sleeping Beauty (SB) has become a tool that can cause gain and loss of gene function with its ability to insert a DNA transposon at any TA dinucleotide site and has become a tool used to recapitulate hematopoietic and epithelial malignancies that mimic human tumors. Many tools have been created to identify which mutational insertions are driving these generated malignancies, deemed Common Insertion Sites (CIS), however, none currently implement case-control comparisons in a direct data-driven way. We have developed NetCIS: a Network-based Common Insertion Site analysis tool to robustly identify CIS in a case-control phenotype selection screen.

NetCIS uses mathematical graphs to identify pseudo-CIS (pCIS) in case and control samples. These pCIS are used to determine a CIS, where a case pCIS overlaps with a control pCIS, or in a situation where there is no overlap in a pCIS. Overlapping pCISes are tested for significance in their read counts per insertion sites with the Wilcoxon rank sums test as well as the number of sample cases and controls that are found in the pCIS using the Fisher exact test. We prioritize a data-driven approach and do not require that insertions are annotated to a reference TA site nor that a CIS is annotated to a gene or other genomic feature. 

We tested NetCIS on a previously published SB screen that used mutagenized T Cells to identify CIS mutations that can improve intratumoral trafficking. We were able to recapitulate in our top ranked CIS the previously validated genes Ehhadh, Sprr1b, and Aak1, in addition to discovering novel annotated CIS to microRNA 101c and the protein tyrosine phosphatase receptor type D (Ptprd), along with many unannotated CIS that would have been missed by other methods.

NetCIS code and accompanying tutorial can be found at https://github.com/FischyM/NetCIS.

## Requirements

A Linux x86_64 computer with a minimum of 8 cores and 32 GB RAM (we recommended 32 cores and 196 GB RAM)

Users should have a basic knowledge of Linux command line and will benefit with knowledge of Conda and Python.

NetCIS has been tested on Ubuntu 22.04.1, with 32 cores and 192 GB RAM. It is possible to run NetCIS with less computational resources, however, more resources (i.e., more cores) will lead to faster run times. There is no gurantee for NetCIS to run on other operating systems.

## Setup

If not already installed, download and install Conda following this documentation <https://docs.anaconda.com/miniconda/install/#quick-command-line-install>

Then, make a directory to download NetCIS and related tutorial files to. NetCIS can be obtained with the following command. Once downloaded, you can recreate the netcis environment with conda.

```bash
mkdir netcis_tutorial
cd netcis_tutorial
git clone git@github.com:FischyM/NetCIS.git
conda env create --name netcis --file NetCIS/requirements.yml
conda activate netcis
```

Download input files from Zenodo

```bash

```

Download fastq files from NCBI BioProject. The conda environment already has NCBI's command line tools installed for this next step.

```bash

```
