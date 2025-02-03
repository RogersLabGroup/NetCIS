#!/bin/bash
set -e
set -u
set -o pipefail

# must run this before: conda activate netcis
conda env export > requirements.yml

# create environment in conda with 
# conda env create --name netcis --file requirements.yml


# conda config --add channels bioconda
# conda config --add channels conda-forge
# conda create -n netcis python=3.11 numpy pandas seaborn biopython tqdm docopt pysam scipy scikit-learn networkx adjusttext pygenometracks gseapy jupyterlab samtools bowtie2=2.5.4 cutadapt=5.0 
# pip install ranky
