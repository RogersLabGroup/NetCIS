#!/bin/bash
set -e
set -u
set -o pipefail

# must run this before: conda activate netcis
conda env export > requirements.yml

# create environment in conda with 
# conda env create --name netcis --file requirements.yml



# conda install
# python=3.11 numpy pandas seaborn biopython tqdm docopt pysam scipy 
# scikit-learn networkx adjusttext gseapy jupyterlab entrez-direct sra-tools

# pip install
# ranky
