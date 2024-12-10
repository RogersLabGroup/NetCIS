#!/bin/bash
set -e
set -u
set -o pipefail

# must run this before: conda activate netcis
conda env export > requirements.yml

# create environment in conda with 
# conda env create --name netcis --file requirements.yml
