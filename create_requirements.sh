#!/bin/bash
set -e
set -u
set -o pipefail

# must run this before: conda activate netcis
conda list -e > requirements.txt

# create environment in conda with 
# conda create --name <env> --file requirements.txt
