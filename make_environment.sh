#!/bin/bash

# Anaconda is required for this script
# Set up bioconda and conda-forge prior to running 

# Exit on error
set -e

# Check if conda is installed
if ! command -v conda &> /dev/null; then
    echo "Error: conda could not be found. Please install Anaconda."
    exit 1
fi

# Set the name of the conda environment
ENV_NAME="sc_analysis_env"

# Source conda in the current shell for 'conda activate'
source "$(conda info --base)/etc/profile.d/conda.sh"

# Creates a new Conda environment with the specified name
conda create -y -n $ENV_NAME -c conda-forge r-base=4.3.3 \
r-essentials r-devtools 

# Install all required conda-forge packages
conda install -y -n $ENV_NAME -c conda-forge r-seurat r-soupx r-doparallel 

# Install all required bioconductor packages
conda install -y -n $ENV_NAME -c bioconda \
bioconductor-dropletutils bioconductor-deseq2 bioconductor-glmgampoi \
r-pheatmap

# Activate the new environment
conda activate sc_analysis_env

# Install packages from GITHUB
Rscript -e 'devtools::install_github("chris-mcginnis-ucsf/DoubletFinder")'
Rscript -e 'devtools::install_github("immunogenomics/presto")'