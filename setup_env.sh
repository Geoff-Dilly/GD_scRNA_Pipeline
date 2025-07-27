#!/bin/bash

# setup_env.sh — Set up conda environment for single-cell RNA-seq analysis pipeline
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

# Create the conda environment from the YAML file
conda env create -f sc_analysis_env.yaml --name $ENV_NAME
conda activate $ENV_NAME

# Install packages from Github
Rscript -e 'devtools::install_github("chris-mcginnis-ucsf/DoubletFinder")'
Rscript -e 'devtools::install_github("immunogenomics/presto")'

# Install glmGamPoi in R
# glmGamPoi causes issues in MacOS otherwise
Rscript -e 'BiocManager::install("glmGamPoi")'