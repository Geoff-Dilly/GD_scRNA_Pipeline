#!/bin/bash

# install.sh — Set up single-cell RNA-seq analysis pipeline directory structure
# Author: Geoff Dilly

set -e

# Get the folder containing install.sh script
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd $SCRIPT_DIR

scRNA_home_dir=$(pwd)
echo -e "\nSetting up scRNA analysis pipeline in: \n$scRNA_home_dir\n"

# Create the directories that the single-cell pipeline needs to run
mkdir -p "csv_results/cluster_counts"
mkdir -p "csv_results/degs_all"
mkdir -p "csv_results/marker_genes_all"
mkdir -p "raw_data"
mkdir -p "plots/clustering_plots/marker_feature_plots/"
mkdir -p "plots/clustering_plots/marker_violin_plots/"
mkdir -p "plots/deseq_plots/dispersion_plots"
mkdir -p "plots/deseq_plots/heatmaps"
mkdir -p "plots/deseq_plots/ma_plots"
mkdir -p "plots/deseq_plots/pca_plots"
mkdir -p "plots/deseq_plots/volcano_plots"
mkdir -p "plots/quality_control"
mkdir -p "r_data"
mkdir -p "logs"

# Add a here file to the home directory
touch .here

# Create an r script that logs the home directory
if [ ! -f "analysis_home_dir.R" ]; then
echo 'scRNA_home_dir <- "'${scRNA_home_dir}'"' >> "analysis_home_dir.R"
fi

# Set up config files if they don't already exist
if [ ! -f "sc_sample_metadata.csv" ]; then
cat > sc_sample_metadata.csv <<'EOF'
Sample_name,Sex,Treatment
<SAMPLE>,<M>,<CONTROL>
EOF
fi

if [ ! -f "sc_experiment_config.yaml" ]; then
cat > sc_experiment_config.yaml <<'EOF'
# Experiment configuration file for Messing/Mayfield scRNA-seq analysis using Seurat
# Written by Geoffrey A. Dilly in June 2024

# Project_name - Datatype: String
# Name of the project for Seurat
project_name: GD_scRNA_Experiment

# Project_name - Datatype: String
# Prefix for filenames of R objects
prefix: GD_scRNA

# Home_folder - Datatype: String (Path)
# Home directory where all analysis directories will be located
# This is generally not needed but good to note
home_folder: __HOME_DIR__

# Raw_data_folder - Datatype: String (Path)
# Location of the 10x cellranger runs that will be used as raw data
# Should contain data that can be read by the Read10X() function in Seurat
# The pipeline does not use this variable, it only used the folder in metadata
raw_data_folder: raw_data

# mito_pattern - Datatype: String
# Str pattern that identifies mitochondrial genes
# Mouse = '^mt-'; Rat = '^Mt-'; Human = '^MT-'
mito_pattern: '^mt-'

# ribo_pattern - Datatype: String
# Str pattern that identifies ribosomal genes
# For rat/mouse use '^Rp[ls][[:digit:]]|^Rn[[:digit:]]'
ribo_pattern: '^Rp[ls][[:digit:]]|^Rn[[:digit:]]'

# nFeature_RNA_cutoff - Datatype: Int 
# Minimum number of nCount detected to include cell
# Recommended = ≥200
nFeature_RNA_cutoff: 200

# percent_mito_cutoff - Datatype: Int
# Maximum percentage of mitochondrial reads to include cell
# Recommended = 5
percent_mito_cutoff: 5

# percent_ribo_cutoff - Datatype: Int
# Maximum percentage of ribosomal reads to include cell
# Recommended = 5
percent_ribo_cutoff: 5

# expected_doublet_pct - Datatype: Int
# Expected percentage of doublets in the dataset
# Recommended = 7.5 (5 to 10%)
expct_doublet_pct: 7.5

# compute_soupx - Datatype: Bool (true/false)
# If true, soupX will be run on the raw data
compute_soupx: false

# soupx_adjust - Datatype: Bool (true/false)
# If true, use SoupX to adjust the raw counts for ambient RNA contamination
soupx_adjust: false

# remove_doublets - Datatype: Bool (true/false)
# If true, doublets identified by DoubletFinder will be removed prior to SCT normalization
remove_doublets: false

# remove_mito_genes - Datatype: Bool (true/false)
# If true, mitochondrial genes will be removed prior to clustering
remove_mito_genes: false

# remove_ribo_genes - Datatype: Bool (true/false)
# If true, ribosomal genes will be removed prior to clustering
remove_ribo_genes: false

# remove_top_nUMIs - Datatype: Bool (true/false)
# If true, the 25% of cells with the highest UMI count will be removed from the data
remove_top_nUMIs: false

# clustering_PCAs - Datatype: Int
# Default is 25 but a better estimate can be obtained by running ElbowPlot()
clustering_PCAs: 25

# clustering_resolution - Datatype: Float
# Resolution may be estimated with Clustree
clustering_resolution: 0.5

# cluster_plot_ident - Datatype: String
# Identity used by cluster and DGE plots
cluster_plot_ident: seurat_clusters

# exogenous_genes - Datatype: List of Strings
# List of exogenous genes to be excluded from clustering
exogenous_genes:
  - eGFP
  - Cre

# marker_gene_reference - Datatype: String
# reference from marker_genes_db.csv to use for marker gene plots
marker_gene_reference: Dilly_et_al_2022

EOF
fi

sed 's|"path/to/home/dir"|"'"$scRNA_home_dir"'"|' sc_experiment_config.R > temp && mv temp sc_experiment_config.R

if [ ! -f ".gitignore" ]; then
cat > .gitignore <<'EOF'
# Ignore generated data and logs
r_data/
raw_data/
logs/
plots/
csv_results/
snRNA_Log.txt

# Ignore home folder declaration
analysis_home_dir.R

EOF
fi

echo -e "\nSetup Complete"