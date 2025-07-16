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
mkdir -p "CSV_Results/Cluster_Counts"
mkdir -p "CSV_Results/DEGs_All"
mkdir -p "CSV_Results/Marker_Genes_All"
mkdir -p "Raw_Data"
mkdir -p "Plots/Clustering_Plots/Marker_Feature_Plots/"
mkdir -p "Plots/Clustering_Plots/Marker_Violin_Plots/"
mkdir -p "Plots/DESEQ_Plots/Dispersion_Plots"
mkdir -p "Plots/DESEQ_Plots/Heatmaps"
mkdir -p "Plots/DESEQ_Plots/MA_Plots"
mkdir -p "Plots/DESEQ_Plots/PCAs"
mkdir -p "Plots/DESEQ_Plots/Volcano_Plots"
mkdir -p "Plots/Quality_Control"
mkdir -p "R_Data"
mkdir -p "Logs"

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

if [ ! -f "sc_experiment_config.R" ]; then
cat > sc_experiment_config.R <<'EOF'
# Experiment configuration file for Messing/Mayfield scRNA-seq analysis using Seurat
# Written by Geoffrey A. Dilly in June 2024

# Project_name - Datatype: String
# Name of the project for Seurat
scConfig.Project_name <- "GD_scRNA_Experiment"

# Project_name - Datatype: String
# Prefix for filenames of R objects
scConfig.Prefix <- "GD_scRNA"

# Home_folder - Datatype: String (Path)
# Home directory where all analysis directories will be located
# This is generally not needed but good to note
scConfig.Home_folder <-   "path/to/home/dir"

# Raw_data_folder - Datatype: String (Path)
# Location of the 10x cellranger runs that will be used as raw data
# Should contain data that can be read by the Read10X() function in Seurat
# Data should be located at <Raw_data_folder>/<Sample_name>/<Sample_name>/outs/filtered_feature_bc_matrix
# In order to run SoupX, .../unfiltered_feature_bc_matrix is also necessary
scConfig.Raw_data_folder <- "Raw_Data"

# mito_pattern - Datatype: String
# Str pattern that identifies mitochondrial genes
# Mouse = "^mt-"; Rat = "^Mt-"; Human = "^MT-"
scConfig.mito_pattern <- "^mt-"

# ribo_pattern - Datatype: String
# Str pattern that identifies ribosomal genes
scConfig.ribo_pattern <- "^Rp[ls]"

# nFeature_RNA_cutoff - Datatype: Int 
# Minimum number of nCount detected to include cell
# Recommended = ≥200
scConfig.nFeature_RNA_cutoff <- 200

# percent_mito_cutoff - Datatype: Int
# Maximum percentage of mitochondrial reads to include cell
# Recommended = 5
scConfig.percent_mito_cutoff <- 5

# percent_ribo_cutoff - Datatype: Int
# Maximum percentage of ribosomal reads to include cell
# Recommended = 5
scConfig.percent_ribo_cutoff <- 5

# expected_doublet_pct - Datatype: Int
# Expected percentage of doublets in the dataset
# Recommended = 7.5 (5 to 10%)
scConfig.expct_doublet_pct <- 7.5

# compute_soupx - Datatype: Bool (TRUE/FALSE)
# If TRUE, soupX will be run on the raw data
scConfig.compute_soupx <- FALSE

# soupx_adjust - Datatype: Bool (TRUE/FALSE)
# If TRUE, use SoupX to adjust the raw counts for ambient RNA contamination
scConfig.soupx_adjust <- FALSE

# remove_doublets - Datatype: Bool (TRUE/FALSE)
# If TRUE doublets identified by DoubletFinder be removed prior to SCT normalization
scConfig.remove_doublets <- FALSE

# remove_mito_genes - Datatype: Bool (TRUE/FALSE)
# If TRUE mitochondrial genes will be removed prior to clustering
scConfig.remove_mito_genes <- FALSE

# remove_ribo_genes - Datatype: Bool (TRUE/FALSE)
# If TRUE ribosomal genes will be removed prior to clustering
scConfig.remove_ribo_genes <- FALSE

# remove_top_nUMIs - Datatype: Bool (TRUE/FALSE)
# If TRUE the 25% of cells with the highest UMI count will be removed from the data
scConfig.remove_top_nUMIs <- FALSE

# clustering_PCAs - Datatype: Int
# Default is 25 but a better estimate can be obtained by running ElbowPlot()
scConfig.clustering_PCAs <- 25

# clustering_resolution - Datatype: Int
# Resolution may be estimated with Clustree
scConfig.clustering_resolution <- 0.5

# cluster_plot_ident - Datatype: String
# Identity used by cluster and DGE plots
scConfig.cluster_plot_ident <- "seurat_clusters"

# exogenous_genes - Datatype: List of Strings
# List of exogenous genes to be excluded from clustering
scConfig.exogenous_genes <- c()

EOF
fi

sed 's|"path/to/home/dir"|"'"$scRNA_home_dir"'"|' sc_experiment_config.R > temp && mv temp sc_experiment_config.R

echo -e "\nSetup Complete"