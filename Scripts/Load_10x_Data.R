# Load_10x_Data.R
# Purpose: Load a list of 10x objects into Seurat and save them as RDS
# Author: Geoff Dilly

library(here)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(DropletUtils)
snRNA_home_dir <- here()
setwd(snRNA_home_dir)

# Log the start time and a timestamped copy of the script
write(paste0("Load_10x_Data - Start: ", Sys.time()),file = "snRNA_Log.txt", append = TRUE)
file.copy("Scripts/Load_10x_Data.R", paste0("Logs/Time_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", "Load_10x_Data.R"), overwrite = FALSE)

# Load the configuration file and metadata
source("sc_experiment_config.R")
scConfig.Sample_metadata <- read.csv("sc_sample_metadata.csv")

# Check for required metadata columns
if(!("Sample_name" %in% colnames(scConfig.Sample_metadata))) {
  stop("Mandatory metadata column <Sample_name> is not present")
}

if(!("Treatment" %in% colnames(scConfig.Sample_metadata))) {
  stop("Mandatory metadata column <Treatment> is not present")
}

if(!("Sex" %in% colnames(scConfig.Sample_metadata))) {
  stop("Mandatory metadata column <Sex> is not present")
}

# Loads the Cell Ranger output into Seurat
for (i in seq_len(nrow(scConfig.Sample_metadata))) {
  sample <- scConfig.Sample_metadata[i, ]
  sample_seurat.data <- Read10X(paste0(scConfig.Raw_data_folder, "/", sample$Sample_name, "/", sample$Sample_name, "/outs/filtered_feature_bc_matrix"))
  sample_seurat <- CreateSeuratObject(counts = sample_seurat.data, project = scConfig.Project_name, min.cells = 1, min.features = 1)
  sample_seurat  <- PercentageFeatureSet(sample_seurat, pattern = scConfig.mito_pattern, col.name = "percent_mito")
  sample_seurat  <- PercentageFeatureSet(sample_seurat, pattern = scConfig.ribo_pattern, col.name = "percent_ribo")
  sample_seurat <- subset(sample_seurat, subset = nFeature_RNA > scConfig.nFeature_RNA_cutoff & percent_mito < scConfig.percent_mito_cutoff)

  # Loads each column of the metadata into the seurat object
  for(col in colnames(sample)) {
    sample_seurat[[col]] <- sample[[col]]
  }

  # Save the Seurat object
  saveRDS(sample_seurat, file = paste0("R_Data/", sample$Sample_name, "_seurat.rds"))
}

# Log the completion time
write(paste0("Load_10x_Data - Finish: ", Sys.time()), file = "snRNA_Log.txt", append = TRUE)
