# 01_load_data.R
# Purpose: Load a list of 10x objects into Seurat and save them as RDS
# Author: Geoff Dilly

library(here)
library(Seurat)
library(foreach)
library(doParallel)
snRNA_home_dir <- here()
setwd(snRNA_home_dir)

# Load custom functions
source("R/modules/log_utils.R")

# Log the start time and a timestamped copy of the script
write(paste0("01_load_data - Start: ", Sys.time()), file = "snRNA_Log.txt", append = TRUE)
write_script_log("R/01_load_data.R")

# Load the configuration file and metadata
source("sc_experiment_config.R")
scConfig.Sample_metadata <- read.csv("sc_sample_metadata.csv")

# Check for required metadata columns
if (!("Sample_name" %in% colnames(scConfig.Sample_metadata))) {
  stop("Mandatory metadata column <Sample_name> is not present")
}

if (!("Treatment" %in% colnames(scConfig.Sample_metadata))) {
  stop("Mandatory metadata column <Treatment> is not present")
}

if (!("Sex" %in% colnames(scConfig.Sample_metadata))) {
  stop("Mandatory metadata column <Sex> is not present")
}

if (!("Raw_data_dir" %in% colnames(scConfig.Sample_metadata))) {
  stop("Mandatory metadata column <Raw_data_dir> is not present")
}

# Setup parallel backend
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores/2)
registerDoParallel(cl)

sample_list <- split(scConfig.Sample_metadata, seq_len(nrow(scConfig.Sample_metadata)))

foreach(sample = sample_list, .packages = c("Seurat")) %dopar% {
  sample_seurat_data <- Read10X(paste0(sample$Raw_data_dir, "/filtered_feature_bc_matrix"))
  sample_seurat <- CreateSeuratObject(counts = sample_seurat_data, project = scConfig.Project_name, min.cells = 1, min.features = 1)
  sample_seurat <- PercentageFeatureSet(sample_seurat, pattern = scConfig.mito_pattern, col.name = "percent_mito")
  sample_seurat <- PercentageFeatureSet(sample_seurat, pattern = scConfig.ribo_pattern, col.name = "percent_ribo")
  sample_seurat <- subset(sample_seurat, subset = nFeature_RNA > scConfig.nFeature_RNA_cutoff & percent_mito < scConfig.percent_mito_cutoff)
  for (col in setdiff(colnames(sample), "Raw_data_dir")) {
    sample_seurat[[col]] <- sample[[col]]
  }
  saveRDS(sample_seurat, file = paste0("R_Data/", sample$Sample_name, "_seurat.rds"))
  NULL
}

stopCluster(cl)
# Log the completion time
write(paste0("01_load_data - Finish: ", Sys.time()), file = "snRNA_Log.txt", append = TRUE)
