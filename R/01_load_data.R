# 01_load_data.R
# Purpose: Load a list of 10x objects into Seurat and save them as RDS
# Author: Geoff Dilly

library(here)
library(Seurat)
library(dplyr)
library(foreach)
library(doParallel)
library(SoupX)
scRNA_home_dir <- here()
setwd(scRNA_home_dir)

# For testing only!!!!
scConfig.compute_soupx <- TRUE

# Setup ####
# Load custom functions
source("R/modules/log_utils.R")
source("R/modules/soupx_utils.R")

# Log the start time and a timestamped copy of the script
write(paste0("01_load_data - Start: ", Sys.time()), file = "scRNA_Log.txt", append = TRUE)
write_script_log("R/01_load_data.R")

# Load the configuration file and metadata
source("sc_experiment_config.R")
scConfig.Sample_metadata <- read.csv("sc_sample_metadata.csv")

# Setup parallel backend
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores / 2)
registerDoParallel(cl)
on.exit(stopCluster(cl))

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

# Load data and run SoupX (optional) ####
sample_list <- split(scConfig.Sample_metadata, seq_len(nrow(scConfig.Sample_metadata)))

top_ambient_genes <- foreach(sample = sample_list, .packages = c("Seurat", "SoupX", "dplyr")) %dopar% {
  top_ambient <- NULL
  filt_matrix <- Read10X(file.path(sample$Raw_data_dir, "filtered_feature_bc_matrix"))
  sample_seurat <- CreateSeuratObject(counts = filt_matrix, project = scConfig.Project_name, min.cells = 1, min.features = 1)

  if (scConfig.compute_soupx == TRUE) {
    soupx_results <- run_soupx_correction(sample, sample_seurat)
    sample_seurat <- soupx_results$seurat_obj
    top_ambient <- soupx_results$top_ambient
  }

  # Filter the Seurat object for QC and add metadata
  sample_seurat <- PercentageFeatureSet(sample_seurat, pattern = scConfig.mito_pattern, col.name = "percent_mito")
  sample_seurat <- PercentageFeatureSet(sample_seurat, pattern = scConfig.ribo_pattern, col.name = "percent_ribo")
  sample_seurat <- subset(sample_seurat, subset = nFeature_RNA > scConfig.nFeature_RNA_cutoff &
                            percent_mito < scConfig.percent_mito_cutoff &
                            percent_ribo < scConfig.percent_ribo_cutoff)
  for (col in setdiff(colnames(sample), "Raw_data_dir")) {
    sample_seurat[[col]] <- sample[[col]]
  }
  saveRDS(sample_seurat, file = paste0("R_Data/", sample$Sample_name, "_seurat.rds"))

  return(top_ambient)
}

stopCluster(cl)

if (scConfig.compute_soupx == TRUE) {
  # Combine and write SoupX summary
  summary_df <- bind_rows(top_ambient_genes)
  write.csv(summary_df, "CSV_Results/Ambient_genes_summary.csv", row.names = FALSE)
}

# Log the completion time
write(paste0("01_load_data - Finish: ", Sys.time()), file = "scRNA_Log.txt", append = TRUE)
