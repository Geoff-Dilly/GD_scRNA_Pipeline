# 01b_soupx.R
# Purpose: Run SoupX on single sample Seurat objects
# Author: Geoff Dilly

library(here)
library(Seurat)
library(dplyr)
library(SoupX)
library(foreach)
library(doParallel)
snRNA_home_dir <- here()
setwd(snRNA_home_dir)

# Setup ####
# Load custom functions
source("R/modules/log_utils.R")

# Log the start time and a timestamped copy of the script
write(paste0("01b_soupx - Start: ", Sys.time()), file = "snRNA_Log.txt", append = TRUE)
write_script_log("R/01b_soupx.R")

# Load the configuration file and metadata
source("sc_experiment_config.R", local = TRUE)
scConfig.Sample_metadata <- read.csv("sc_sample_metadata.csv")

# Setup parallel backend
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores / 2)
registerDoParallel(cl)

# Data Preprocessing ####
# Place each sample in a list for further processing
str_sample_list <- scConfig.Sample_metadata$Sample_name

# Code adapted from https://cellgeni.github.io/notebooks/html/new-10kPBMC-SoupX.html

# Run SoupX ####
# Process each sample
top_ambient_genes <- foreach(sample_name = str_sample_list, .packages = c("Seurat", "SoupX")) %dopar% {

  # Read the sample Seurat object
  sample_seurat <- LoadSeuratRds(paste0("R_Data/", sample_name, "_seurat.rds"))

  # Read the raw and filtered matrices for SoupX
  filt_matrix <- Read10X(paste0(sample$Raw_data_dir, "filtered_feature_bc_matrix/"))
  raw_matrix <- Read10X(paste0(sample$Raw_data_dir, "raw_feature_bc_matrix/"))

  # Cluster the Seurat object
  sample_seurat <- SCTransform(sample_seurat, verbose = TRUE)
  sample_seurat <- RunPCA(sample_seurat, verbose = TRUE)
  sample_seurat <- RunUMAP(sample_seurat, dims = 1:30, verbose = TRUE)
  sample_seurat <- FindNeighbors(sample_seurat, dims = 1:30, verbose = TRUE)
  sample_seurat <- FindClusters(sample_seurat, verbose = TRUE)

  # Create the SoupX channel
  soup_channel <- SoupChannel(raw_matrix, filt_matrix)

  # Extract metadata and UMAP coordinates
  meta <- sample_seurat@meta.data
  umap <- sample_seurat@reductions$umap@cell.embeddings

  # Assign clustering to SoupX
  soup_channel <- setClusters(soup_channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup_channel <- setDR(soup_channel, umap)

  # Estimate and apply the SoupX correction
  soup_channel <- autoEstCont(soup_channel)
  adj_matrix <- adjustCounts(soup_channel, roundToInt = TRUE)

  # Add as new assay
  sample_seurat[["SoupX"]] <- CreateAssayObject(counts = adj_matrix)

  # Sort soupProfile by estimated contamination and take the top 10
  top_ambient <- head(soup_channel$soupProfile[order(soup_channel$soupProfile$est, decreasing = TRUE), ], n = 10)
  top_ambient$Sample <- sample_name

  # Reload the original Seurat object (before clustering/SCT/UMAP)
  orig_seurat <- readRDS(paste0("R_Data/", sample_name, "_seurat.rds"))

  # Add the SoupX assay
  orig_seurat[["SoupX"]] <- CreateAssayObject(counts = adj_matrix)

  # Save the Seurat object with SoupX results
  saveRDS(orig_seurat, file = paste0("R_Data/", sample_name, "_seurat.rds"))

  top_ambient
}

# Combine and write summary
summary_df <- bind_rows(top_ambient_genes)
write.csv(summary_df, "CSV_Results/DEGs_All/Ambient_genes_summary.csv", row.names = FALSE)

# Log the completion time
write(paste0("01b_soupx - Finish: ", Sys.time()), file = "snRNA_Log.txt", append = TRUE)
