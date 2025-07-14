# 01b_soupx.R
# Purpose: Run SoupX on single sample Seurat objects
# Author: Geoff Dilly

library(here)
library(Seurat)
library(dplyr)
library(SoupX)
library(foreach)
library(doParallel)
scRNA_home_dir <- here()
setwd(scRNA_home_dir)

# Setup ####
# Load custom functions
source("R/modules/log_utils.R")

# Log the start time and a timestamped copy of the script
write(paste0("01b_soupx - Start: ", Sys.time()), file = "scRNA_Log.txt", append = TRUE)
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
top_ambient_genes <- foreach(sample = sample_list, .packages = c("Seurat", "SoupX", "dplyr")) %dopar% {

  sample_name <- sample$Sample_name
  raw_dir     <- sample$Raw_data_dir

  # Read the raw and filtered matrices for SoupX
  filt_matrix <- Read10X(file.path(raw_dir, "filtered_feature_bc_matrix"))
  raw_matrix  <- Read10X(file.path(raw_dir, "raw_feature_bc_matrix"))

  # Create Seurat object from *all* cells in filtered matrix
  seurat_all <- CreateSeuratObject(counts = filt_matrix, project = sample_name)

  # Perform minimal clustering (on all barcodes) for SoupX if you wish
  seurat_all <- NormalizeData(seurat_all)
  seurat_all <- FindVariableFeatures(seurat_all)
  seurat_all <- ScaleData(seurat_all)
  seurat_all <- RunPCA(seurat_all)
  seurat_all <- FindNeighbors(seurat_all, dims = 1:10)
  seurat_all <- FindClusters(seurat_all, resolution = 0.5)
  meta <- seurat_all@meta.data

  # Create SoupX channel
  soup_channel <- SoupChannel(tod = raw_matrix, toc = filt_matrix)

  # Robust barcode matching for clusters
  valid_barcodes <- intersect(rownames(meta), soupChannelCells(soup_channel))
  clusters_for_soupx <- setNames(meta$seurat_clusters[valid_barcodes], valid_barcodes)
  soup_channel <- setClusters(soup_channel, clusters_for_soupx)

  # Optionally, add DR coordinates
  # seurat_all <- RunUMAP(seurat_all, dims = 1:10)
  # umap <- seurat_all@reductions$umap@cell.embeddings
  # soup_channel <- setDR(soup_channel, umap[valid_barcodes, ])

  # Estimate and apply SoupX correction
  soup_channel <- autoEstCont(soup_channel)
  adj_matrix <- adjustCounts(soup_channel, roundToInt = TRUE)

  # Save the top ambient genes
  top_ambient <- head(soup_channel$soupProfile[order(soup_channel$soupProfile$est, decreasing = TRUE), ], n = 10)
  top_ambient$Sample <- sample_name

  # Save a Seurat object with the SoupX assay (no filtering yet!)
  seurat_all[["SoupX"]] <- CreateAssayObject(counts = adj_matrix)
  saveRDS(seurat_all, file = paste0("R_Data/", sample_name, "_seurat_SoupX.rds"))

  top_ambient
}

# Combine and write summary
summary_df <- bind_rows(top_ambient_genes)
write.csv(summary_df, "CSV_Results/DEGs_All/Ambient_genes_summary.csv", row.names = FALSE)

# Log the completion time
write(paste0("01b_soupx - Finish: ", Sys.time()), file = "scRNA_Log.txt", append = TRUE)
