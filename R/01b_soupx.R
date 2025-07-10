# 01b_soupx.R
# Purpose: Run SoupX on single sample Seurat objects
# Author: Geoff Dilly

library(here)
library(Seurat)
library(SoupX)
library(DropletUtils)
library(foreach)
library(doParallel)
snRNA_home_dir <- here()
setwd(snRNA_home_dir)

# Log the start time and a timestamped copy of the script
write(paste0("Run_SoupX - Start: ", Sys.time()), file = "snRNA_Log.txt", append = TRUE)
file.copy("Scripts/Run_SoupX.R", paste0("Logs/Time_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", "Run_SoupX.R"), overwrite = FALSE)

# Code adapted from https://cellgeni.github.io/notebooks/html/new-10kPBMC-SoupX.html

# Load the configuration file and metadata
source("sc_experiment_config.R", local = TRUE)
scConfig.Sample_metadata <- read.csv("sc_sample_metadata.csv")

# Setup parallel backend
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores/2)
registerDoParallel(cl)

# Place each sample in a list for further processing
str_sample_list <- scConfig.Sample_metadata$Sample_name

# Process each sample
foreach(sample_name = str_sample_list, .packages = c("Seurat", "SoupX", "DropletUtils")) %dopar% {
  filt_matrix <- Read10X(paste0(scConfig.Raw_data_folder, "/", sample_name, "/", sample_name, "/outs/filtered_feature_bc_matrix"))
  raw_matrix <- Read10X(paste0(scConfig.Raw_data_folder, "/", sample_name, "/", sample_name, "/outs/raw_feature_bc_matrix"))
  str(raw_matrix)
  str(filt_matrix)

  # Create a Seurat object with the filtered data
  sample_seurat <- CreateSeuratObject(counts = filt_matrix)

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
  head(meta)

  # Estimate and apply the SoupX correction
  soup_channel <- autoEstCont(soup_channel)
  adj_matrix <- adjustCounts(soup_channel, roundToInt = TRUE)

  # Display the top ambient RNA
  head(soup_channel$soupProfile[order(soup_channel$soupProfile$est, decreasing = TRUE), ], n = 10)

  # Save the corrected matrix
  write10xCounts(paste0("Raw_Data/", sample_name, "_filtered_feature_bc_matrix_SoupX_adjusted.h5"), adj_matrix)

  NULL
}

# Log the completion time
write(paste0("Run_SoupX - Finish: ", Sys.time()), file = "snRNA_Log.txt", append = TRUE)
