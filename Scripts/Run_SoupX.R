# Run_SoupX.R
# Purpose: Run SoupX on single sample Seurat objects
# Author: Geoff Dilly

library(here)
library(Seurat)
library(SoupX)
library(DropletUtils)
library(ggplot2)
library(DoubletFinder)
library(knitr)
snRNA_home_dir <- here()
setwd(snRNA_home_dir)

# Log the start time and a timestamped copy of the script
write(paste0("Run_SoupX - Start: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
file.copy("Scripts/Run_SoupX.R", paste0("Logs/Time_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", "Run_SoupX.R"), overwrite = FALSE)

# Code adapted from https://cellgeni.github.io/notebooks/html/new-10kPBMC-SoupX.html

# Load the configuration file and metadata
source("sc_experiment_config.R", local = TRUE)
scConfig.Sample_metadata <- read.csv("sc_sample_metadata.csv")

# Place each sample in a list for further processing 
str_sample_list <- scConfig.Sample_metadata$Sample_name


# Process each sample
for (sample_name in str_sample_list) {
  filt.matrix <- Read10X(paste0(scConfig.Raw_data_folder, "/", sample$Sample_name, "/", sample$Sample_name, "/outs/filtered_feature_bc_matrix"))
  raw.matrix <- Read10X(paste0(scConfig.Raw_data_folder, "/", sample$Sample_name, "/", sample$Sample_name, "/outs/raw_feature_bc_matrix"))
  str(raw.matrix)
  str(filt.matrix)

  # Create a Seurat object with the filtered data
  GD_10x_Sample <- CreateSeuratObject(counts = filt.matrix)

  # Create the SoupX channel
  soup.channel <- SoupChannel(raw.matrix, filt.matrix)

  # Cluster the Seurat object
  GD_10x_Sample <- SCTransform(GD_10x_Sample, verbose = T)
  GD_10x_Sample <- RunPCA(GD_10x_Sample, verbose = T)
  GD_10x_Sample <- RunUMAP(GD_10x_Sample, dims = 1:30, verbose = T)
  GD_10x_Sample <- FindNeighbors(GD_10x_Sample, dims = 1:30, verbose = T)
  GD_10x_Sample <- FindClusters(GD_10x_Sample, verbose = T)

  # Extract metadata and UMAP coordinates
  meta <- GD_10x_Sample@meta.data
  umap <- GD_10x_Sample@reductions$umap@cell.embeddings

  # Assign clustering to SoupX
  soup.channel <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup.channel <- setDR(soup.channel, umap)
  head(meta)

  # Estimate and apply the SoupX correction
  soup.channel <- autoEstCont(soup.channel)
  adj.matrix <- adjustCounts(soup.channel, roundToInt = T)
  
  # Display the top ambient RNA
  head(soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ], n = 50)
  
  # Save the corrected matrix
  write10xCounts(paste0("Raw_Data/", sample_name, "_filtered_feature_bc_matrix_SoupX_adjusted.h5"), adj.matrix)

  # Optional visualization for SoupX correction
  #!topsoupgenes <- head(soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ], n = 100)
  #!neatgene <- "Apoe"
  #!soupchangeplot <- plotChangeMap(soup.channel, adj.matrix, neatgene)
  #!soupchangeplot
}

# Log the completion time
write(paste0("Run_SoupX - Finish: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
