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
    # Load the unfiltered matrix for SoupX
    raw_matrix <- Read10X(file.path(sample$Raw_data_dir, "raw_feature_bc_matrix"))

    # Minimally cluster the seurat object
    sample_seurat <- SCTransform(sample_seurat, verbose = FALSE)
    sample_seurat <- RunPCA(sample_seurat, verbose = FALSE)
    sample_seurat <- FindNeighbors(sample_seurat, dims = 1:10)
    sample_seurat <- FindClusters(sample_seurat)

    # Run SoupX
    soup_channel <- SoupChannel(tod = raw_matrix, toc = filt_matrix)

    # Assign clusters to SoupX
    cluster_labels <- sample_seurat$seurat_clusters
    cell_names     <- rownames(sample_seurat@meta.data)
    soup_channel   <- setClusters(soup_channel, setNames(cluster_labels, cell_names))

    soup_channel <- autoEstCont(soup_channel)
    adj_matrix   <- adjustCounts(soup_channel, roundToInt = TRUE)

    # Add the SoupX results as a new assay
    sample_seurat[["SoupX"]] <- CreateAssayObject(counts = adj_matrix)

    # Save the top ambient genes
    top_ambient <- head(soup_channel$soupProfile[order(soup_channel$soupProfile$est, decreasing = TRUE), ], n = 10)
    top_ambient$Sample <- sample$Sample_name
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

  top_ambient
}

stopCluster(cl)

if (scConfig.compute_soupx == TRUE) {
  # Combine and write SoupX summary
  summary_df <- bind_rows(top_ambient_genes)
  write.csv(summary_df, "CSV_Results/DEGs_All/Ambient_genes_summary.csv", row.names = FALSE)
}

# Log the completion time
write(paste0("01_load_data - Finish: ", Sys.time()), file = "scRNA_Log.txt", append = TRUE)
