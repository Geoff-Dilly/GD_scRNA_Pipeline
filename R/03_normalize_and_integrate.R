# 03_normalize_and_integrate.R
# Purpose: Merges Seurat objects from multiple samples and normalized the resulting object with SCTransform
# Author: Geoff Dilly

library(here)
library(Seurat)
scRNA_home_dir <- here()
setwd(scRNA_home_dir)

# Setup ####
# Load custom functions
source("R/modules/log_utils.R")

## Log the start time and a timestamped copy of the script
write(paste0("03_normalize_and_integrate - Start: ", Sys.time()), file = "scRNA_Log.txt", append = TRUE)
write_script_log("R/03_normalize_and_integrate.R")

# Set 'R_MAX_VSIZE' to maximum RAM usage
Sys.setenv("R_MAX_VSIZE" = 32000000000)

# Load the configuration file and metadata
source("sc_experiment_config.R")
scConfig.Sample_metadata <- read.csv("sc_sample_metadata.csv")

# Get sample names in a list of strings
str_sample_list <- scConfig.Sample_metadata$Sample_name

# Initialize a list to store the Seurat objects
seurat_objects <- list()

# Filter data ####
# Read RDS files and assign them to variables dynamically
for (sample in str_sample_list) {
  sample_seurat <- readRDS(paste0("R_Data/", sample, "_seurat_Doublets.rds"))

  if (scConfig.soupx_adjust == TRUE) {
    default_assay(sample_seurat) <- "SoupX"
  }

  # Remove called doublets if specified
  if (scConfig.remove_doublets == TRUE) {
    sample_seurat <- subset(sample_seurat, subset = Doublet_Call == "Singlet")
  }

  # Remove mitochondrial genes if specified
  if (scConfig.remove_mito_genes == TRUE) {
    mito_genes <- grep(scConfig.mito_pattern, rownames(sample_seurat), value = TRUE)
    sample_seurat <- subset(
      sample_seurat,
      features = setdiff(rownames(sample_seurat), mito_genes)
    )
  }

  # Remove ribosomal genes if specified
  if (scConfig.remove_ribo_genes == TRUE) {
    ribo_genes <- grep(scConfig.ribo_pattern, rownames(sample_seurat), value = TRUE)
    sample_seurat <- subset(
      sample_seurat,
      features = setdiff(rownames(sample_seurat), ribo_genes)
    )
  }

  # Remove the top quartile high UMI cells if specified
  if (scConfig.remove_top_nUMIs == TRUE) {
    umi_threshold <- quantile(sample_seurat$nCount_RNA, 0.75)
    sample_seurat <- subset(sample_seurat, subset = nCount_RNA < umi_threshold)
  }

  seurat_objects[[sample]] <- sample_seurat
}

# Normalize with SCTransform ####
seurat_objects <- lapply(seurat_objects, function(obj) {
  SCTransform(obj, assay = "RNA", verbose = TRUE)
})

# Integrate the datasets ####
# Select integration features
features <- SelectIntegrationFeatures(object.list = seurat_objects, nfeatures = 3000)

# Prep objects for integration
seurat_objects <- PrepSCTIntegration(object.list = seurat_objects, anchor.features = features)
anchors <- FindIntegrationAnchors(
  object.list = seurat_objects,
  normalization.method = "SCT",
  anchor.features = features
)

# Integrate the datasets
integrated_seurat <- IntegrateData(
  anchorset = anchors,
  normalization.method = "SCT"
)

saveRDS(integrated_seurat, paste0("R_Data/", scConfig.Prefix, "_SCT_integrated.rds"))

# Log the completion time
write(paste0("03_normalize_and_integrate - Finish: ", Sys.time()), file = "scRNA_Log.txt", append = TRUE)
