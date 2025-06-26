# Merge_and_Normalize.R
# Purpose: Merges Seurat objects from multiple samples and normalized the resulting object with SCTransform
# Author: Geoff Dilly

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
snRNA_home_dir <- "/Users/gad479/Documents/GD_scRNA_Pipeline"
setwd(snRNA_home_dir)

# Log the start time and a timestamped copy of the script
write(paste0("Merge_and_Normalize - Start: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
file.copy("Scripts/Merge_and_Normalize.R", paste0("Logs/Time_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", "Merge_and_Normalize.R"), overwrite = FALSE)

# Set 'R_MAX_VSIZE' to maximum RAM usage
Sys.setenv('R_MAX_VSIZE'=32000000000)

# Load the configuration file and metadata
source("sc_experiment_config.R")
scConfig.Sample_metadata <- read.csv("sc_sample_metadata.csv")

# Get sample names in a list of strings
str_sample_list <- scConfig.Sample_metadata$Sample_name

# Initialize a list to store the Seurat objects
seurat_objects <- list()

# Read RDS files and assign them to variables dynamically
for (sample in str_sample_list) {
  sample_seurat <- readRDS(paste0("R_Data/", sample, "_seurat_Doublets.rds"))
  seurat_objects[[sample]] <- sample_seurat
}

# Merge samples into one Seurat object without integration
# Convert the list of Seurat objects to a list of arguments for the merge function
#! MergeNorm_Combined_Seurat <- merge, c(seurat_objects[[1]], seurat_objects[-1]))
MergeNorm_Combined_Seurat <- merge(x = seurat_objects[[1]], y = seurat_objects[-1], add.cell.ids = str_sample_list)

saveRDS(MergeNorm_Combined_Seurat, paste0("R_Data/",scConfig.Prefix ,"_combined_merged.rds"))
rm(seurat_objects)
MergeNorm_Combined_Seurat <- readRDS(paste0("R_Data/",scConfig.Prefix ,"_combined_merged.rds"))

# Remove called doublets if you want to
if(scConfig.remove_doublets == TRUE) {
  MergeNorm_Combined_Seurat <- subset(MergeNorm_Combined_Seurat, subset = Doublet_Call == "Singlet")}

# Remove mitochondrial genes if you want to
if(scConfig.remove_mito_genes == TRUE) {
  MergeNorm_Combined_Seurat <- JoinLayers(MergeNorm_Combined_Seurat)
  mito.genes <- grep(pattern = scConfig.mito_pattern, x = rownames(x = MergeNorm_Combined_Seurat@assays$RNA$counts), value = TRUE)
  counts <- GetAssayData(MergeNorm_Combined_Seurat, assay = "RNA", layer = "counts")
  counts <- counts[-(which(rownames(counts) %in% mito.genes)),]
  MergeNorm_Combined_Seurat <- subset(MergeNorm_Combined_Seurat, features = rownames(counts))
  MergeNorm_Combined_Seurat[["RNA"]] <- split(MergeNorm_Combined_Seurat[["RNA"]], f = MergeNorm_Combined_Seurat$Sample_name)
  Layers(MergeNorm_Combined_Seurat)
  rm(counts)
}

# Run SCTransform and save the normalized Seurat object
MergeNorm_Combined_Seurat <- SCTransform(MergeNorm_Combined_Seurat, verbose = TRUE, conserve.memory=TRUE)
saveRDS(MergeNorm_Combined_Seurat, paste0("R_Data/",scConfig.Prefix ,"_combined_SCT.rds"))

# Log the completion time
write(paste0("Merge_and_Normalize - Finish: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
