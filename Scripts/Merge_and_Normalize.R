# Merge_and_Normalize.R
# Purpose: Merges Seurat objects from multiple samples and normalized the resulting object with SCTransform
# Author: Geoff Dilly

library(here)
library(Seurat)
library(patchwork)
library(ggplot2)
snRNA_home_dir <- here()
setwd(snRNA_home_dir)

# Log the start time and a timestamped copy of the script
write(paste0("Merge_and_Normalize - Start: ", Sys.time()), file = "snRNA_Log.txt", append = TRUE)
file.copy("Scripts/Merge_and_Normalize.R", paste0("Logs/Time_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", "Merge_and_Normalize.R"), overwrite = FALSE)

# Set 'R_MAX_VSIZE' to maximum RAM usage
Sys.setenv("R_MAX_VSIZE" = 32000000000)

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
#! mergenorm_combined_seurat <- merge, c(seurat_objects[[1]], seurat_objects[-1]))
mergenorm_combined_seurat <- merge(x = seurat_objects[[1]], y = seurat_objects[-1], add.cell.ids = str_sample_list)

saveRDS(mergenorm_combined_seurat, paste0("R_Data/", scConfig.Prefix, "_combined_merged.rds"))
rm(seurat_objects)
mergenorm_combined_seurat <- readRDS(paste0("R_Data/", scConfig.Prefix, "_combined_merged.rds"))

# Remove called doublets if you want to
if (scConfig.remove_doublets == TRUE) {
  mergenorm_combined_seurat <- subset(mergenorm_combined_seurat, subset = Doublet_Call == "Singlet")
}

# Remove mitochondrial genes if you want to
if (scConfig.remove_mito_genes == TRUE) {
  mergenorm_combined_seurat <- JoinLayers(mergenorm_combined_seurat)
  mito_genes <- grep(pattern = scConfig.mito_pattern, x = rownames(x = mergenorm_combined_seurat@assays$RNA$counts), value = TRUE)
  counts <- GetAssayData(mergenorm_combined_seurat, assay = "RNA", layer = "counts")
  counts <- counts[-(which(rownames(counts) %in% mito_genes)),]
  mergenorm_combined_seurat <- subset(mergenorm_combined_seurat, features = rownames(counts))
  mergenorm_combined_seurat[["RNA"]] <- split(mergenorm_combined_seurat[["RNA"]], f = mergenorm_combined_seurat$Sample_name)
  Layers(mergenorm_combined_seurat)
  rm(counts)
}

# Run SCTransform and save the normalized Seurat object
mergenorm_combined_seurat <- SCTransform(mergenorm_combined_seurat, verbose = TRUE, conserve.memory = TRUE)
saveRDS(mergenorm_combined_seurat, paste0("R_Data/", scConfig.Prefix, "_combined_SCT.rds"))

# Log the completion time
write(paste0("Merge_and_Normalize - Finish: ", Sys.time()), file = "snRNA_Log.txt", append = TRUE)
