# 03_normalize_and_integrate.R
# Purpose: Merges Seurat objects from multiple samples and normalizes the resulting object with SCTransform
# Author: Geoff Dilly

library(here)
library(yaml)
library(Seurat)

# Setup ####
# Load custom functions
source(here::here("R/modules/run_utils.R"))
source(here::here("R/modules/qc_utils.R"))
source(here::here("R/modules/plot_utils.R"))

# Load the configuration file and metadata
scConfig <- yaml::read_yaml(here::here("sc_experiment_config.yaml"))
scConfig$Sample_metadata <- read.csv(here::here("sc_sample_metadata.csv"))

# Check for required directories
check_required_dirs()

## Log the start time and a timestamped copy of the script
write(paste0("03_normalize_and_integrate - Start: ", Sys.time()), file = here::here("scRNA_Log.txt"), append = TRUE)
log_connection <- write_script_log(here::here("R/03_normalize_and_integrate.R"))

# Log all output to the end of the log file
sink(log_connection, append = TRUE)
sink(log_connection, type = "message", append = TRUE)
on.exit({
  sink(NULL)
  sink(NULL, type = "message")
})

# Set 'R_MAX_VSIZE' to maximum RAM usage
Sys.setenv("R_MAX_VSIZE" = 32000000000)

# Get sample names in a list of strings
str_sample_list <- scConfig$Sample_metadata$Sample_name

# Initialize a list to store the Seurat objects
seurat_objects <- list()

# Filter data ####
for (sample in str_sample_list) {
  sample_seurat <- readRDS(here::here("R_Data", paste0(sample, "_seurat_Doublets.rds")))

  # Filter cells and genes based on scConfig params
  sample_seurat <- qc_filter_seurat(sample_seurat)

  seurat_objects[[sample]] <- sample_seurat
}

# Merge and normalize with SCTransform ####
# Merge the Seurat objects
seurat_combined <- merge(x = seurat_objects[[1]], y = seurat_objects[-1])
rm(seurat_objects)

# Perform SCTransform and PCA
seurat_combined <- SCTransform(seurat_combined,
                               assay = "RNA",
                               verbose = TRUE,
                               conserve.memory = TRUE,
                               return.only.var.genes = TRUE)

seurat_combined <- RunPCA(seurat_combined)

# Save the object with SCT normalization and PCA
sct_obj_path <- here::here("R_Data", paste0(scConfig$prefix, "_object_norm.rds"))
saveRDS(seurat_combined, sct_obj_path)

# Integrate the samples ####
seurat_combined <- IntegrateLayers(object = seurat_combined, 
                                   method = RPCAIntegration, 
                                   orig.reduction = "pca", 
                                   new.reduction = "integrated.rpca",
                                   assay = "SCT",
                                   normalization.method = "SCT",
                                   verbose = TRUE)

saveRDS(seurat_combined, here::here("R_Data", paste0(scConfig$prefix, "_SCT_integrated.rds")))

# Examine QC metrics by animal ####
Idents(seurat_combined) <- seurat_combined$Sample_name

nFeature_vln_by_animal <- VlnPlot(seurat_combined, features = "nFeature_RNA", pt.size = 0)
save_plot_pdf(nFeature_vln_by_animal, here::here("Plots/Quality_Control", "nFeature_ViolinPlot_byAnimal.pdf"), height = 4, width = 6)

nCount_vln_by_animal <- VlnPlot(seurat_combined, features = "nCount_RNA", pt.size = 0)
save_plot_pdf(nCount_vln_by_animal, here::here("Plots/Quality_Control", "nCount_ViolinPlot_byAnimal.pdf"), height = 4, width = 6)

percent_mito_vln_by_animal <- VlnPlot(seurat_combined, features = "percent_mito", pt.size = 0)
save_plot_pdf(percent_mito_vln_by_animal, here::here("Plots/Quality_Control", "percentMito_ViolinPlot_byAnimal.pdf"), height = 4, width = 6)

# Log the completion time
write(paste0("03_normalize_and_integrate - Finish: ", Sys.time()), file = here::here("scRNA_Log.txt"), append = TRUE)
