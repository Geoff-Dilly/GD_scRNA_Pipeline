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

# Function to get results directory based on run time 
output_dir <- get_results_dir(run_time = Sys.getenv("RUN_TIME"), 
                            prefix = scConfig$prefix)

## Log the start time and a timestamped copy of the script
write(paste0("03_normalize_and_integrate - Start: ", Sys.time()), file = here::here("scRNA_Log.txt"), append = TRUE)
log_connection <- write_script_log(here::here("R/03_normalize_and_integrate.R"), log_dir = here::here(output_dir, "logs"))

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
  sample_seurat <- readRDS(here::here("data", "R_Data", paste0(sample, "_seurat_Doublets.rds")))

  # Filter cells and genes based on scConfig params
  sample_seurat <- qc_filter_seurat(sample_seurat)

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

saveRDS(integrated_seurat, here::here("data", "R_Data", paste0(scConfig$prefix, "_SCT_integrated.rds")))

# Examine QC metrics by animal ####
Idents(integrated_seurat) <- integrated_seurat$Sample_name

nFeature_vln_by_animal <- VlnPlot(integrated_seurat, features = "nFeature_RNA", pt.size = 0)
save_plot_pdf(nFeature_vln_by_animal, here::here(output_dir, "plots", "Quality_Control", "nFeature_ViolinPlot_byAnimal.pdf"), height = 4, width = 6)

nCount_vln_by_animal <- VlnPlot(integrated_seurat, features = "nCount_RNA", pt.size = 0)
save_plot_pdf(nCount_vln_by_animal, here::here(output_dir, "plots", "Quality_Control", "nCount_ViolinPlot_byAnimal.pdf"), height = 4, width = 6)

percent_mito_vln_by_animal <- VlnPlot(integrated_seurat, features = "percent_mito", pt.size = 0)
save_plot_pdf(percent_mito_vln_by_animal, here::here(output_dir, "plots", "Quality_Control", "percentMito_ViolinPlot_byAnimal.pdf"), height = 4, width = 6)

# Log the completion time
write(paste0("03_normalize_and_integrate - Finish: ", Sys.time()), file = here::here("scRNA_Log.txt"), append = TRUE)
