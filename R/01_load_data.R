# 01_load_data.R
# Purpose: Load a list of 10x objects into Seurat and save them as RDS
# Author: Geoff Dilly

library(here)
library(yaml)
library(Seurat)
library(dplyr)
library(foreach)
library(doParallel)
library(SoupX)

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

# Setup parallel backend
n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Log the start time and a timestamped copy of the script
write(paste0("01_load_data - Start: ", Sys.time()), file = here::here("scRNA_Log.txt"), append = TRUE)
log_connection <- write_script_log(here::here("R/01_load_data.R"))

# Log all output to the end of the log file
sink(log_connection, append = TRUE)
sink(log_connection, type = "message", append = TRUE)
on.exit({
  sink(NULL)
  sink(NULL, type = "message")
  stopCluster(cl)
})

# Check for required metadata columns
mandatory_metadata_columns <- c("Sample_name", "Treatment", "Sex", "Raw_data_dir")
for (col in mandatory_metadata_columns) {
  if (!col %in% colnames(scConfig$Sample_metadata)) {
    stop(paste("Missing mandatory metadata column:", col))
  }
}

# Load data and run SoupX (optional) ####
sample_list <- split(scConfig$Sample_metadata, seq_len(nrow(scConfig$Sample_metadata)))

top_ambient_genes <- foreach(sample = sample_list, .packages = c("Seurat", "SoupX", "dplyr")) %dopar% {
  top_ambient <- NULL

  # Ensure matrix directories are extracted and ready
  prepare_matrix_dir(
    sample$Raw_data_dir,
    if (dir.exists(file.path(sample$Raw_data_dir, "sample_filtered_feature_bc_matrix")))
      "sample_filtered_feature_bc_matrix"
    else
      "filtered_feature_bc_matrix"
  )

  if (scConfig$compute_soupx) {
    prepare_matrix_dir(
      sample$Raw_data_dir,
      if (dir.exists(file.path(sample$Raw_data_dir, "sample_raw_feature_bc_matrix")))
        "sample_raw_feature_bc_matrix"
      else
        "raw_feature_bc_matrix"
    )
  }

  # Load the 10x data into a Seurat object
filt_matrix <- Read10X(
  file.path(
    sample$Raw_data_dir,
    if (dir.exists(file.path(sample$Raw_data_dir, "sample_filtered_feature_bc_matrix")))
      "sample_filtered_feature_bc_matrix"
    else
      "filtered_feature_bc_matrix"
  )
)  sample_seurat <- CreateSeuratObject(counts = filt_matrix, project = scConfig$project_name, min.cells = 1, min.features = 1)

  # Compute SoupX assay (optional)
  if (scConfig$compute_soupx) {
    soupx_results <- run_soupx_correction(sample, sample_seurat)
    sample_seurat <- soupx_results$seurat_obj
    top_ambient <- soupx_results$top_ambient
  }

  # Filter the Seurat object for QC and add metadata
  sample_seurat <- PercentageFeatureSet(sample_seurat, pattern = scConfig$mito_pattern, col.name = "percent_mito")
  sample_seurat <- PercentageFeatureSet(sample_seurat, pattern = scConfig$ribo_pattern, col.name = "percent_ribo")
  sample_seurat <- subset(sample_seurat, subset = nFeature_RNA > scConfig$nFeature_RNA_cutoff &
                            percent_mito < scConfig$percent_mito_cutoff &
                            percent_ribo < scConfig$percent_ribo_cutoff)
  for (col in setdiff(colnames(sample), "Raw_data_dir")) {
    sample_seurat[[col]] <- sample[[col]]
  }
  saveRDS(
    sample_seurat,
    file = here::here("r_data", paste0(sample$Sample_name, "_seurat.rds"))
  )

  return(top_ambient)
}

if (scConfig$compute_soupx) {
  # Combine and write SoupX summary
  summary_df <- bind_rows(top_ambient_genes)
  write.csv(
    summary_df,
    here::here("csv_results", "Ambient_genes_summary.csv"),
    row.names = FALSE
  )
}

# Log the completion time
write(paste0("01_load_data - Finish: ", Sys.time()), file = here::here("scRNA_Log.txt"), append = TRUE)
