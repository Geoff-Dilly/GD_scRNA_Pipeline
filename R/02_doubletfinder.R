# 02_doubletfinder.R
# Purpose: Run DoubletFinder on single sample Seurat objects
# Author: Geoff Dilly

library(here)
library(yaml)
library(Seurat)
library(stringr)
library(DoubletFinder)
library(doParallel)
library(foreach)

# Setup ####
# Load custom functions
source(here::here("R/modules/plot_utils.R"))
source(here::here("R/modules/log_utils.R"))
source(here::here("R/modules/qc_utils.R"))

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
write(paste0("02_doubletfinder - Start: ", Sys.time()), file = here::here("scRNA_Log.txt"), append = TRUE)
log_connection <- write_script_log(here::here("R/02_doubletfinder.R"))

# Log all output to the end of the log file
sink(log_connection, append = TRUE)
sink(log_connection, type = "message", append = TRUE)
on.exit({
  sink(NULL)
  sink(NULL, type = "message")
  stopCluster(cl)
})

# Run DoubletFinder ####
# Place each sample in a list for further processing
str_sample_list <- scConfig$Sample_metadata$Sample_name

# Run DoubletFinder in a parallel processing loop
foreach(sample_name = str_sample_list, .packages = c("Seurat", "DoubletFinder", "stringr")) %dopar% {
  sample_seurat <- readRDS(here::here("R_Data", paste0(sample_name, "_seurat.rds")))

  sample_seurat <- NormalizeData(sample_seurat)
  sample_seurat <- FindVariableFeatures(sample_seurat, selection.method = "vst", nfeatures = 2000)
  sample_seurat <- ScaleData(sample_seurat)
  sample_seurat <- RunPCA(sample_seurat)
  sample_seurat <- RunUMAP(sample_seurat, dims = 1:10)

  sweep_res_list <- paramSweep(sample_seurat, PCs = 1:10, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_res_list, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)

  pK_value <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric), drop = TRUE]))

  homotypic_prop <- modelHomotypic(sample_seurat@meta.data$seurat_clusters)
  nExp_poi <- round((scConfig$expct_doublet_pct/100) * nrow(sample_seurat@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic_prop))

  sample_seurat <- doubletFinder(sample_seurat,
                                 PCs = 1:10,
                                 pN = 0.25,
                                 pK = pK_value,
                                 nExp = nExp_poi,
                                 sct = FALSE)

  meta_cols <- colnames(sample_seurat@meta.data)
  score <- str_subset(meta_cols, "^pANN")
  call <- str_subset(meta_cols, "^DF.cl")
  sample_seurat$Doublet_Score <- sample_seurat[[score]]
  sample_seurat$Doublet_Call <- sample_seurat[[call]]
  sample_seurat[[call]] <- NULL
  sample_seurat[[score]] <- NULL

  # Remove large unnecessary assays
  sample_seurat[["RNA"]]$scale.data <- NULL
  sample_seurat[["RNA"]]$data <- NULL

  saveRDS(
    sample_seurat,
    file = here::here("R_Data", paste0(sample_name, "_seurat_Doublets.rds"))
  )
  NULL
}

# Log the completion time
write(paste0("02_doubletfinder - Finish: ", Sys.time()), file = here::here("scRNA_Log.txt"), append = TRUE)
