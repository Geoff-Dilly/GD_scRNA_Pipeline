# 02_doubletfinder.R
# Purpose: Run DoubletFinder on single sample Seurat objects
# Author: Geoff Dilly

library(here)
library(Seurat)
library(stringr)
library(DoubletFinder)
library(doParallel)
library(foreach)
snRNA_home_dir <- here()
setwd(snRNA_home_dir)

# Load custom functions
source("R/modules/plot_utils.R")

# Log the start time and a timestamped copy of the script
write(paste0("02_doubletfinder - Start: ", Sys.time()), file = "snRNA_Log.txt", append = TRUE)
file.copy("R/02_doubletfinder.R", paste0("Logs/Time_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", "02_doubletfinder.R"), overwrite = FALSE)

# Read the sample metadata file
scConfig.Sample_metadata <- read.csv("sc_sample_metadata.csv")

# Setup parallel backend
n_cores <- parallel::detectCores() - 1  # Or set n_cores manually if desired
cl <- makeCluster(n_cores)
registerDoParallel(cl)

str_sample_list <- scConfig.Sample_metadata$Sample_name

# Run doubletFinder --------------------
# Parallel processing loop
foreach(sample_name = str_sample_list, .packages = c("Seurat", "DoubletFinder", "stringr")) %dopar% {
  sample_seurat <- LoadSeuratRds(paste0("R_Data/", sample_name, "_seurat.rds"))

  sample_seurat <- NormalizeData(sample_seurat)
  sample_seurat <- FindVariableFeatures(sample_seurat, selection.method = "vst", nfeatures = 2000)
  sample_seurat <- ScaleData(sample_seurat)
  sample_seurat <- RunPCA(sample_seurat)
  sample_seurat <- RunUMAP(sample_seurat, dims = 1:10)

  sweep_res_list <- paramSweep(sample_seurat, PCs = 1:10, sct = FALSE)
  sweep_stats <- summarizeSweep(sweep_res_list, GT = FALSE)
  bcmvn <- find.pK(sweep_stats)

  pK_value <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric), drop = TRUE]))

  homotypic_prop <- modelHomotypic(sample_seurat@meta.data$Sample_name)
  nExp_poi <- round(0.075 * nrow(sample_seurat@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic_prop))

  sample_seurat <- doubletFinder(sample_seurat, PCs = 1:10, pN = 0.25, pK = pK_value, nExp = nExp_poi, sct = FALSE)

  meta_cols <- colnames(sample_seurat@meta.data)
  score <- str_subset(meta_cols, "^pANN")
  call <- str_subset(meta_cols, "^DF.cl")
  sample_seurat$Doublet_Score <- sample_seurat[[score]]
  sample_seurat$Doublet_Call <- sample_seurat[[call]]
  sample_seurat[[call]] <- NULL
  sample_seurat[[score]] <- NULL

  sample_seurat[["RNA"]]$scale.data <- NULL
  sample_seurat[["RNA"]]$data <- NULL

  saveRDS(sample_seurat, file = paste0("R_Data/", sample_name, "_seurat_Doublets.rds"))
  NULL
}

stopCluster(cl)

# Log the completion time
write(paste0("02_doubletfinder - Finish: ", Sys.time()), file = "snRNA_Log.txt", append = TRUE)
