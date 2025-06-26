# Run_DoubletFinder.R
# Purpose: Run DoubletFinder on single sample Seurat objects
# Author: Geoff Dilly

library(here)
library(Seurat)
library(stringr)
library(DoubletFinder)
library(data.table)
snRNA_home_dir <- here()
setwd(snRNA_home_dir)

# Log the start time and a timestamped copy of the script
write(paste0("Run_DoubletFinder - Start: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
file.copy("Scripts/Run_DoubletFinder.R", paste0("Logs/Time_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", "Run_DoubletFinder.R.R"), overwrite = FALSE)

# Read the sample metadata file
scConfig.Sample_metadata <- read.csv("sc_sample_metadata.csv")

str_sample_list <- scConfig.Sample_metadata$Sample_name

# Run doubletFinder ---------------------------------------------------------------------------------------
for (sample_name in str_sample_list) {
  # Pre-process Seurat object (standard)
  GD_10x_Sample <- LoadSeuratRds(paste0("R_Data/", sample_name, "_seurat.rds"))

  GD_10x_Sample <- NormalizeData(GD_10x_Sample)
  GD_10x_Sample <- FindVariableFeatures(GD_10x_Sample, selection.method = "vst", nfeatures = 2000)
  GD_10x_Sample <- ScaleData(GD_10x_Sample)
  GD_10x_Sample <- RunPCA(GD_10x_Sample)
  GD_10x_Sample <- RunUMAP(GD_10x_Sample, dims = 1:10)

  # pK Identification (no ground-truth)
  sweep.res.list <- paramSweep(GD_10x_Sample, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)

  pK_value <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric), drop=T]))
  print(pK_value)

  # Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(GD_10x_Sample@meta.data$Sample_name)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.075*nrow(GD_10x_Sample@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  # Run DoubletFinder
  GD_10x_Sample <- doubletFinder(GD_10x_Sample, PCs = 1:10, pN = 0.25, pK = pK_value, nExp = nExp_poi, sct = FALSE)
  
  # Rename the doubletFinder results
  meta_cols <- colnames(GD_10x_Sample@meta.data)
  score <- str_subset(meta_cols, "^pANN")
  call <- str_subset(meta_cols, "^DF.cl")
  GD_10x_Sample$Doublet_Score <- GD_10x_Sample[[score]]
  GD_10x_Sample$Doublet_Call <- GD_10x_Sample[[call]]
  GD_10x_Sample[[call]] <- NULL
  GD_10x_Sample[[score]] <- NULL
  
  # Remove unnecessary data layers
  GD_10x_Sample[["RNA"]]$scale.data <- NULL
  GD_10x_Sample[["RNA"]]$data <- NULL
  
  # Save the Seurat object with doubletFinder Results
  saveRDS(GD_10x_Sample, file = paste0("R_Data/", sample_name, "_seurat_Doublets.rds"))

}

# Examine the results
#!CeA_sample <- LoadSeuratRds("/Volumes/users/geoff_scratch/astrocyte_reanalysis/data/r_data/GD_2A_Alcohol_doublets.rds")
#!DimPlot(CeA_sample)
#!FeaturePlot(CeA_sample, features = "pANN_0.25_0.24_709", cols = c("blue", "yellow"), max.cutoff = 0.8)
#!FeaturePlot(CeA_sample, features = "DF.classifications_0.25_0.24_709")
#!doublet_freqs <- as.data.frame(table(CeA_sample@meta.data$DF.classifications_0.25_0.24_709))
#!print(doublet_freqs)
#!percentage_doublets <- doublet_freqs$Freq[1] / (doublet_freqs$Freq[2] + doublet_freqs$Freq[1])*100
#!percentage_doublets

# Log the completion time
write(paste0("Run_DoubletFinder - Finish: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
