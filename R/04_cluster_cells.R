# 04_cluster_cells.R
# Purpose: Performs the standard Seurat clustering workflow on a normalized Seurat object
# Author: Geoff Dilly

library(here)
library(dplyr)
library(Seurat)
library(ggplot2)
snRNA_home_dir <- here()
setwd(snRNA_home_dir)

# Load custom functions
source("R/modules/plot_utils.R")

# Log the start time and a timestamped copy of the script
write(paste0("04_cluster_cells - Start: ", Sys.time()), file = "snRNA_Log.txt", append = TRUE)
file.copy("R/04_cluster_cells.R", paste0("Logs/Time_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", "04_cluster_cells.R"), overwrite = FALSE)

# Load the configuration file and metadata
source("sc_experiment_config.R")
scConfig.Sample_metadata <- read.csv("sc_sample_metadata.csv")

# Load normalized data
combined_seurat <- readRDS(paste0("R_Data/", scConfig.Prefix, "_SCT_integrated.rds"))

# Get the number of cells per each sample
cells_per_sample <- table(combined_seurat$Sample_name)
write.csv(cells_per_sample, "CSV_Results/Cells_per_sample.csv")

# Check quality metrics for each cell
qc_violins_plot <- VlnPlot(combined_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "Doublet_Score"), ncol = 4, pt.size = 0)
save_plot_pdf(qc_violins_plot, "Plots/Quality_Control/QC_Violins.pdf", height = 4, width = 12)

# Scale data and run UMAP ----------------------
# Perform PCA
combined_seurat <- RunPCA(combined_seurat, npcs = 100, verbose = TRUE)

# Visualize the dimensionality of the PCs and pick the number of PCs
dimheatmaps_plot <- DimHeatmap(combined_seurat, dims = 1:9, cells = 500, balanced = TRUE)
save_plot_pdf(dimheatmaps_plot, "Plots/Quality_Control/DimHeatmaps.pdf", height = 4, width = 6)

elbow_plot <- ElbowPlot(combined_seurat, ndims = 100)
save_plot_pdf(elbow_plot, "Plots/Quality_Control/ElbowPlot.pdf", height = 4, width = 6)

# Perform UMAP dimensional reduction on the data
# Default is 25 dimensions
combined_seurat <- RunUMAP(combined_seurat, reduction = "pca", dims = 1:25)

# Examine the resulting UMAP-------------
raw_umap_plot <- DimPlot(combined_seurat)
save_plot_pdf(raw_umap_plot, "Plots/Quality_Control/Raw_UMAP.pdf", height = 4, width = 6)

# Examine the UMAP Plot for quality control and viability
qc_umap_plot <- FeaturePlot(combined_seurat, features = c("percent_mito", "nFeature_RNA", "Doublet_Score"), ncol = 3)
save_plot_pdf(qc_umap_plot, "Plots/Quality_Control/QC_UMAP.pdf", height = 4, width = 12)

# Perform clustering -----------------------
# Identifies clusters of cells within the UMAP
combined_seurat <- FindNeighbors(combined_seurat, reduction = "pca", dims = 1:25)
combined_seurat <- FindClusters(combined_seurat, resolution = scConfig.clustering_resolution)

# Get the number of cells per each cluster
cells_per_cluster <- table(combined_seurat$seurat_clusters)
write.csv(cells_per_cluster, "CSV_Results/Cells_per_cluster.csv")

# Save the clustered Seurat object
saveRDS(combined_seurat, paste0("R_Data/", scConfig.Prefix, "_combined_clustered.rds"))

# Examine the resulting UMAP-------------
clustered_umap_plot <- DimPlot(combined_seurat, label = TRUE)
save_plot_pdf(clustered_umap_plot, "Plots/Clustering_Plots/Clustered_UMAP.pdf", height = 4, width = 6)

# Visualize QC metrics in each cluster
percent_mito_vln_plot <- VlnPlot(combined_seurat, features = "percent_mito", pt.size = 0)
save_plot_pdf(percent_mito_vln_plot, "Plots/Quality_Control/percentMito_Violin.pdf", height = 4, width = 6)

nFeature_vln_plot <- VlnPlot(combined_seurat, features = "nFeature_RNA", pt.size = 0) # nolint
save_plot_pdf(nFeature_vln_plot, "Plots/Clustering_Plots/nFeature_Violin.pdf", height = 4, width = 6)

nCount_vln_plot <- VlnPlot(combined_seurat, features = "nCount_RNA", pt.size = 0) # nolint
save_plot_pdf(nCount_vln_plot, "Plots/Clustering_Plots/nCount_Violin.pdf", height = 4, width = 6)

doublet_score_vln_plot <- VlnPlot(combined_seurat, features = "Doublet_Score", pt.size = 0)
save_plot_pdf(doublet_score_vln_plot, "Plots/Clustering_Plots/Doublet_Score_Violin.pdf", height = 4, width = 6)

# Visualize marker gene expression in each cluster
# Using markers from Dilly et al. 2022
id_features <- c("Mbp", "Mobp", "Plp1", "Gad1", "Gad2",
                 "Ndrg2", "Slc1a2", "Slc4a4",
                 "Slc17a7", "Satb1", "Neurod6", "Vcan",
                 "Pdgfra", "Pcdh15", "Csf1r", "Apbb1ip", "P2ry12",
                 "Flt1", "B2m", "Bmp4", "Cnp", "Ccdc153",
                 "Rsph1", "Tmem212", "Rbfox3")

major_cells_dotplot <- DotPlot(combined_seurat, features = id_features) + RotatedAxis()
save_plot_pdf(major_cells_dotplot, "Plots/Clustering_Plots/Major_Cell_Types_DotPlot.pdf", height = 6, width = 8)

# Examine QC metrics by animal
Idents(combined_seurat) <- combined_seurat$Sample_name

nFeature_vln_by_animal <- VlnPlot(combined_seurat, features = "nFeature_RNA", pt.size = 0) # nolint
save_plot_pdf(nFeature_vln_by_animal, "Plots/Quality_Control/nFeature_ViolinPlot_byAnimal.pdf", height = 4, width = 6)

nCount_vln_by_animal <- VlnPlot(combined_seurat, features = "nCount_RNA", pt.size = 0) # nolint
save_plot_pdf(nCount_vln_by_animal, "Plots/Quality_Control/nCount_ViolinPlot_byAnimal.pdf", height = 4, width = 6)

percent_mito_vln_by_animal <- VlnPlot(combined_seurat, features = "percent_mito", pt.size = 0)
save_plot_pdf(percent_mito_vln_by_animal, "Plots/Quality_Control/percentMito_ViolinPlot_byAnimal.pdf", height = 4, width = 6)

# Log the completion time
write(paste0("04_cluster_cells - Finish: ", Sys.time()), file = "snRNA_Log.txt", append = TRUE)
