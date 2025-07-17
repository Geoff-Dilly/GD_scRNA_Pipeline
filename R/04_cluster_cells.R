# 04_cluster_cells.R
# Purpose: Performs the standard Seurat clustering workflow on a normalized Seurat object
# Author: Geoff Dilly

library(here)
library(dplyr)
library(Seurat)
library(ggplot2)
scRNA_home_dir <- here()
setwd(scRNA_home_dir)

# Setup ####
# Load custom functions
source("R/modules/plot_utils.R")
source("R/modules/log_utils.R")

# Check for required directories
check_required_dirs()

# Log the start time and a timestamped copy of the script
write(paste0("04_cluster_cells - Start: ", Sys.time()), file = "scRNA_Log.txt", append = TRUE)
log_file <- write_script_log("R/04_cluster_cells.R")

# Log all output to the end of the log file
sink(log_file, append = TRUE)
sink(log_file, type = "message", append = TRUE)
on.exit({
  sink(NULL)
  sink(NULL, type = "message")
})

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

# Remove exogenous genes from variable features
VariableFeatures(combined_seurat) <- setdiff(VariableFeatures(combined_seurat), scConfig.exogenous_genes)

# Scale data and run UMAP ####
# Perform PCA
combined_seurat <- RunPCA(combined_seurat, features = VariableFeatures(combined_seurat), npcs = 100, verbose = TRUE)

# Visualize the dimensionality of the PCs and pick the number of PCs
dimheatmaps_plot <- DimHeatmap(combined_seurat, dims = 1:9, cells = 500, balanced = TRUE)
save_plot_pdf(dimheatmaps_plot, "Plots/Quality_Control/DimHeatmaps.pdf", height = 10, width = 12)

elbow_plot <- ElbowPlot(combined_seurat, ndims = 100)
save_plot_pdf(elbow_plot, "Plots/Quality_Control/ElbowPlot.pdf", height = 4, width = 6)

# Perform UMAP dimensional reduction on the data
# Default is 25 dimensions
combined_seurat <- RunUMAP(combined_seurat, reduction = "pca", dims = 1:25)

# Examine the UMAP Plot for quality control and viability
qc_umap_plot <- FeaturePlot(combined_seurat, features = c("percent_mito", "nFeature_RNA", "Doublet_Score"), ncol = 3)
save_plot_pdf(qc_umap_plot, "Plots/Quality_Control/QC_UMAP.pdf", height = 4, width = 12)

# Perform clustering ####
# Identifies clusters of cells within the UMAP
combined_seurat <- FindNeighbors(combined_seurat, reduction = "pca", dims = 1:25)
combined_seurat <- FindClusters(combined_seurat, resolution = scConfig.clustering_resolution)

# Get the number of cells per each cluster
cells_per_cluster <- table(combined_seurat$seurat_clusters)
write.csv(cells_per_cluster, "CSV_Results/Cells_per_cluster.csv")

# Save the clustered Seurat object
saveRDS(combined_seurat, paste0("R_Data/", scConfig.Prefix, "_combined_clustered.rds"))

# Examine the resulting UMAP
clustered_umap_plot <- DimPlot(combined_seurat, label = TRUE)
save_plot_pdf(clustered_umap_plot, "Plots/Clustering_Plots/Clustered_UMAP.pdf", height = 4, width = 6)

# Visualize QC metrics in each cluster
percent_mito_vln_plot <- VlnPlot(combined_seurat, features = "percent_mito", pt.size = 0)
save_plot_pdf(percent_mito_vln_plot, "Plots/Quality_Control/percentMito_Violin.pdf", height = 4, width = 6)

nFeature_vln_plot <- VlnPlot(combined_seurat, features = "nFeature_RNA", pt.size = 0) # nolint
save_plot_pdf(nFeature_vln_plot, "Plots/Quality_Control/nFeature_Violin.pdf", height = 4, width = 6)

nCount_vln_plot <- VlnPlot(combined_seurat, features = "nCount_RNA", pt.size = 0) # nolint
save_plot_pdf(nCount_vln_plot, "Plots/Quality_Control/nCount_Violin.pdf", height = 4, width = 6)

doublet_score_vln_plot <- VlnPlot(combined_seurat, features = "Doublet_Score", pt.size = 0)
save_plot_pdf(doublet_score_vln_plot, "Plots/Quality_Control/Doublet_Score_Violin.pdf", height = 4, width = 6)

# Visualize marker gene expression ####
# Using markers from Dilly et al. 2022
marker_tbl <- read.csv("reference/marker_gene_db.csv", stringsAsFactors = FALSE)
marker_genes <- marker_tbl %>% filter(reference == scConfig.marker_gene_reference) %>% pull(gene)

major_cells_dotplot <- DotPlot(combined_seurat, features = marker_genes) + RotatedAxis()
save_plot_pdf(major_cells_dotplot, "Plots/Clustering_Plots/Major_Cell_Types_DotPlot.pdf", height = 6, width = 8)

# Log the completion time
write(paste0("04_cluster_cells - Finish: ", Sys.time()), file = "scRNA_Log.txt", append = TRUE)
