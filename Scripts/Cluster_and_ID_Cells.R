# Cluster_and_ID_Cells.R
# Purpose: Performs the standard Seurat clustering workflow on a normalized Seurat object
# Author: Geoff Dilly

library(here)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
snRNA_home_dir <- here()
setwd(snRNA_home_dir)

# Log the start time and a timestamped copy of the script
write(paste0("Cluster_and_ID_Cells - Start: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
file.copy("Scripts/Cluster_and_ID_Cells.R", paste0("Logs/Time_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", "Cluster_and_ID_Cells.R"), overwrite = FALSE)

# Load the configuration file and metadata
source("sc_experiment_config.R")
scConfig.Sample_metadata <- read.csv("sc_sample_metadata.csv")

# Load normalized data
Combined_Seurat <- readRDS(paste0("R_Data/",scConfig.Prefix ,"_combined_SCT.rds"))

# Get the number of cells per each sample
cells_per_sample <- table(Combined_Seurat$Sample_name)
write.csv(cells_per_sample,"CSV_Results/Cells_per_sample.csv")

# Check quality metrics for each cell
QC_Violins_Plot <- VlnPlot(Combined_Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "Doublet_Score"), ncol = 4, pt.size = 0)
pdf("Plots/Quality_Control/QC_VlnPlot.pdf", height = 6, width = 10)
  print(QC_Violins_Plot)
  dev.off()

# Scale data and run UMAP ----------------------
# Perform PCA
Combined_Seurat <- RunPCA(Combined_Seurat, npcs = 100 ,verbose = TRUE)

# Visualize the dimensionality of the PCs and pick the number of PCs
DimHeatmaps_Plot <- DimHeatmap(Combined_Seurat, dims = 1:9, cells = 500, balanced = TRUE)
pdf("Plots/Quality_Control/Dim_Heatmap.pdf", height = 6, width = 8)
  print(DimHeatmaps_Plot)
  dev.off()

Elbow_Plot <- ElbowPlot(Combined_Seurat, ndims = 100)
pdf("Plots/Quality_Control/Elbowplot.pdf")
  print(Elbow_Plot)
  dev.off()
# I chose 25 PCs

# Perform UMAP dimensional reduction on the data
Combined_Seurat <- RunUMAP(Combined_Seurat, reduction = "pca", dims = 1:25)

# Examine the resulting UMAP-------------
Raw_UMAP_Plot <- DimPlot(Combined_Seurat)
pdf("Plots/Clustering_Plots/Raw_UMAP.pdf", height = 4, width = 6)
  print(Raw_UMAP_Plot)
  dev.off()

# Examine the UMAP Plot for quality control and viability
QC_UMAP_Plot <- FeaturePlot(Combined_Seurat, features = c("percent_mito","nFeature_RNA", "Doublet_Score"), ncol=3)
pdf("Plots/Quality_Control/QC_UMAP.pdf", height = 4, width = 12)
  print(QC_UMAP_Plot)
  dev.off()

# Perform clustering -----------------------
# Identifies clusters of cells within the UMAP
# resolutions and dims values were found in Eric's paper
Combined_Seurat <- FindNeighbors(Combined_Seurat, reduction = "pca", dims = 1:25)
Combined_Seurat <- FindClusters(Combined_Seurat, resolution = scConfig.clustering_resolution)

# Get the number of cells per each cluster
cells_per_cluster <- table(Combined_Seurat$seurat_clusters)
write.csv(cells_per_cluster,"CSV_Results/Cells_per_cluster.csv")

# Save the clustered Seurat object
saveRDS(Combined_Seurat, paste0("R_Data/",scConfig.Prefix ,"_combined_clustered.rds"))

# Examine the resulting UMAP-------------
Clustered_UMAP_Plot <- DimPlot(Combined_Seurat, label = TRUE)
pdf("Plots/Clustering_Plots/Clustered_UMAP.pdf", height = 4, width = 6)
  print(Clustered_UMAP_Plot)
  dev.off()

# Visualize QC metrics in each cluster
percent_mito_vln_Plot <- VlnPlot(Combined_Seurat, features = "percent_mito", pt.size = 0)
pdf("Plots/Clustering_Plots/PercentMito_Violin.pdf")
  print(percent_mito_vln_Plot)
  dev.off()

nFeature_vln_Plot <- VlnPlot(Combined_Seurat, features = "nFeature_RNA", pt.size = 0)
pdf("Plots/Clustering_Plots/overall_nFeature_Violin.pdf")
  print(nFeature_vln_Plot)
  dev.off()

nCount_vln_Plot <- VlnPlot(Combined_Seurat, features = "nCount_RNA", pt.size = 0)
pdf("Plots/Clustering_Plots/nCount_Violin.pdf")
  print(nCount_vln_Plot)
  dev.off()

Doublet_score_vln_Plot <- VlnPlot(Combined_Seurat, features = "Doublet_Score", pt.size = 0)
pdf("Plots/Quality_Control/Doublet_Score_Violin.pdf")
  print(Doublet_score_vln_Plot)
  dev.off()

# Visualize marker gene expression in each cluster
# Using markers from my paper
id_features <- c("Mbp", "Mobp", "Plp1", "Gad1", "Gad2",
                 "Ndrg2", "Slc1a2", "Slc4a4",
                 "Slc17a7", "Satb1", "Neurod6","Vcan", 
                 "Pdgfra", "Pcdh15", "Csf1r" , "Apbb1ip", "P2ry12", 
                 "Flt1", "B2m", "Bmp4", "Cnp", "Ccdc153", 
                 "Rsph1","Tmem212", "Rbfox3")

Major_cells_dotplot <- DotPlot(Combined_Seurat, features = id_features)+RotatedAxis()
pdf("Plots/Clustering_Plots/Major_cells_dotplot.pdf")
  print(Major_cells_dotplot)
  dev.off()

# Examine QC metrics by animal
Idents(Combined_Seurat) <- Combined_Seurat$Sample_name

nFeature_ViolinPlot_byAnimal <- VlnPlot(Combined_Seurat, features = "nFeature_RNA", pt.size = 0)
pdf("Plots/Quality_Control/nFeature_ViolinPlot_byAnimal.pdf", height = 4, width = 6)
  print(nFeature_ViolinPlot_byAnimal)
  dev.off()

nCount_ViolinPlot_byAnimal <- VlnPlot(Combined_Seurat, features = "nCount_RNA", pt.size = 0)
pdf("Plots/Quality_Control/nCount_ViolinPlot_byAnimal.pdf", height = 4, width = 6)
  print(nCount_ViolinPlot_byAnimal)
  dev.off()

percent_mito_ViolinPlot_byAnimal <- VlnPlot(Combined_Seurat, features = "percent_mito", pt.size = 0)
pdf("Plots/Quality_Control/pctMito_ViolinPlot_byAnimal.pdf", height = 4, width = 6)
  print(percent_mito_ViolinPlot_byAnimal)
  dev.off()

# Log the completion time
write(paste0("Cluster_and_ID_Cells - Finish: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
