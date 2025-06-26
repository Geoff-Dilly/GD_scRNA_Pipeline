# Identify_Marker_Genes.R
# Purpose: Performs FindAllMarkers on a clustered Seurat objects and saves lists of marker genes
# Author: Geoff Dilly

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
snRNA_home_dir <- "__HOME_DIR__"
setwd(snRNA_home_dir)

# Log the start time and a timestamped copy of the script
write(paste0("Identify_Marker_Genes - Start: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
file.copy("Scripts/Identify_Marker_Genes.R", paste0("Logs/Time_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", "Identify_Marker_Genes.R"), overwrite = FALSE)

# Load the configuration file and metadata
source("sc_experiment_config.R")
scConfig.Sample_metadata <- read.csv("sc_sample_metadata.csv")

# It is reccommended that you install the Presto package before running this script
#> install.packages("devtools")
#> devtools::install_github("immunogenomics/presto")


# Load the previously clustered Seurat object
Combined_Seurat <- readRDS(paste0("R_Data/",scConfig.Prefix ,"_combined_clustered.rds"))

# Identify markers for each cluster
Combined_Seurat <- PrepSCTFindMarkers(Combined_Seurat)
All_Markers <- FindAllMarkers(object = Combined_Seurat)
All_Markers$PCT_Fold <- All_Markers$pct.1/All_Markers$pct.2
All_Markers$PCT_Delta <- All_Markers$pct.1-All_Markers$pct.2
write.csv(All_Markers, "CSV_Results/Marker_Genes_All/All_marker_genes.csv")

# Make more manageable lists of the top markers and save as CSV
All_Markers <- read.csv("CSV_Results/Marker_Genes_All/All_marker_genes.csv")

Top30_Cell_type_markers <- All_Markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.csv(Top30_Cell_type_markers,"CSV_Results/Marker_Genes_All/Marker_Genes_Top30.csv")

Top10_Cell_type_markers <- All_Markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(Top10_Cell_type_markers,"CSV_Results/Marker_Genes_All/Marker_Genes_Top10.csv")

# Visualize the top marker gene in each cluster
top_markers <- All_Markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
TopMarker_DotPlot <- DotPlot(Combined_Seurat, features = unique(top_markers$gene))+RotatedAxis()
pdf("Plots/Clustering_Plots/Top_Marker_DotPlot.pdf", height = 8, width = 12)
  print(TopMarker_DotPlot)
  dev.off()

# Visualize the top 2 marker genes in each cluster
top_markers2 <- All_Markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
Top2Markers_DotPlot <- DotPlot(Combined_Seurat, features = unique(top_markers2$gene))+RotatedAxis()
pdf("Plots/Clustering_Plots/Top2_Markers_DotPlot.pdf", height = 8, width = 12)
  print(Top2Markers_DotPlot)
  dev.off()

# Log the completion time
write(paste0("Identify_Marker_Genes - Finish: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
