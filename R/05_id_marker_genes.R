# 05_id_marker_genes.R
# Purpose: Performs FindAllMarkers on a clustered Seurat objects and saves lists of marker genes
# Author: Geoff Dilly

library(here)
library(dplyr)
library(Seurat)
library(ggplot2)
snRNA_home_dir <- here()
setwd(snRNA_home_dir)

# Log the start time and a timestamped copy of the script
write(paste0("05_id_marker_genes - Start: ", Sys.time()), file = "snRNA_Log.txt", append = TRUE)
file.copy("Scripts/05_id_marker_genes.R", paste0("Logs/Time_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", "05_id_marker_genes.R"), overwrite = FALSE)

# Load the configuration file and metadata
source("sc_experiment_config.R")
scConfig.Sample_metadata <- read.csv("sc_sample_metadata.csv")

# It is reccommended that you install the Presto package before running this script
#> install.packages("devtools")
#> devtools::install_github("immunogenomics/presto")

# Load the previously clustered Seurat object
combined_seurat <- readRDS(paste0("R_Data/", scConfig.Prefix, "_combined_clustered.rds"))

# Identify markers for each cluster
combined_seurat <- PrepSCTFindMarkers(combined_seurat)
all_markers <- FindAllMarkers(object = combined_seurat)
all_markers$PCT_Fold <- all_markers$pct.1 / all_markers$pct.2
all_markers$PCT_Delta <- all_markers$pct.1 - all_markers$pct.2
write.csv(all_markers, "CSV_Results/Marker_Genes_All/All_marker_genes.csv")

# Make more manageable lists of the top markers and save as CSV
all_markers <- read.csv("CSV_Results/Marker_Genes_All/All_marker_genes.csv")

top30_cell_type_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.csv(top30_cell_type_markers, "CSV_Results/Marker_Genes_All/Marker_Genes_Top30.csv")

top10_cell_type_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(top10_cell_type_markers, "CSV_Results/Marker_Genes_All/Marker_Genes_Top10.csv")

# Visualize the top marker gene in each cluster
top_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
topmarker_dotplot <- DotPlot(combined_seurat, features = unique(top_markers$gene)) + RotatedAxis()
pdf("Plots/Clustering_Plots/Top_Marker_DotPlot.pdf", height = 8, width = 12)
print(topmarker_dotplot)
dev.off()

# Visualize the top 2 marker genes in each cluster
top_markers2 <- all_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2markers_dotplot <- DotPlot(combined_seurat, features = unique(top_markers2$gene)) + RotatedAxis()
pdf("Plots/Clustering_Plots/Top2_Markers_DotPlot.pdf", height = 8, width = 12)
print(top2markers_dotplot)
dev.off()

# Log the completion time
write(paste0("05_id_marker_genes - Finish: ", Sys.time()), file = "snRNA_Log.txt", append = TRUE)
