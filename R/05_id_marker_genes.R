# 05_id_marker_genes.R
# Purpose: Performs FindAllMarkers on a clustered Seurat objects and saves lists of marker genes
# Author: Geoff Dilly

library(here)
library(yaml)
library(dplyr)
library(Seurat)
library(ggplot2)

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

# Log the start time and a timestamped copy of the script
write(paste0("05_id_marker_genes - Start: ", Sys.time()), file = here::here("scRNA_Log.txt"), append = TRUE)
log_connection <- write_script_log(here::here("R/05_id_marker_genes.R"), log_dir = here::here(output_dir, "logs"))

# Log all output to the end of the log file
sink(log_connection, append = TRUE)
sink(log_connection, type = "message", append = TRUE)
on.exit({
  sink(NULL)
  sink(NULL, type = "message")
})

# It is recommended that you install the Presto package before running this script
#> install.packages("devtools")
#> devtools::install_github("immunogenomics/presto")

# Load data ####
# Load the previously clustered Seurat object
combined_seurat <- readRDS(here::here("data", "R_Data", paste0(scConfig$prefix, "_combined_clustered.rds")))

# Find marker genes and save results ####
# Identify markers for each cluster
combined_seurat <- PrepSCTFindMarkers(combined_seurat)
all_markers <- FindAllMarkers(object = combined_seurat, random.seed = 42)
all_markers$PCT_Fold <- all_markers$pct.1 / all_markers$pct.2
all_markers$PCT_Delta <- all_markers$pct.1 - all_markers$pct.2
write.csv(all_markers, here::here(output_dir, "csv_results", "Marker_Genes_All", "All_marker_genes.csv"))
write.csv(all_markers, here::here("reference", "All_marker_genes.csv"))

# Make more manageable lists of the top markers and save as CSV
top30_cell_type_markers <- all_markers %>% group_by(cluster) %>% slice_max(n = 30, order_by = p_val_adj)
write.csv(top30_cell_type_markers, here::here(output_dir, "csv_results", "Marker_Genes_All", "Marker_Genes_Top30.csv"))

top10_cell_type_markers <- all_markers %>% group_by(cluster) %>% slice_max(n = 10, order_by = p_val_adj)
write.csv(top10_cell_type_markers, here::here(output_dir, "csv_results", "Marker_Genes_All", "Marker_Genes_Top10.csv"))

# Make plots ####
# Visualize the top marker gene in each cluster
top_markers <- all_markers %>% group_by(cluster) %>% slice_max(n = 1, order_by = avg_log2FC)
topmarker_dotplot <- DotPlot(combined_seurat, features = unique(top_markers$gene)) + RotatedAxis()
save_plot_pdf(topmarker_dotplot, here::here(output_dir, "plots", "Cluster_Plots", "Top_Markers_DotPlot.pdf"), height = 8, width = 12)

# Visualize the top 2 marker genes in each cluster
top_markers2 <- all_markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
top2markers_dotplot <- DotPlot(combined_seurat, features = unique(top_markers2$gene)) + RotatedAxis()
save_plot_pdf(top2markers_dotplot, here::here(output_dir, "plots", "Cluster_Plots", "Top2_Markers_DotPlot.pdf"), height = 8, width = 12)

# Log the completion time
write(paste0("05_id_marker_genes - Finish: ", Sys.time()), file = here::here("scRNA_Log.txt"), append = TRUE)
