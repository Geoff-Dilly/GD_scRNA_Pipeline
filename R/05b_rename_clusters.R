# 05b_rename_clusters.R
# Purpose: Insert cell type labels into a Seurat object
# Author: Geoff Dilly

library(here)
library(yaml)
library(Seurat)

# Setup ####
# Load custom functions
source(here::here("R/modules/log_utils.R"))
source(here::here("R/modules/qc_utils.R"))
source(here::here("R/modules/plot_utils"))

# Load the configuration file and metadata
scConfig <- yaml::read_yaml(here::here("sc_experiment_config.yaml"))
scConfig$Sample_metadata <- read.csv(here::here("sc_sample_metadata.csv"))

# Check for required directories
check_required_dirs()

# Log the start time and a timestamped copy of the script
write(paste0("05b_rename_clusters - Start: ", Sys.time()), file = here::here("scRNA_Log.txt"), append = TRUE)
log_file <- write_script_log(here::here("R/05b_rename_clusters.R"))

# Log all output to the end of the log file
sink(log_file, append = TRUE)
sink(log_file, type = "message", append = TRUE)
on.exit({
  sink(NULL)
  sink(NULL, type = "message")
})

# Load data ####
# Load clustered Seurat object
combined_seurat <- readRDS(here::here("R_Data", paste0(scConfig$prefix, "_combined_clustered.rds")))

# Rename Seurat clusters ####
# Note: Manually edit this field to the correct length
combined_seurat <- RenameIdents(object = combined_seurat,
                                "0" = "Cluster_0",
                                "1" = "Cluster_1",
                                "2" = "Cluster_2",
                                "3" = "Cluster_3",
                                "4" = "Cluster_4",
                                "5" = "Cluster_5",
                                "6" = "Cluster_6",
                                "7" = "Cluster_7",
                                "8" = "Cluster_8",
                                "9" = "Cluster_9",
                                "10" = "Cluster_10",
                                "11" = "Cluster_11",
                                "12" = "Cluster_12",
                                "13" = "Cluster_13",
                                "14" = "Cluster_14",
                                "15" = "Cluster_15",
                                "16" = "Cluster_16",
                                "17" = "Cluster_17",
                                "18" = "Cluster_18",
                                "19" = "Cluster_19",
                                "20" = "Cluster_20",
                                "21" = "Cluster_21",
                                "22" = "Cluster_22",
                                "23" = "Cluster_23",
                                "24" = "Cluster_24",
                                "25" = "Cluster_25",
                                "26" = "Cluster_26",
                                "27" = "Cluster_27",
                                "28" = "Cluster_28",
                                "29" = "Cluster_29",
                                "30" = "Cluster_30",
                                "31" = "Cluster_31",
                                "32" = "Cluster_32",
                                "33" = "Cluster_33",
                                "34" = "Cluster_34",
                                "35" = "Cluster_35",
                                "36" = "Cluster_36",
                                "37" = "Cluster_37",
                                "38" = "Cluster_38",
                                "39" = "Cluster_39",
                                "40" = "Cluster_40")

# Save the new names to combined_seurat$CellType
combined_seurat$CellType <- Idents(combined_seurat)

# Reset the Idents to seurat_clusters
Idents(combined_seurat) <- combined_seurat$seurat_clusters

# Save the clustered Seurat object
saveRDS(combined_seurat, here::here("R_Data", paste0(scConfig$prefix, "_combined_clustered.rds")))

# Log the completion time
write(paste0("05b_rename_clusters - Finish: ", Sys.time()), file = here::here("scRNA_Log.txt"), append = TRUE)
