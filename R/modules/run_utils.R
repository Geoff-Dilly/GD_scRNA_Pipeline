#' @title Check Required Directories
#' @description This function checks if the required directories exist and suggests install.sh if not.
check_required_dirs <- function() {
  required_dirs <- c("CSV_Results/Cluster_Counts", "CSV_Results/DEGs_All",
                     "CSV_Results/Marker_Genes_All", "Raw_Data", "Plots/Clustering_Plots/Marker_Feature_Plots/",
                     "Plots/Clustering_Plots/Marker_Violin_Plots/", "Plots/DESEQ_Plots/Dispersion_Plots",
                     "Plots/DESEQ_Plots/Heatmaps", "Plots/DESEQ_Plots/MA_Plots", "Plots/DESEQ_Plots/PCAs",
                     "Plots/DESEQ_Plots/Volcano_Plots", "Plots/Quality_Control", "R_Data", "Logs")

  required_dirs <- c("Results", "Logs", "Data/R_Data", "Data/Raw_Data", "R")

  missing_dirs <- required_dirs[!dir.exists(required_dirs)]

  if (length(missing_dirs) > 0) {
    message("The following required directories are missing:\n",
            paste(missing_dirs, collapse = "\n"))
    stop("Required directories are missing. Please run install.sh to set up the directory structure.")
  }
}

#' @title Make Results Directory
#' @description This function creates a results directory structure for the scRNA analysis pipeline.
Get_results_dir <- function(run_time = Sys.getenv("RUN_TIME"), 
                           prefix = scConfig$prefix) {

  # Create a run name based on the project name and run time
  run_name <- paste0(prefix, "_", run_time)
  run_out_dir <- file.path("Results", run_name)

# Base folders
  dir.create(here::here(run_out_dir), showWarnings = FALSE, recursive = TRUE)
  dir.create(here::here(run_out_dir, "CSV_Results"), showWarnings = FALSE)
  dir.create(here::here(run_out_dir, "Plots"), showWarnings = FALSE)
  dir.create(here::here(run_out_dir, "Logs"), showWarnings = FALSE)

  # QC and plot subfolders
  dir.create(here::here(run_out_dir, "Plots", "Quality_Control"), showWarnings = FALSE)
  dir.create(here::here(run_out_dir, "Plots", "Clustering_Plots", "Marker_Feature_Plots"), showWarnings = FALSE, recursive = TRUE)
  dir.create(here::here(run_out_dir, "Plots", "Clustering_Plots", "Marker_Violin_Plots"), showWarnings = FALSE, recursive = TRUE)
  dir.create(here::here(run_out_dir, "Plots", "DESEQ_Plots", "Dispersion_Plots"), showWarnings = FALSE, recursive = TRUE)
  dir.create(here::here(run_out_dir, "Plots", "DESEQ_Plots", "Heatmaps"), showWarnings = FALSE, recursive = TRUE)
  dir.create(here::here(run_out_dir, "Plots", "DESEQ_Plots", "MA_Plots"), showWarnings = FALSE, recursive = TRUE)
  dir.create(here::here(run_out_dir, "Plots", "DESEQ_Plots", "PCAs"), showWarnings = FALSE, recursive = TRUE)
  dir.create(here::here(run_out_dir, "Plots", "DESEQ_Plots", "Volcano_Plots"), showWarnings = FALSE, recursive = TRUE)

  return(run_out_dir)
}

#' @title Write Script Log
#' @description This function writes a log of the script execution,
#' @description including the script content, configuration file, and session info.
write_script_log <- function(script_file, config_file = "sc_experiment_config.yaml", log_dir = "Logs/") {
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

  # Open a log file in the Logs directory
  log_file <- file.path(log_dir, paste0("Time_", timestamp, "_Script_Log_",
                                        tools::file_path_sans_ext(basename(script_file)), ".txt"))

  log_connection <- file(log_file, open = "wt")

  # Write the script to the log file
  writeLines("Script Log:\n", log_connection)
  writeLines(readLines(script_file), log_connection)
  writeLines("\n\n====================================================================\n\n", log_connection)

  # Write the configuration file and session info
  writeLines("Config Log:\n", log_connection)
  writeLines(readLines(config_file), log_connection)
  writeLines("\n\n====================================================================\n\n", log_connection)

  writeLines("Session info:\n", log_connection)
  capture.output(sessionInfo(), file = log_connection)
  writeLines("\n\n====================================================================\n\n", log_connection)

  return(log_connection)
}

#' @title Check Required Directories
#' @description This function checks if the required directories exist and suggests install.sh if not.
check_required_dirs <- function() {
  required_dirs <- c("CSV_Results/Cluster_Counts", "CSV_Results/DEGs_All",
                     "CSV_Results/Marker_Genes_All", "Raw_Data", "Plots/Clustering_Plots/Marker_Feature_Plots/",
                     "Plots/Clustering_Plots/Marker_Violin_Plots/", "Plots/DESEQ_Plots/Dispersion_Plots",
                     "Plots/DESEQ_Plots/Heatmaps", "Plots/DESEQ_Plots/MA_Plots", "Plots/DESEQ_Plots/PCAs",
                     "Plots/DESEQ_Plots/Volcano_Plots", "Plots/Quality_Control", "R_Data", "Logs")

  missing_dirs <- required_dirs[!dir.exists(required_dirs)]

  if (length(missing_dirs) > 0) {
    message("The following required directories are missing:\n",
            paste(missing_dirs, collapse = "\n"))
    stop("Required directories are missing. Please run install.sh to set up the directory structure.")
  }
}