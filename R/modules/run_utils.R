#' @title Check Required Directories
#' @description This function checks if the required directories exist and suggests install.sh if not.
check_required_dirs <- function() {
  required_dirs <- c("csv_results/Cluster_Counts", "csv_results/DEGs_All",
                     "csv_results/Marker_Genes_All", "Raw_Data", "Plots/Cluster_Plots/Marker_Feature_Plots/",
                     "Plots/Cluster_Plots/Marker_Violin_Plots/", "Plots/DGE_Plots/Dispersion_Plots",
                     "Plots/DGE_Plots/Heatmaps", "Plots/DGE_Plots/MA_Plots", "Plots/DGE_Plots/PCAs",
                     "Plots/DGE_Plots/Volcano_Plots", "Plots/Quality_Control", "R_Data", "Logs")

  required_dirs <- c("results", "data/R_Data", "data/Raw_Data", "R")

  missing_dirs <- required_dirs[!dir.exists(required_dirs)]

  if (length(missing_dirs) > 0) {
    message("The following required directories are missing:\n",
            paste(missing_dirs, collapse = "\n"))
    stop("Required directories are missing. Please run install.sh to set up the directory structure.")
  }
}

#' @title Make Results Directory
#' @description This function creates a results directory structure for the scRNA analysis pipeline.
get_results_dir <- function(run_time = Sys.getenv("RUN_TIME"), 
                           prefix = scConfig$prefix) {

  # Create a run name based on the project name and run time
  run_name <- paste0(prefix, "_", run_time)
  run_out_dir <- file.path("results", run_name)

# Base folders
  dir.create(here::here(run_out_dir), showWarnings = FALSE, recursive = TRUE)
  dir.create(here::here(run_out_dir, "csv_results"), showWarnings = FALSE)
  dir.create(here::here(run_out_dir, "plots"), showWarnings = FALSE)
  dir.create(here::here(run_out_dir, "logs"), showWarnings = FALSE)

  # QC and plot subfolders
  dir.create(here::here(run_out_dir, "plots", "Quality_Control"), showWarnings = FALSE)
  dir.create(here::here(run_out_dir, "plots", "Cluster_Plots", "Marker_Feature_Plots"), showWarnings = FALSE, recursive = TRUE)
  dir.create(here::here(run_out_dir, "plots", "Cluster_Plots", "Marker_Violin_Plots"), showWarnings = FALSE, recursive = TRUE)
  dir.create(here::here(run_out_dir, "plots", "DGE_Plots", "Dispersion_Plots"), showWarnings = FALSE, recursive = TRUE)
  dir.create(here::here(run_out_dir, "plots", "DGE_Plots", "Heatmaps"), showWarnings = FALSE, recursive = TRUE)
  dir.create(here::here(run_out_dir, "plots", "DGE_Plots", "MA_Plots"), showWarnings = FALSE, recursive = TRUE)
  dir.create(here::here(run_out_dir, "plots", "DGE_Plots", "PCAs"), showWarnings = FALSE, recursive = TRUE)
  dir.create(here::here(run_out_dir, "plots", "DGE_Plots", "Volcano_Plots"), showWarnings = FALSE, recursive = TRUE)
  dir.create(here::here(run_out_dir, "csv_results", "Marker_Genes_All"), showWarnings = FALSE)
  dir.create(here::here(run_out_dir, "csv_results", "Cluster_Counts"), showWarnings = FALSE)
  dir.create(here::here(run_out_dir, "csv_results", "DEGs_All"), showWarnings = FALSE)


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