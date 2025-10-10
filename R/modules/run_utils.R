# log_utils.R

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

#' @title Get Matrix Path
#' @description Extract the matrix from .zip/.tar if necessary
prepare_matrix_dir <- function(raw_data_dir, req_matrix) {

  matrix_dir <- file.path(raw_data_dir, req_matrix)
  tar_file   <- paste0(matrix_dir, ".tar.gz")
  zip_file   <- paste0(matrix_dir, ".zip")

  if (dir.exists(matrix_dir)) {
    message("Matrix dir detected: No extraction needed")

  } else if (file.exists(tar_file)) {
    message("Extracting matrix from .tar file: ", tar_file)
    utils::untar(tar_file, exdir = raw_data_dir)

  } else if (file.exists(zip_file)) {
    message("Extracting matrix from .zip file: ", zip_file)
    utils::unzip(zip_file, exdir = raw_data_dir)

  } else {
    stop(paste("No compatible 10x matrix found in: ", raw_data_dir))
  }
}