#' @title Write Script Log
#' @description This function writes a log of the script execution,
#' @description including the script content, configuration file, and session info.
write_script_log <- function(script_file) {
  config_file <- "sc_experiment_config.R"
  timestamp <- format(Sys.time(), "%Y-%m-%d_%H-%M-%S")

  # Open a log file in the Logs directory
  log_file <- paste0("Logs/Time_", timestamp, "_Script_Log_",
                     tools::file_path_sans_ext(basename(script_file)), ".txt")

  log_connection <- file(log_file, open = "wt")

  # Write the script to the log file
  writeLines(paste("Script Log:\n"), log_connection)
  writeLines(readLines(script_file), log_connection)

  # Write the configuration file and session info
  writeLines(paste("\nConfig Log:\n"), log_connection)
  writeLines(readLines(config_file), log_connection)

  writeLines("\nSession info:\n", log_connection)
  capture.output(sessionInfo(), file = log_connection)

  # Close the log file connection
  close(log_connection)
}

write_script_log("R/01_load_data.R")