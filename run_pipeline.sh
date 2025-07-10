#!/bin/bash

# Find the absolute path to the project directory and call it
SCRIPTS_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PROJECT_HOME="$SCRIPTS_DIR/.." 
cd "$PROJECT_HOME"

# Run the R scripts in order
# Completed or optional scripts can commented out

R --quiet --no-restore --file="R/01_load_data.R";
#R --quiet --no-restore --file="R/01b_soupx.R";
R --quiet --no-restore --file="R/02_doubletfinder.R";
R --quiet --no-restore --file="R/03_normalize_and_integrate.R";
R --quiet --no-restore --file="R/04_cluster_cells.R"; 
R --quiet --no-restore --file="R/05_id_marker_genes.R";
#R --quiet --no-restore --file="R/05b_rename_clusters.R";
R --quiet --no-restore --file="R/06_make_plots.R";
#R --quiet --no-restore --file="R/07_dge_1var.R";
