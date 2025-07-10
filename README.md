# Single-cell Analysis Pipeline 

**Author:** Geoff Dilly, PhD  
**Lab:** Messing Lab, UT Austin

## Overview

This repository provides a streamlined pipeline for analyzing 10x Chromium single-cell data using R and Seurat. The pipeline performs QC, clusters cells, identifies marker genes, and generates plots and CSV files in a recommended directory structure. This pipeline is intended for research in the Messing Lab at UT Austin, and can be adapted for similar use. 

## Contents

- `install.sh`: Bash script to generate the recommended file structure in a local directory. 
- `make_environment.sh`: Bash script to make an Anaconda environment with the appropriate dependencies.
- `sc_analysis_env.yaml`: YAML containing Anaconda dependencies.
- `R/`: R scripts for data processing, clustering, cell identification, visualization, and DGE.
- `run_pipeline.sh`: Editable Bash script to run the R scripts in order.
- `sc_experiment_config.R`: R configuration file to set pipeline parameters.
- `sc_sample_metadata.csv`: Metadata CSV file for identifying and labeling samples.

## Usage

1. **Clone the repository:**
    ```sh
    git clone https://github.com/Geoff-Dilly/GD_scRNA_Pipeline
    cd <repo directory>
    ```
2. **Install file structure:**
    ```sh
    bash install.sh
    ```
3. **Setup metadata and configuration:**
   Edit `sc_experiment_config.R` and `sc_sample_metadata.csv` as needed
4. **Run the analysis pipeline:**
    ```sh
    bash run_pipeline.sh
    ```
5. **Examine outputs:**
   View plots, results, and logs.

## Scripts

### Bash Scripts
| Script | Description |
|---|---|
| install.sh | Generate the recommended file structure | 
| make_environment.sh | Creates an Anaconda environment | 
| run_pipeline.sh | Run the R scripts in order | 

### R Scripts
| Script | Description | Output |
|---|---|---|
| 01_load_data.R | Load sample data, metadata, and basic QC | Sample-level Seurat object |
| 01b_soupx.R | *Optional:* Run SoupX on each sample | Sample-level Seurat object |
| 02_doubletfinder.R | Run DoubletFinder on each sample | Sample-level Seurat object |
| 03_normalize_and_integrate.R | Integrate samples and normalize with scTransform | Experiment-level Seurat objects |
| 04_cluster_cells.R | Perform dimensional reduction and clustering analysis | Experiment-level Seurat object |
| 05_id_marker_genes.R | Identify marker genes by cluster | Experiment-level Seurat object |
| 05b_rename_clusters.R | *Optional:* Add cell type names to metadata | Experiment-level Seurat object |
| 06_make_plots.R | Makes various plots | PDF plots |
| 07_dge_1var.R | Differential gene expression analysis (1 variable: Treatment) | CSVs and PDF plots|

## Directory Structure

- `Raw_Data/` — Raw input files (10x count matrices)
- `R_Data/` — RDS files of Seurat objects 
- `Plots/` — Output figures and QC plots
- `CSV_Results/` - CSVs of cell counts, markers, and DGE results
- `Logs/` — Run logs and script backups
- `R/` — R scripts


## Metadata Format

- Metadata is a CSV that must contain four required columns: Sample_name, Sex, Treatment, Raw_data_dir
- Raw_data_dir should direct to the cellranger outs/ folder containing 
- Additional columns are automatically read into the Seurat object as sample-level metadata

**Example:**

| Sample_name | Sex | Treatment | Age | Raw_data_dir |
|---|---|---|---|---|
| Subject_1 | M | Drug | 24 | "Raw_Data/Subject_1/outs" |
| Subject_2 | F | Ctrl | 27 | "Raw_Data/Subject_2/outs" |
| Subject_3 | M | Ctrl | 27 | "Raw_Data/Subject_3/outs" |
| Subject_4 | F | Drug | 26 | "Raw_Data/Subject_4/outs" |

Saved as a CSV: 
```
Sample_name,Sex,Treatment,Age
Subject_1,M,Drug,24,"Raw_Data/Subject_1/outs"
Subject_2,F,Ctrl,27,"Raw_Data/Subject_2/outs"
Subject_3,M,Ctrl,27,"Raw_Data/Subject_3/outs"
Subject_4,F,Drug,26,"Raw_Data/Subject_4/outs"
```

## Requirements

- R = 4.3.3
- Seurat ≥ 5.0
- Anaconda
- See `make_environment.sh` or `sc_analysis_env.yaml` for environment setup.

## Acknowledgements

This repository is maintained for experiments in the Messing Lab at UT Austin. 