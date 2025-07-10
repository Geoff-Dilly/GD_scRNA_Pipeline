# Single-cell Analysis Pipeline 

## Author: Geoff Dilly

## Overview:

This repository provides a streamlined pipeline for analyzing 10x Chromium single-cell data using R and Seurat. The pipeline performs QC, clusters cells, identifies marker genes, and generates plots and CSV files in a recommended directory structure. This pipeline is intended for research in the Messing Lab at UT Austin, and can be adapted for similar use. 

## Contents:

- install.sh
	Bash script to generate the recommended file structure in a local directory. 

- R Scripts
	- Proceessing 10x data for analysis with Seurat
	- Cluster and marker gene identification  
	- Data visualization

- Configuration and Metadata
	- R configuration file to set pipeline parameters
	- Metadata CSV file for identifying and labeling samples

- run_pipeline.sh
	Editable Bash script to run the R scripts in order

- Directory Structure
	Includes recommended folders for results, raw data, plots, and logs

## Usage:

1. Clone the repository:

    ```shell
    git clone https://github.com/Geoff-Dilly/GD_scRNA_Pipeline
    cd <repo directory>
    ```

2. Install file structure:

    ```shell
    bash install.sh
    ```

3. Setup metadata and configuration:  
   Modify the configuration file and metadata CSV to match your data

4. Run the analysis pipeline:

    ```shell
    bash run_pipeline.sh
    ```

5. Examine outputs:  
   View plots, results, and logs

## Scripts:

### Bash Scripts:
| Script | Description |
|---|---|
| install.sh | Generate the recommended file structure | 
| make_environment.sh | Creates an Anaconda environment | 
| run_pipeline.sh | Run the R scripts in order | 

### R Scripts
| Script | Description | Output |
|---|---|---|
| 01_load_data.R | Load sample data, metadata, and basic QC | Sample-level seurat object |
| 01b_soupx.R | *Optional:* Run SoupX on each sample | Sample-level seurat object *Optional*|
| 02_doubletfinder.R | Run DoubletFinder on each sample | Sample-level seurat object |
| 03_normalize_and_integrate.R | Integrate samples and normalize with scTransform | Experiment-level seurat objects |
| 04_cluster_cells.R | Perform dim reduction and clustering analysis | Experiment-level seurat object |
| 05_id_marker_genes.R | Identify marker genes by cluster | Experiment-level seurat object |
| 05b_rename_clusters.R | *Optional:* Add cell type names to metadata | Experiment-level seurat object |
| 06_make_plots.R | Makes various plots | Plots as PDFs |
| 07_dge_1var.R | Differential gene expression analysis (1 variable) | CSVs and Plots as PDFs |

## Directory Structure

- `Raw_Data` — Raw input files (fastq, count matrices, etc.)
- `R_Data/` — RDS files of Seurat objects at each stage
- `Plots/` — Output figures and QC plots
- `CSV_Results` - CSVs of cell counts, markers, and DGE results
- `Logs/` — Run logs and script backups
- `R/` — R scripts


## Metadata Format:

- Metadata is read from a CSV and included in the seurat object (see example below)

- Three columns are mandatory: Sample_name, Sex, Treatment

- Additionally columns are automatically read into the Seurat object as sample-level metadata

### Example:

As a table:

| Sample_name | Sex | Treatment | Age |
|---|---|---|---|
| Subject_1 | M | Drug | 24 |
| Subject_2 | F | Ctrl | 27 |
| Subject_3 | M | Ctrl | 27 |
| Subject_4 | F | Ctrl | 26 |

As a CSV:
    
```
Sample_name,Sex,Treatment,Age
Subject_1,M,Drug,24
Subject_2,F,Ctrl,27
Subject_3,M,Ctrl,27
Subject_4,F,Ctrl,26
```

## Acknowledgements:

This repository is maintained for current and future experiments in the Messing Lab at UT Austin. 