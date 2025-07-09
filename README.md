# Single-cell analysis pipeline 

## Author: Geoff Dilly

## Overview:

This repository provides a streamlined pipeline for analyzing 10x Chromium single-cell data using R and Seurat. The pipeline performs QC, clusters cells, identifies marker genes, and generates plots and CSV files in a recommended directory structure. This pipeline is intended for research in the Messing Lab at UT Austin, and can be adapted for similar use. 

## Contents:

- Install.sh
	Bash script to generate the recommended file structure in a local directory. 

- R Scripts
	- Proceessing 10x data for analysis with Seurat
	- Cluster and marker gene identification  
	- Data visualization

- Configuration and Metadata
	- R configuration file to set pipeline parameters
	- Metadata CSV file for identifying and labeling samples

- Run_Analysis_Verbose.sh
	Editable Bash script to run the R scripts in order

- Directory Structure
	Includes recommended folders for results, raw data, plots, and logs

## Usage:
	1. Clone the repository

	'''bash
	git clone https://github.com/Geoff-Dilly/GD_scRNA_Pipeline
	cd <repo directory>
	'''

	2. Install file structure
	'''bash
	bash install.sh
	'''

	3. Setup metadata and configuration
	Modify the configuration file and metadata CSV to match your data

	4. Run the analysis pipeline
	'''bash
	bash Run_Analysis_Verbose.sh
	'''

	5. Examine outputs
	View plots, results, and logs

## Acknowledgements:
This repository is maintained for current and future experiments in the Messing Lab at UT Austin. 