# Experiment configuration file for Messing/Mayfield scRNA-seq analysis using Seurat
# Written by Geoffrey A. Dilly in June 2024

# Project_name - Datatype: String
# Name of the project for Seurat
scConfig$Project_name <- "GD_scRNA_Experiment"

# Project_name - Datatype: String
# Prefix for filenames of R objects
scConfig$Prefix <- "GD_scRNA"

# Home_folder - Datatype: String (Path)
# Home directory where all analysis directories will be located
# This is generally not needed but good to note
scConfig$Home_folder <-   "__HOME_DIR__"

# Raw_data_folder - Datatype: String (Path)
# Location of the 10x cellranger runs that will be used as raw data
# Should contain data that can be read by the Read10X() function in Seurat
scConfig$Raw_data_folder <- "Raw_Data"

# mito_pattern - Datatype: String
# Str pattern that identifies mitochondrial genes
# Mouse = "^mt-"; Rat = "^Mt-"; Human = "^MT-"
scConfig$mito_pattern <- "^mt-"

# ribo_pattern - Datatype: String
# Str pattern that identifies ribosomal genes
scConfig$ribo_pattern <- "^Rp[ls]"

# nFeature_RNA_cutoff - Datatype: Int 
# Minimum number of nCount detected to include cell
# Recommended = ≥200
scConfig$nFeature_RNA_cutoff <- 200

# percent_mito_cutoff - Datatype: Int
# Maximum percentage of mitochondrial reads to include cell
# Recommended = 5
scConfig$percent_mito_cutoff <- 5

# percent_ribo_cutoff - Datatype: Int
# Maximum percentage of ribosomal reads to include cell
# Recommended = 5
scConfig$percent_ribo_cutoff <- 5

# expected_doublet_pct - Datatype: Int
# Expected percentage of doublets in the dataset
# Recommended = 7.5 (5 to 10%)
scConfig$expct_doublet_pct <- 7.5

# compute_soupx - Datatype: Bool (TRUE/FALSE)
# If TRUE, soupX will be run on the raw data
scConfig$compute_soupx <- FALSE

# soupx_adjust - Datatype: Bool (TRUE/FALSE)
# If TRUE, use SoupX to adjust the raw counts for ambient RNA contamination
scConfig$soupx_adjust <- FALSE

# remove_doublets - Datatype: Bool (TRUE/FALSE)
# If TRUE doublets identified by DoubletFinder be removed prior to SCT normalization
scConfig$remove_doublets <- FALSE

# remove_mito_genes - Datatype: Bool (TRUE/FALSE)
# If TRUE mitochondrial genes will be removed prior to clustering
scConfig$remove_mito_genes <- FALSE

# remove_ribo_genes - Datatype: Bool (TRUE/FALSE)
# If TRUE ribosomal genes will be removed prior to clustering
scConfig$remove_ribo_genes <- FALSE

# remove_top_nUMIs - Datatype: Bool (TRUE/FALSE)
# If TRUE the 25% of cells with the highest UMI count will be removed from the data
scConfig$remove_top_nUMIs <- FALSE

# clustering_PCAs - Datatype: Int
# Default is 25 but a better estimate can be obtained by running ElbowPlot()
scConfig$clustering_PCAs <- 25

# clustering_resolution - Datatype: Int
# Resolution may be estimated with Clustree
scConfig$clustering_resolution <- 0.5

# cluster_plot_ident - Datatype: String
# Identity used by cluster and DGE plots
scConfig$cluster_plot_ident <- "seurat_clusters"

# exogenous_genes - Datatype: List of Strings
# List of exogenous genes to be excluded from clustering
scConfig$exogenous_genes <- c("eGFP", "Cre")

# marker_gene_reference - Datatype: String
# reference from marker_genes_db.csv to use for marker gene plots
scConfig$marker_gene_reference <- "Dilly_et_al_2022"