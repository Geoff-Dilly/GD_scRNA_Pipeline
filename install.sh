#!/bin/bash

# Get file contaning sn_directoy_setup script
SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
cd $SCRIPT_DIR

# Create the directories that the single-cell pipeline needs to run
mkdir -p "CSV_Results/Cluster_Counts"
mkdir -p "CSV_Results/DEGs_All"
mkdir -p "CSV_Results/Marker_Genes_All"
mkdir -p "Raw_Data"
mkdir -p "Plots/Clustering_Plots/Marker_Feature_Plots/"
mkdir -p "Plots/Clustering_Plots/Marker_Violin_Plots/"
mkdir -p "Plots/DESEQ_Plots/dispersion_plots"
mkdir -p "Plots/DESEQ_Plots/Heatmaps"
mkdir -p "Plots/DESEQ_Plots/MA_Plots"
mkdir -p "Plots/DESEQ_Plots/PCAs"
mkdir -p "Plots/DESEQ_Plots/Volcano_Plots"
mkdir -p "Plots/Quality_Control"
mkdir -p "R_Data"
mkdir -p "Scripts"
mkdir -p "Logs"

snRNA_home_dir=$(pwd)
echo $snRNA_home_dir

echo 'snRNA_home_dir <- "'${snRNA_home_dir}'"' >> "analysis_home_dir.R"
cat "analysis_home_dir.R"

set_home_line=$'snRNA_home_dir <- "'${snRNA_home_dir}'"'

cat > Scripts/Load_10x_Data.R <<'EOF'
# Load_10x_Data.R
# Purpose: Load a list of 10x objects into Seurat and save them as RDS
# Author: Geoff Dilly

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(DropletUtils)
snRNA_home_dir <- "__HOME_DIR__"
setwd(snRNA_home_dir)

# Log the start time and a timestamped copy of the script
write(paste0("Load_10x_Data - Start: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
file.copy("Scripts/Load_10x_Data.R", paste0("Logs/Time_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", "Load_10x_Data.R"), overwrite = FALSE)

# Load the configuration file and metadata
source("sc_experiment_config.R")
scConfig.Sample_metadata <- read.csv("sc_sample_metadata.csv")

# Check for required metadata columns
if(!("Sample_name" %in% colnames(scConfig.Sample_metadata))) {
  stop("Mandatory metadata column <Sample_name> is not present")
}

if(!("Treatment" %in% colnames(scConfig.Sample_metadata))) {
  stop("Mandatory metadata column <Treatment> is not present")
}

if(!("Sex" %in% colnames(scConfig.Sample_metadata))) {
  stop("Mandatory metadata column <Sex> is not present")
}

# Loads the Cell Ranger output into Seurat
for (i in 1:nrow(scConfig.Sample_metadata)) {
sample <- scConfig.Sample_metadata[i, ]
sample_seurat.data <- Read10X(paste0(scConfig.Raw_data_folder, "/", sample$Sample_name, "/", sample$Sample_name, "/outs/filtered_feature_bc_matrix"))
sample_seurat <- CreateSeuratObject(counts = sample_seurat.data, project = scConfig.Project_name, min.cells = 1, min.features = 1)
sample_seurat  <- PercentageFeatureSet(sample_seurat, pattern = scConfig.mito_pattern, col.name = "percent_mito")
sample_seurat  <- PercentageFeatureSet(sample_seurat, pattern = scConfig.ribo_pattern, col.name = "percent_ribo")
sample_seurat <- subset(sample_seurat, subset = nFeature_RNA > scConfig.nFeature_RNA_cutoff & percent_mito < scConfig.percent_mito_cutoff)

# Loads each column of the metadata into the seurat object 
for(col in colnames(sample)) {
  sample_seurat[[col]] <- sample[[col]]
}

# Save the Seurat object
saveRDS(sample_seurat, file = paste0("R_Data/", sample$Sample_name, "_seurat.rds"))
}

# Log the completion time
write(paste0("Load_10x_Data - Finish: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
EOF

sed 's|__HOME_DIR__|'"$snRNA_home_dir"'|' Scripts/Load_10x_Data.R > temp && mv temp Scripts/Load_10x_Data.R

cat > Scripts/Merge_and_Normalize.R <<'EOF'
# Merge_and_Normalize.R
# Purpose: Merges Seurat objects from multiple samples and normalized the resulting object with SCTransform
# Author: Geoff Dilly

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
snRNA_home_dir <- "__HOME_DIR__"
setwd(snRNA_home_dir)

# Log the start time and a timestamped copy of the script
write(paste0("Merge_and_Normalize - Start: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
file.copy("Scripts/Merge_and_Normalize.R", paste0("Logs/Time_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", "Merge_and_Normalize.R"), overwrite = FALSE)

# Set 'R_MAX_VSIZE' to maximum RAM usage
Sys.setenv('R_MAX_VSIZE'=32000000000)

# Load the configuration file and metadata
source("sc_experiment_config.R")
scConfig.Sample_metadata <- read.csv("sc_sample_metadata.csv")

# Get sample names in a list of strings
str_sample_list <- scConfig.Sample_metadata$Sample_name

# Initialize a list to store the Seurat objects
seurat_objects <- list()

# Read RDS files and assign them to variables dynamically
for (sample in str_sample_list) {
  sample_seurat <- readRDS(paste0("R_Data/", sample, "_seurat_Doublets.rds"))
  seurat_objects[[sample]] <- sample_seurat
}

# Merge samples into one Seurat object without integration
# Convert the list of Seurat objects to a list of arguments for the merge function
#! MergeNorm_Combined_Seurat <- merge, c(seurat_objects[[1]], seurat_objects[-1]))
MergeNorm_Combined_Seurat <- merge(x = seurat_objects[[1]], y = seurat_objects[-1], add.cell.ids = str_sample_list)

saveRDS(MergeNorm_Combined_Seurat, paste0("R_Data/",scConfig.Prefix ,"_combined_merged.rds"))
rm(seurat_objects)
MergeNorm_Combined_Seurat <- readRDS(paste0("R_Data/",scConfig.Prefix ,"_combined_merged.rds"))

# Remove called doublets if you want to
if(scConfig.remove_doublets == TRUE) {
  MergeNorm_Combined_Seurat <- subset(MergeNorm_Combined_Seurat, subset = Doublet_Call == "Singlet")}

# Remove mitochondrial genes if you want to
if(scConfig.remove_mito_genes == TRUE) {
  MergeNorm_Combined_Seurat <- JoinLayers(MergeNorm_Combined_Seurat)
  mito.genes <- grep(pattern = scConfig.mito_pattern, x = rownames(x = MergeNorm_Combined_Seurat@assays$RNA$counts), value = TRUE)
  counts <- GetAssayData(MergeNorm_Combined_Seurat, assay = "RNA", layer = "counts")
  counts <- counts[-(which(rownames(counts) %in% mito.genes)),]
  MergeNorm_Combined_Seurat <- subset(MergeNorm_Combined_Seurat, features = rownames(counts))
  MergeNorm_Combined_Seurat[["RNA"]] <- split(MergeNorm_Combined_Seurat[["RNA"]], f = MergeNorm_Combined_Seurat$Sample_name)
  Layers(MergeNorm_Combined_Seurat)
  rm(counts)
}

# Run SCTransform and save the normalized Seurat object
MergeNorm_Combined_Seurat <- SCTransform(MergeNorm_Combined_Seurat, verbose = TRUE, conserve.memory=TRUE)
saveRDS(MergeNorm_Combined_Seurat, paste0("R_Data/",scConfig.Prefix ,"_combined_SCT.rds"))

# Log the completion time
write(paste0("Merge_and_Normalize - Finish: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
EOF

sed 's|__HOME_DIR__|'"$snRNA_home_dir"'|' Scripts/Merge_and_Normalize.R > temp && mv temp Scripts/Merge_and_Normalize.R

cat > Scripts/Cluster_and_ID_Cells.R <<'EOF'
# Cluster_and_ID_Cells.R
# Purpose: Performs the standard Seurat clustering workflow on a normalized Seurat object
# Author: Geoff Dilly

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
snRNA_home_dir <- "__HOME_DIR__"
setwd(snRNA_home_dir)

# Log the start time and a timestamped copy of the script
write(paste0("Cluster_and_ID_Cells - Start: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
file.copy("Scripts/Cluster_and_ID_Cells.R", paste0("Logs/Time_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", "Cluster_and_ID_Cells.R"), overwrite = FALSE)

# Load the configuration file and metadata
source("sc_experiment_config.R")
scConfig.Sample_metadata <- read.csv("sc_sample_metadata.csv")

# Load normalized data
Combined_Seurat <- readRDS(paste0("R_Data/",scConfig.Prefix ,"_combined_SCT.rds"))

# Get the number of cells per each sample
cells_per_sample <- table(Combined_Seurat$Sample_name)
write.csv(cells_per_sample,"CSV_Results/Cells_per_sample.csv")

# Check quality metrics for each cell
QC_Violins_Plot <- VlnPlot(Combined_Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "Doublet_Score"), ncol = 4, pt.size = 0)
pdf("Plots/Quality_Control/QC_VlnPlot.pdf", height = 6, width = 10)
  print(QC_Violins_Plot)
  dev.off()

# Scale data and run UMAP ----------------------
# Perform PCA
Combined_Seurat <- RunPCA(Combined_Seurat, npcs = 100 ,verbose = TRUE)

# Visualize the dimensionality of the PCs and pick the number of PCs
DimHeatmaps_Plot <- DimHeatmap(Combined_Seurat, dims = 1:9, cells = 500, balanced = TRUE)
pdf("Plots/Quality_Control/Dim_Heatmap.pdf", height = 6, width = 8)
  print(DimHeatmaps_Plot)
  dev.off()

Elbow_Plot <- ElbowPlot(Combined_Seurat, ndims = 100)
pdf("Plots/Quality_Control/Elbowplot.pdf")
  print(Elbow_Plot)
  dev.off()
# I chose 25 PCs

# Perform UMAP dimensional reduction on the data
Combined_Seurat <- RunUMAP(Combined_Seurat, reduction = "pca", dims = 1:25)

# Examine the resulting UMAP-------------
Raw_UMAP_Plot <- DimPlot(Combined_Seurat)
pdf("Plots/Clustering_Plots/Raw_UMAP.pdf", height = 4, width = 6)
  print(Raw_UMAP_Plot)
  dev.off()

# Examine the UMAP Plot for quality control and viability
QC_UMAP_Plot <- FeaturePlot(Combined_Seurat, features = c("percent_mito","nFeature_RNA", "Doublet_Score"), ncol=3)
pdf("Plots/Quality_Control/QC_UMAP.pdf", height = 4, width = 12)
  print(QC_UMAP_Plot)
  dev.off()

# Perform clustering -----------------------
# Identifies clusters of cells within the UMAP
# resolutions and dims values were found in Eric's paper
Combined_Seurat <- FindNeighbors(Combined_Seurat, reduction = "pca", dims = 1:25)
Combined_Seurat <- FindClusters(Combined_Seurat, resolution = scConfig.clustering_resolution)

# Get the number of cells per each cluster
cells_per_cluster <- table(Combined_Seurat$seurat_clusters)
write.csv(cells_per_cluster,"CSV_Results/Cells_per_cluster.csv")

# Save the clustered Seurat object
saveRDS(Combined_Seurat, paste0("R_Data/",scConfig.Prefix ,"_combined_clustered.rds"))

# Examine the resulting UMAP-------------
Clustered_UMAP_Plot <- DimPlot(Combined_Seurat, label = TRUE)
pdf("Plots/Clustering_Plots/Clustered_UMAP.pdf", height = 4, width = 6)
  print(Clustered_UMAP_Plot)
  dev.off()

# Visualize QC metrics in each cluster
percent_mito_vln_Plot <- VlnPlot(Combined_Seurat, features = "percent_mito", pt.size = 0)
pdf("Plots/Clustering_Plots/PercentMito_Violin.pdf")
  print(percent_mito_vln_Plot)
  dev.off()

nFeature_vln_Plot <- VlnPlot(Combined_Seurat, features = "nFeature_RNA", pt.size = 0)
pdf("Plots/Clustering_Plots/overall_nFeature_Violin.pdf")
  print(nFeature_vln_Plot)
  dev.off()

nCount_vln_Plot <- VlnPlot(Combined_Seurat, features = "nCount_RNA", pt.size = 0)
pdf("Plots/Clustering_Plots/nCount_Violin.pdf")
  print(nCount_vln_Plot)
  dev.off()

Doublet_score_vln_Plot <- VlnPlot(Combined_Seurat, features = "Doublet_Score", pt.size = 0)
pdf("Plots/Quality_Control/Doublet_Score_Violin.pdf")
  print(Doublet_score_vln_Plot)
  dev.off()

# Visualize marker gene expression in each cluster
# Using markers from my paper
id_features <- c("Mbp", "Mobp", "Plp1", "Gad1", "Gad2",
                 "Ndrg2", "Slc1a2", "Slc4a4",
                 "Slc17a7", "Satb1", "Neurod6","Vcan", 
                 "Pdgfra", "Pcdh15", "Csf1r" , "Apbb1ip", "P2ry12", 
                 "Flt1", "B2m", "Bmp4", "Cnp", "Ccdc153", 
                 "Rsph1","Tmem212", "Rbfox3")

Major_cells_dotplot <- DotPlot(Combined_Seurat, features = id_features)+RotatedAxis()
pdf("Plots/Clustering_Plots/Major_cells_dotplot.pdf")
  print(Major_cells_dotplot)
  dev.off()

# Examine QC metrics by animal
Idents(Combined_Seurat) <- Combined_Seurat$Sample_name

nFeature_ViolinPlot_byAnimal <- VlnPlot(Combined_Seurat, features = "nFeature_RNA", pt.size = 0)
pdf("Plots/Quality_Control/nFeature_ViolinPlot_byAnimal.pdf", height = 4, width = 6)
  print(nFeature_ViolinPlot_byAnimal)
  dev.off()

nCount_ViolinPlot_byAnimal <- VlnPlot(Combined_Seurat, features = "nCount_RNA", pt.size = 0)
pdf("Plots/Quality_Control/nCount_ViolinPlot_byAnimal.pdf", height = 4, width = 6)
  print(nCount_ViolinPlot_byAnimal)
  dev.off()

percent_mito_ViolinPlot_byAnimal <- VlnPlot(Combined_Seurat, features = "percent_mito", pt.size = 0)
pdf("Plots/Quality_Control/pctMito_ViolinPlot_byAnimal.pdf", height = 4, width = 6)
  print(percent_mito_ViolinPlot_byAnimal)
  dev.off()

# Log the completion time
write(paste0("Cluster_and_ID_Cells - Finish: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
EOF

sed 's|__HOME_DIR__|'"$snRNA_home_dir"'|' Scripts/Cluster_and_ID_Cells.R > temp && mv temp Scripts/Cluster_and_ID_Cells.R

cat > Scripts/Identify_Marker_Genes.R <<'EOF'
# Identify_Marker_Genes.R
# Purpose: Performs FindAllMarkers on a clustered Seurat objects and saves lists of marker genes
# Author: Geoff Dilly

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
snRNA_home_dir <- "__HOME_DIR__"
setwd(snRNA_home_dir)

# Log the start time and a timestamped copy of the script
write(paste0("Identify_Marker_Genes - Start: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
file.copy("Scripts/Identify_Marker_Genes.R", paste0("Logs/Time_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", "Identify_Marker_Genes.R"), overwrite = FALSE)

# Load the configuration file and metadata
source("sc_experiment_config.R")
scConfig.Sample_metadata <- read.csv("sc_sample_metadata.csv")

# It is reccommended that you install the Presto package before running this script
#> install.packages("devtools")
#> devtools::install_github("immunogenomics/presto")


# Load the previously clustered Seurat object
Combined_Seurat <- readRDS(paste0("R_Data/",scConfig.Prefix ,"_combined_clustered.rds"))

# Identify markers for each cluster
Combined_Seurat <- PrepSCTFindMarkers(Combined_Seurat)
All_Markers <- FindAllMarkers(object = Combined_Seurat)
All_Markers$PCT_Fold <- All_Markers$pct.1/All_Markers$pct.2
All_Markers$PCT_Delta <- All_Markers$pct.1-All_Markers$pct.2
write.csv(All_Markers, "CSV_Results/Marker_Genes_All/All_marker_genes.csv")

# Make more manageable lists of the top markers and save as CSV
All_Markers <- read.csv("CSV_Results/Marker_Genes_All/All_marker_genes.csv")

Top30_Cell_type_markers <- All_Markers %>% group_by(cluster) %>% top_n(n = 30, wt = avg_log2FC)
write.csv(Top30_Cell_type_markers,"CSV_Results/Marker_Genes_All/Marker_Genes_Top30.csv")

Top10_Cell_type_markers <- All_Markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(Top10_Cell_type_markers,"CSV_Results/Marker_Genes_All/Marker_Genes_Top10.csv")

# Visualize the top marker gene in each cluster
top_markers <- All_Markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
TopMarker_DotPlot <- DotPlot(Combined_Seurat, features = unique(top_markers$gene))+RotatedAxis()
pdf("Plots/Clustering_Plots/Top_Marker_DotPlot.pdf", height = 8, width = 12)
  print(TopMarker_DotPlot)
  dev.off()

# Visualize the top 2 marker genes in each cluster
top_markers2 <- All_Markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
Top2Markers_DotPlot <- DotPlot(Combined_Seurat, features = unique(top_markers2$gene))+RotatedAxis()
pdf("Plots/Clustering_Plots/Top2_Markers_DotPlot.pdf", height = 8, width = 12)
  print(Top2Markers_DotPlot)
  dev.off()

# Log the completion time
write(paste0("Identify_Marker_Genes - Finish: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
EOF

sed 's|__HOME_DIR__|'"$snRNA_home_dir"'|' Scripts/Identify_Marker_Genes.R > temp && mv temp Scripts/Identify_Marker_Genes.R

cat > Scripts/Rename_Clusters.R <<'EOF'
# Rename_Clusters.R
# Purpose: Insert cell type labels into a Seurat object
# Author: Geoff Dilly

library(Seurat)
snRNA_home_dir <- "__HOME_DIR__"
setwd(snRNA_home_dir)

# Log the start time and a timestamped copy of the script
write(paste0("Rename_Clusters - Start: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
file.copy("Scripts/Rename_Clusters.R", paste0("Logs/Time_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", "Rename_Clusters.R"), overwrite = FALSE)

# Load the configuration file and metadata
source("sc_experiment_config.R")
scConfig.Sample_metadata <- read.csv("sc_sample_metadata.csv")

# Load clustered Seurat object
Combined_Seurat <- readRDS(paste0("R_Data/",scConfig.Prefix ,"_combined_clustered.rds"))

# Rename each Seurat cluster
# Note: You will have to manually edit this field to the correct length
# Running this script will substantially increase the size of the Seurat object
Combined_Seurat <- RenameIdents(object = Combined_Seurat,
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

# Save the new names to Combined_Seurat$CellType
Combined_Seurat$CellType <- Idents(Combined_Seurat)

# Reset the Idents to seurat_clusters
Idents(Combined_Seurat) <- Combined_Seurat$seurat_clusters

# Save the clustered Seurat object
saveRDS(Combined_Seurat, paste0("R_Data/",scConfig.Prefix ,"_combined_clustered.rds"))

# Log the completion time
write(paste0("Rename_Clusters - Finish: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
EOF

sed 's|__HOME_DIR__|'"$snRNA_home_dir"'|' Scripts/Rename_Clusters.R > temp && mv temp Scripts/Rename_Clusters.R

cat > Scripts/Make_Plots.R <<'EOF'
# Make_Plots.R
# Purpose: Make various plots for single-cell analysis
# Author: Geoff Dilly

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
snRNA_home_dir <- "__HOME_DIR__"
setwd(snRNA_home_dir)

# Log the start time and a timestamped copy of the script
write(paste0("Make_Plots - Start: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
file.copy("Scripts/Make_Plots.R", paste0("Logs/Time_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", "Make_Plots.R"), overwrite = FALSE)

# Load the configuration file
source("SC_Experiment_config.R")

# Load the clustered Seurat object
Combined_Seurat <- readRDS(paste0("R_Data/",scConfig.Prefix ,"_combined_clustered.rds"))

# Make all of the figures from previous scripts --------------

# Check quality metrics for each cell
QC_Violins_Plot <- VlnPlot(Combined_Seurat, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3, pt.size = 0)
pdf("Plots/Quality_Control/QC_VlnPlot.pdf")
  print(QC_Violins_Plot)
  dev.off()

# Examine the raw UMAP
Raw_UMAP_Plot <- DimPlot(Combined_Seurat)
pdf("Plots/Clustering_Plots/Raw_UMAP.pdf", height = 4, width = 6)
  print(Raw_UMAP_Plot)
  dev.off()

# Examine QC metrics on the UMAP plot
QC_UMAP_Plot <- FeaturePlot(Combined_Seurat, features = c("percent_mito","nFeature_RNA", "Doublet_Score"))
pdf("Plots/Quality_Control/QC_UMAP.pdf", height = 4, width = 10)
  print(QC_UMAP_Plot)
  dev.off()

# Examine percent ribosomal RNA on a UMAP
Ribo_UMAP_Plot <- FeaturePlot(Combined_Seurat, features = "percent_ribo")
pdf("Plots/Quality_Control/PctRibo_UMAP.pdf", height = 4, width = 10)
  print(Ribo_UMAP_Plot)
  dev.off()

# Examine doublet score in doublets and non doublets
Doublet_Call_UMAP <- FeaturePlot(Combined_Seurat, features = "Doublet_Score", split.by = "Doublet_Call")
pdf("Plots/Quality_Control/Doublet_Call_UMAP.pdf", height = 4, width = 10)
  print(Doublet_Call_UMAP)
  dev.off()

# Examine QC metrics by animal 
Idents(Combined_Seurat) <- Combined_Seurat$Sample_name

# Umap
Sample_UMAP_Plot <- DimPlot(Combined_Seurat)
pdf("Plots/Quality_Control/Sample_UMAP.pdf", height = 4, width = 6)
  print(Sample_UMAP_Plot)
  dev.off()

# Violin plots
nFeature_ViolinPlot_byAnimal <- VlnPlot(Combined_Seurat, features = "nFeature_RNA", pt.size = 0)
pdf("Plots/Quality_Control/nFeature_ViolinPlot_byAnimal.pdf", height = 4, width = 6)
  print(nFeature_ViolinPlot_byAnimal)
  dev.off()

nCount_ViolinPlot_byAnimal <- VlnPlot(Combined_Seurat, features = "nCount_RNA", pt.size = 0)
pdf("Plots/Quality_Control/nCount_ViolinPlot_byAnimal.pdf", height = 4, width = 6)
  print(nCount_ViolinPlot_byAnimal)
  dev.off()

percent_mito_ViolinPlot_byAnimal <- VlnPlot(Combined_Seurat, features = "percent_mito", pt.size = 0)
pdf("Plots/Quality_Control/pctMito_ViolinPlot_byAnimal.pdf", height = 4, width = 6)
  print(percent_mito_ViolinPlot_byAnimal)
  dev.off()

percent_ribo_ViolinPlot_byAnimal <- VlnPlot(Combined_Seurat, features = "percent_ribo", pt.size = 0)
pdf("Plots/Quality_Control/pctribo_ViolinPlot_byAnimal.pdf", height = 4, width = 6)
  print(percent_ribo_ViolinPlot_byAnimal)
  dev.off()

# Examine proportions
#plot samples as proportion or percentage of cluster
Sample_prop_barPlot <- ggplot(Combined_Seurat@meta.data, aes(x=seurat_clusters, fill=Sample_name)) + 
  geom_bar(position = "fill") +
  RotatedAxis() + 
  xlab("Sample Name")
pdf("Plots/Quality_Control/Sample_prop_barPlot.pdf", height = 4, width = 6)
  print(Sample_prop_barPlot)
  dev.off()

#Plot samples as proportion or percentage of cluster
Cluster_prop_barPlot <- ggplot(Combined_Seurat@meta.data, aes(x=Sample_name, fill=seurat_clusters)) + 
  geom_bar(position = "fill") +
  RotatedAxis() + 
  xlab("Sample Name")
pdf("Plots/Quality_Control/Cluster_prop_barPlot.pdf", height = 4, width = 6)
  print(Cluster_prop_barPlot)
  dev.off()

#Plot proportions of each condition in each cluster
#!Cond_prop_barPlot <- ggplot(Combined_Seurat@meta.data, aes(x=seurat_clusters, fill=Virus)) + 
#!  geom_bar(position = "fill") + 
#!  geom_hline(yintercept = 12/33, linetype = "dashed", size = 1) + 
#!  geom_hline(yintercept = 22/33, linetype = "dashed", size = 1) +
#!  RotatedAxis() + 
#!  xlab("Sample Name")
#!pdf("Plots/Quality_Control/Cond_prop_barPlot.pdf", height = 4, width = 6)
#!  print(Cond_prop_barPlot)
#!  dev.off()

#plot proportions of each condition in each cluster
Cond_prop_barPlot <- ggplot(Combined_Seurat@meta.data, aes(x=seurat_clusters, fill=Treatment)) + 
  geom_bar(position = "fill") + 
  geom_hline(yintercept = 12/33, linetype = "dashed", size = 1) + 
  geom_hline(yintercept = 22/33, linetype = "dashed", size = 1) +
  RotatedAxis() + 
  xlab("Sample Name")
pdf("Plots/Quality_Control/Cond_prop_barPlot.pdf", height = 4, width = 6)
  print(Cond_prop_barPlot)
  dev.off()

#plot proportions of each sex in each cluster
Sex_prop_barPlot <- ggplot(Combined_Seurat@meta.data, aes(x=seurat_clusters, fill=Sex)) + 
  geom_bar(position = "fill") + 
  geom_hline(yintercept = 0.5, linetype = "dashed", size = 1) +
  RotatedAxis() + 
  xlab("Cluster")
pdf("Plots/Quality_Control/Sex_prop_barPlot.pdf", height = 4, width = 6)
  print(Sex_prop_barPlot)
  dev.off()

# Clustered QC plots ------------------

# Examine the clustered UMAP
Idents(Combined_Seurat) <- Combined_Seurat$seurat_clusters

Clustered_UMAP_Plot <- DimPlot(Combined_Seurat)
pdf("Plots/Quality_Control/Clustered_UMAP.pdf", height = 4, width = 6)
  print(Clustered_UMAP_Plot)
  dev.off()

# Visualize QC metrics in each cluster
percent_mito_vln_Plot <- VlnPlot(Combined_Seurat, features = "percent_mito", pt.size = 0)
pdf("Plots/Quality_Control/PercentMito_Violin.pdf")
  print(percent_mito_vln_Plot)
  dev.off()

nFeature_vln_Plot <- VlnPlot(Combined_Seurat, features = "nFeature_RNA", pt.size = 0)
pdf("Plots/Quality_Control/nFeature_Violin.pdf", height = 4, width = 8)
  print(nFeature_vln_Plot)
  dev.off()

nCount_vln_Plot <- VlnPlot(Combined_Seurat, features = "nCount_RNA", pt.size = 0)
pdf("Plots/Quality_Control/nCount_Violin.pdf", height = 4, width = 8)
  print(nCount_vln_Plot)
  dev.off()

nCount_vln_Plot <- VlnPlot(Combined_Seurat, features = "nCount_RNA", pt.size = 0)
pdf("Plots/Quality_Control/nCount_Violin.pdf", height = 4, width = 8)
  print(nCount_vln_Plot)
  dev.off()

Doublet_score_vln_Plot <- VlnPlot(Combined_Seurat, features = "Doublet_Score", pt.size = 0)
pdf("Plots/Quality_Control/Doublet_Score_Violin.pdf", height = 4, width = 8)
  print(Doublet_score_vln_Plot)
  dev.off()

Percent_ribo_vln_Plot <- VlnPlot(Combined_Seurat, features = "percent_ribo", pt.size = 0)
pdf("Plots/Quality_Control/percent_ribo_Violin.pdf", height = 4, width = 8)
  print(Percent_ribo_vln_Plot)
  dev.off()

# Look at the marker genes identified by Seurat -------------
All_Markers <- read.csv("CSV_results/Marker_Genes_All/All_marker_genes.csv")

# Visualize the top marker gene expression in each cluster
top_markers <- All_Markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
TopMarker_DotPlot <- DotPlot(Combined_Seurat, features = unique(top_markers$gene))+RotatedAxis()
pdf("Plots/Clustering_Plots/Top_Marker_DotPlot.pdf", height = 8, width = 12)
  print(TopMarker_DotPlot)
  dev.off()

# Visualize the top 2 marker genes expression in each cluster
top_markers2 <- All_Markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
Top2Markers_DotPlot <- DotPlot(Combined_Seurat, features = unique(top_markers2$gene))+RotatedAxis()
pdf("Plots/Clustering_Plots/Top2_Markers_DotPlot.pdf", height = 8, width = 12)
  print(Top2Markers_DotPlot)
  dev.off()

# Visualize the top 2 marker genes expression in each cluster
top_markers2_diff <- All_Markers %>% group_by(cluster) %>% top_n(n = 2, wt = PCT_Delta)
Top2Markers_diff_DotPlot <- DotPlot(Combined_Seurat, features = unique(top_markers2_diff$gene))+RotatedAxis()
pdf("Plots/Clustering_Plots/Top2_Markers_diff_DotPlot.pdf", height = 8, width = 12)
print(Top2Markers_diff_DotPlot)
dev.off()
  
  
# Marker gene plots ----------------
# Using markers from Dilly et al 2022 ---------------
id_features <- c("Mbp", "Mobp", "Plp1", "Gad1", "Gad2",
                 "Ndrg2", "Slc1a2", "Slc4a4",
                 "Slc17a7", "Satb1", "Neurod6","Vcan", 
                 "Pdgfra", "Pcdh15", "Csf1r" , "Apbb1ip", "P2ry12", 
                 "Flt1", "B2m", "Bmp4", "Cnp", "Ccdc153", 
                 "Rsph1","Tmem212", "Rbfox3")

Major_cells_dotplot <- DotPlot(Combined_Seurat, features = id_features)+RotatedAxis()
pdf("Plots/Clustering_Plots/Major_cells_dotplot.pdf")
  print(Major_cells_dotplot)
  dev.off()

# Relevant marker genes
# Major Cell Types ---------------------

# Neurons
All_Neurons <- c("Rbfox3", "Syn1", "Gria1")
Inhibitory_Neurons <- c("Gad1", "Gad2")
Excitatory_Neurons <- c("Neurod6", "Satb1", "Slc17a7") # Can be divided into 2 classes
Cholinergic_Neurons <- c("Slc5a7", "Prima1")

# Glia
Astrocytes <- c("Slc1a2", "Slc1a3", "Ndrg2")
Oligodendrocytes <- c("Mog", "Mobp", "Hapln2")
Microglia <- c("Csf1r","Apbb1ip", "Inpp5d")
OPCs <- c("Vcan", "Pdgfra", "Pcdh15")
Immature_Oligos <- c("Bmp4", "Cnp", "Tcf7l2") # - These might benefit from a more conservative name

# Other Cells
Endothelial_Cells <- c("Flt1", "B2m", "Rgs5")
Ependymal_Cells <- c("Ccdc153", "Rsph1", "Tmem212")

# Lists for programming
DS_cell_types <- c("Astrocytes", "Oligodendrocytes", "Microglia", "OPCs", 
                      "Immature_Oligos", "All_Neurons", "Inhibitory_Neurons", "Excitatory_Neurons",
                      "Cholinergic_Neurons", "Endothelial_Cells", "Ependymal_Cells")

Major_marker_genes <- c(Astrocytes, Oligodendrocytes, Microglia, OPCs, 
                        Immature_Oligos, All_Neurons, Inhibitory_Neurons, Excitatory_Neurons,
                        Cholinergic_Neurons, Endothelial_Cells, Ependymal_Cells)

Major_marker_genes_unique <- unique(Major_marker_genes)

Major_marker_genes_list <- list(Astrocytes, Oligodendrocytes, Microglia, OPCs, 
                                Immature_Oligos, All_Neurons, Inhibitory_Neurons, Excitatory_Neurons,
                                Cholinergic_Neurons, Endothelial_Cells, Ependymal_Cells)


# DS Cell Subtypes ----------------------
# From
dSPN <- c("Drd1", "Ebf1", "Pdyn") #Pdyn over Ebf1?
iSPN <- c("Drd2", "Adora2a", "Penk")
eSPN <- c("Otof", "Col11a1")
PVALB_IN <- c("Kit", "Pvalb")
Chol_IN <- c("Tacr1", "Chat")
SST_IN <- c("Sst", "Npy")
Astro <- c("Slc1a3", "Rorb")
Oligo <- c("Mog", "Aspa")
OPC <- c("Pdgfra")
Endo <- c("Flt1", "Slco1a4")
Ependy <- c("Dnah12", "Rsph1", "Gm973")
Mural <- c("Pdgfrb", "Rgs5")
Microglia <- c("C1qc", "Cx3cr1")

DS_cell_types <- c(dSPN, iSPN, eSPN, PVALB_IN, Chol_IN,
                   SST_IN, Astro, Oligo, OPC, Endo, 
                   Ependy, Mural, Microglia)


DS_marker_genes_list <- list(dSPN, iSPN, eSPN, PVALB_IN, Chol_IN,
                   SST_IN, Astro, Oligo, OPC, Endo, 
                   Ependy, Mural, Microglia)

DS_marker_genes_unique <- unique(DS_marker_genes_list)

# Examine the resulting UMAP-------------
Clustered_UMAP_Plot <- DimPlot(Combined_Seurat, label = TRUE)
pdf("Plots/Clustering_Plots/Clustered_UMAP.pdf", height = 6, width = 8)
  print(Clustered_UMAP_Plot)
  dev.off()

all_markers_list <- c(Major_marker_genes_unique, DS_marker_genes_unique)

# Generate featureplots and violinplots for genes of interest
for (marker in all_markers_list){
  gene_feature_Plot <- FeaturePlot(Combined_Seurat, features = marker)
  pdf(paste0("Plots/Clustering_Plots/Marker_Feature_Plots/", marker, "_FeaturePlot.pdf"), height = 4, width = 6)
  print(gene_feature_Plot)
  dev.off()

  gene_violin_Plot <- VlnPlot(Combined_Seurat, features = marker, pt.size = 0)
  pdf(paste0("Plots/Clustering_Plots/Marker_Violin_Plots/", marker, "_ViolinPlot.pdf"), height = 4, width = 6)
  print(gene_violin_Plot)
  dev.off()
}

# Log the completion time
write(paste0("Make_Plots - Finish: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
EOF

sed 's|__HOME_DIR__|'"$snRNA_home_dir"'|' Scripts/Make_Plots.R > temp && mv temp Scripts/Make_Plots.R

cat > Scripts/Run_SoupX.R <<'EOF'
# Run_SoupX.R
# Purpose: Run SoupX on single sample Seurat objects
# Author: Geoff Dilly

library(Seurat)
library(SoupX)
library(DropletUtils)
library(ggplot2)
library(DoubletFinder)
library(knitr)
snRNA_home_dir <- "__HOME_DIR__"
setwd(snRNA_home_dir)

# Log the start time and a timestamped copy of the script
write(paste0("Run_SoupX - Start: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
file.copy("Scripts/Run_SoupX.R", paste0("Logs/Time_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", "Run_SoupX.R"), overwrite = FALSE)

# Code adapted from https://cellgeni.github.io/notebooks/html/new-10kPBMC-SoupX.html

# Load the configuration file and metadata
source("sc_experiment_config.R", local = TRUE)
scConfig.Sample_metadata <- read.csv("sc_sample_metadata.csv")

# Place each sample in a list for further processing 
str_sample_list <- scConfig.Sample_metadata$Sample_name


# Process each sample
for (sample_name in str_sample_list) {
  filt.matrix <- Read10X(paste0(scConfig.Raw_data_folder, "/", sample$Sample_name, "/", sample$Sample_name, "/outs/filtered_feature_bc_matrix"))
  raw.matrix <- Read10X(paste0(scConfig.Raw_data_folder, "/", sample$Sample_name, "/", sample$Sample_name, "/outs/raw_feature_bc_matrix"))
  str(raw.matrix)
  str(filt.matrix)

  # Create a Seurat object with the filtered data
  GD_10x_Sample <- CreateSeuratObject(counts = filt.matrix)

  # Create the SoupX channel
  soup.channel <- SoupChannel(raw.matrix, filt.matrix)

  # Cluster the Seurat object
  GD_10x_Sample <- SCTransform(GD_10x_Sample, verbose = T)
  GD_10x_Sample <- RunPCA(GD_10x_Sample, verbose = T)
  GD_10x_Sample <- RunUMAP(GD_10x_Sample, dims = 1:30, verbose = T)
  GD_10x_Sample <- FindNeighbors(GD_10x_Sample, dims = 1:30, verbose = T)
  GD_10x_Sample <- FindClusters(GD_10x_Sample, verbose = T)

  # Extract metadata and UMAP coordinates
  meta <- GD_10x_Sample@meta.data
  umap <- GD_10x_Sample@reductions$umap@cell.embeddings

  # Assign clustering to SoupX
  soup.channel <- setClusters(soup.channel, setNames(meta$seurat_clusters, rownames(meta)))
  soup.channel <- setDR(soup.channel, umap)
  head(meta)

  # Estimate and apply the SoupX correction
  soup.channel <- autoEstCont(soup.channel)
  adj.matrix <- adjustCounts(soup.channel, roundToInt = T)
  
  # Display the top ambient RNA
  head(soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ], n = 50)
  
  # Save the corrected matrix
  write10xCounts(paste0("Raw_Data/", sample_name, "_filtered_feature_bc_matrix_SoupX_adjusted.h5"), adj.matrix)

  # Optional visualization for SoupX correction
  #!topsoupgenes <- head(soup.channel$soupProfile[order(soup.channel$soupProfile$est, decreasing = T), ], n = 100)
  #!neatgene <- "Apoe"
  #!soupchangeplot <- plotChangeMap(soup.channel, adj.matrix, neatgene)
  #!soupchangeplot
}

# Log the completion time
write(paste0("Run_SoupX - Finish: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
EOF

sed 's|__HOME_DIR__|'"$snRNA_home_dir"'|' Scripts/Run_SoupX.R > temp && mv temp Scripts/Run_SoupX.R

cat > Scripts/Run_DoubletFinder.R <<'EOF'
# Run_DoubletFinder.R
# Purpose: Run DoubletFinder on single sample Seurat objects
# Author: Geoff Dilly

library(Seurat)
library(stringr)
library(DoubletFinder)
library(data.table)
snRNA_home_dir <- "__HOME_DIR__"
setwd(snRNA_home_dir)

# Log the start time and a timestamped copy of the script
write(paste0("Run_DoubletFinder - Start: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
file.copy("Scripts/Run_DoubletFinder.R", paste0("Logs/Time_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", "Run_DoubletFinder.R.R"), overwrite = FALSE)

# Read the sample metadata file
scConfig.Sample_metadata <- read.csv("sc_sample_metadata.csv")

str_sample_list <- scConfig.Sample_metadata$Sample_name

# Run doubletFinder ---------------------------------------------------------------------------------------
for (sample_name in str_sample_list) {
  # Pre-process Seurat object (standard)
  GD_10x_Sample <- LoadSeuratRds(paste0("R_Data/", sample_name, "_seurat.rds"))

  GD_10x_Sample <- NormalizeData(GD_10x_Sample)
  GD_10x_Sample <- FindVariableFeatures(GD_10x_Sample, selection.method = "vst", nfeatures = 2000)
  GD_10x_Sample <- ScaleData(GD_10x_Sample)
  GD_10x_Sample <- RunPCA(GD_10x_Sample)
  GD_10x_Sample <- RunUMAP(GD_10x_Sample, dims = 1:10)

  # pK Identification (no ground-truth)
  sweep.res.list <- paramSweep(GD_10x_Sample, PCs = 1:10, sct = FALSE)
  sweep.stats <- summarizeSweep(sweep.res.list, GT = FALSE)
  bcmvn <- find.pK(sweep.stats)

  pK_value <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric), drop=T]))
  print(pK_value)

  # Homotypic Doublet Proportion Estimate
  homotypic.prop <- modelHomotypic(GD_10x_Sample@meta.data$Sample_name)           ## ex: annotations <- seu_kidney@meta.data$ClusteringResults
  nExp_poi <- round(0.075*nrow(GD_10x_Sample@meta.data))  ## Assuming 7.5% doublet formation rate - tailor for your dataset
  nExp_poi.adj <- round(nExp_poi*(1-homotypic.prop))

  # Run DoubletFinder
  GD_10x_Sample <- doubletFinder(GD_10x_Sample, PCs = 1:10, pN = 0.25, pK = pK_value, nExp = nExp_poi, sct = FALSE)
  
  # Rename the doubletFinder results
  meta_cols <- colnames(GD_10x_Sample@meta.data)
  score <- str_subset(meta_cols, "^pANN")
  call <- str_subset(meta_cols, "^DF.cl")
  GD_10x_Sample$Doublet_Score <- GD_10x_Sample[[score]]
  GD_10x_Sample$Doublet_Call <- GD_10x_Sample[[call]]
  GD_10x_Sample[[call]] <- NULL
  GD_10x_Sample[[score]] <- NULL
  
  # Remove unnecessary data layers
  GD_10x_Sample[["RNA"]]$scale.data <- NULL
  GD_10x_Sample[["RNA"]]$data <- NULL
  
  # Save the Seurat object with doubletFinder Results
  saveRDS(GD_10x_Sample, file = paste0("R_Data/", sample_name, "_seurat_Doublets.rds"))

}

# Examine the results
#!CeA_sample <- LoadSeuratRds("/Volumes/users/geoff_scratch/astrocyte_reanalysis/data/r_data/GD_2A_Alcohol_doublets.rds")
#!DimPlot(CeA_sample)
#!FeaturePlot(CeA_sample, features = "pANN_0.25_0.24_709", cols = c("blue", "yellow"), max.cutoff = 0.8)
#!FeaturePlot(CeA_sample, features = "DF.classifications_0.25_0.24_709")
#!doublet_freqs <- as.data.frame(table(CeA_sample@meta.data$DF.classifications_0.25_0.24_709))
#!print(doublet_freqs)
#!percentage_doublets <- doublet_freqs$Freq[1] / (doublet_freqs$Freq[2] + doublet_freqs$Freq[1])*100
#!percentage_doublets

# Log the completion time
write(paste0("Run_DoubletFinder - Finish: ", Sys.time()),file="snRNA_Log.txt", append = TRUE)
EOF

sed 's|__HOME_DIR__|'"$snRNA_home_dir"'|' Scripts/Run_DoubletFinder.R > temp && mv temp Scripts/Run_DoubletFinder.R

cat > sc_sample_metadata.csv <<'EOF'
Sample_name,Sex,Treatment
<SAMPLE>,<M>,<CONTROL>
EOF

cat > sc_experiment_config.R <<'EOF'
# Experiment configuration file for Messing/Mayfield snRNA-seq analysis using Seurat
# Written by Geoffrey A. Dilly in June 2024

# Project_name - Datatype: String
# Name of the project for Seurat
scConfig.Project_name <- "GD_snRNA_Experiment"

# Project_name - Datatype: String
# Prefix for filenames of R objects
scConfig.Prefix <- "GD_snRNA"

# Home_folder - Datatype: String (Path)
# Home directory where all analysis directories will be located
# This is generally not needed but good to note
scConfig.Home_folder <-   "__HOME_DIR__"

# Raw_data_folder - Datatype: String (Path)
# Location of the 10x cellranger runs that will be used as raw data
# Should contain data that can be read by the Read10X() function in Seurat
# Data should be located at <Raw_data_folder>/<Sample_name>/<Sample_name>/outs/filtered_feature_bc_matrix
# In order to run SoupX, .../unfiltered_feature_bc_matrix is also necessary
scConfig.Raw_data_folder <- "Raw_Data"

# mito_pattern - Datatype: String
# Str pattern that identifies mitochondrial genes
# Mouse = "^mt-"; Rat = "^Mt-"; Human = "^MT-"
scConfig.mito_pattern <- "^mt-"

# ribo_pattern - Datatype: String
# Str pattern that identifies mitochondrial genes
scConfig.ribo_pattern <- "^Rp[ls]"

# nFeature_RNA_cutoff - Datatype: Int 
# Minimum number of nCount detected to include cell
# Recommended = ≥200
scConfig.nFeature_RNA_cutoff <- 200

# percent_mito_cutoff - Datatype: Int
# Maximum percentage of mitochondrial reads to include cell
# Recommended = 5
scConfig.percent_mito_cutoff <- 5

# remove_doublets - Datatype: Bool (TRUE/FALSE)
# If TRUE doublets identified by DoubletFinder be removed prior to SCT normalization
scConfig.remove_doublets <- FALSE

# remove_mito_genes - Datatype: Bool (TRUE/FALSE)
# If TRUE mitochondrial genes will be removed prior to clustering
scConfig.remove_mito_genes <- FALSE

# remove_top_nUMIs - Datatype: Bool (TRUE/FALSE)
# If TRUE the 25% of cells with the highest UMI count will be removed from the data
scConfig.remove_top_nUMIs <- FALSE

# clustering_PCAs - Datatype: Int
# Default is 25 but a better estimate can be obtained by running ElbowPlot()
scConfig.clustering_PCAs <- 25

# clustering_resolution - Datatype: Int
# Resolution may be estimated with Clustree
scConfig.clustering_resolution <- 0.5
EOF

sed 's|__HOME_DIR__|'"$snRNA_home_dir"'|' sc_experiment_config.R > temp && mv temp sc_experiment_config.R

cat > Scripts/Run_Analysis_Verbose.bash <<'EOF'
#!/bin/bash

cd "__HOME_DIR__";

R --quiet --no-restore --file="Scripts/Load_10x_Data.R";
#!R --quiet --no-restore --file="Scripts/Run_SoupX.R";
R --quiet --no-restore --file="Scripts/Run_DoubletFinder.R";
R --quiet --no-restore --file="Scripts/Merge_and_Normalize.R";
R --quiet --no-restore --file="Scripts/Cluster_and_ID_Cells.R"; 
R --quiet --no-restore --file="Scripts/Identify_Marker_Genes.R";
#!R --quiet --no-restore --file="Scripts/Rename_Clusters.R";
R --quiet --no-restore --file="Scripts/Make_Plots.R";
#!R --quiet --no-restore --file="Scripts/One_Var_DGE_Analysis.R";
EOF

sed 's|__HOME_DIR__|'"$snRNA_home_dir"'|' Scripts/Run_Analysis_Verbose.bash > temp && mv temp Scripts/Run_Analysis_Verbose.bash

chmod u+x Scripts/Run_Analysis_Verbose.bash