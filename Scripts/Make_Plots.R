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
