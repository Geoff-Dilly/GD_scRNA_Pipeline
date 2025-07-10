# 06_06_make_plots.R
# Purpose: Make various plots for single-cell analysis
# Author: Geoff Dilly

library(here)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
snRNA_home_dir <- here()
setwd(snRNA_home_dir)

# Log the start time and a timestamped copy of the script
write(paste0("06_make_plots - Start: ", Sys.time()), file = "snRNA_Log.txt", append = TRUE)
file.copy("R/06_make_plots.R", paste0("Logs/Time_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", "06_make_plots.R"), overwrite = FALSE)

# Load the configuration file
source("sc_experiment_config.R")

# Load the clustered Seurat object
combined_seurat <- readRDS(paste0("R_Data/", scConfig.Prefix, "_combined_clustered.rds"))

# Make all of the figures from previous scripts ####

# Check quality metrics for each cell
qc_vlns_plot <- VlnPlot(combined_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent_mito"), ncol = 3, pt.size = 0)
pdf("Plots/Quality_Control/QC_VlnPlot.pdf")
print(qc_vlns_plot)
dev.off()

# Examine the raw UMAP
raw_umap_plot <- DimPlot(combined_seurat)
pdf("Plots/Clustering_Plots/Raw_UMAP.pdf", height = 4, width = 6)
print(raw_umap_plot)
dev.off()

# Examine QC metrics on the UMAP plot
qc_umap_plot <- FeaturePlot(combined_seurat, features = c("percent_mito", "nFeature_RNA", "Doublet_Score"))
pdf("Plots/Quality_Control/QC_UMAP.pdf", height = 4, width = 10)
print(qc_umap_plot)
dev.off()

# Examine percent ribosomal RNA on a UMAP
ribo_umap_plot <- FeaturePlot(combined_seurat, features = "percent_ribo")
pdf("Plots/Quality_Control/PctRibo_UMAP.pdf", height = 4, width = 10)
print(ribo_umap_plot)
dev.off()

# Examine doublet score in doublets and non doublets
doublet_call_umap_plot <- FeaturePlot(combined_seurat, features = "Doublet_Score", split.by = "Doublet_Call")
pdf("Plots/Quality_Control/Doublet_Call_UMAP.pdf", height = 4, width = 10)
print(doublet_call_umap_plot)
dev.off()

# Examine QC metrics by animal
Idents(combined_seurat) <- combined_seurat$Sample_name

# Umap
sample_umap_plot <- DimPlot(combined_seurat)
pdf("Plots/Quality_Control/Sample_UMAP.pdf", height = 4, width = 6)
print(sample_umap_plot)
dev.off()

# Violin plots
nFeature_vln_by_animal <- VlnPlot(combined_seurat, features = "nFeature_RNA", pt.size = 0)
pdf("Plots/Quality_Control/nFeature_vln_byAnimal.pdf", height = 4, width = 6)
print(nFeature_vln_by_animal)
dev.off()

nCount_vln_by_animal <- VlnPlot(combined_seurat, features = "nCount_RNA", pt.size = 0)
pdf("Plots/Quality_Control/nCount_vln_byAnimal.pdf", height = 4, width = 6)
print(nCount_vln_by_animal)
dev.off()

percent_mito_vln_by_animal <- VlnPlot(combined_seurat, features = "percent_mito", pt.size = 0)
pdf("Plots/Quality_Control/pctMito_vln_byAnimal.pdf", height = 4, width = 6)
print(percent_mito_vln_by_animal)
dev.off()

percent_ribo_vln_by_animal <- VlnPlot(combined_seurat, features = "percent_ribo", pt.size = 0)
pdf("Plots/Quality_Control/pctribo_vln_byAnimal.pdf", height = 4, width = 6)
print(percent_ribo_vln_by_animal)
dev.off()

# Examine proportions
#plot samples as proportion or percentage of cluster
sample_prop_bar_plot <- ggplot(combined_seurat@meta.data, aes(x = seurat_clusters, fill = Sample_name)) +
  geom_bar(position = "fill") +
  RotatedAxis() +
  xlab("Sample Name")
pdf("Plots/Quality_Control/Sample_prop_barPlot.pdf", height = 4, width = 6)
print(sample_prop_bar_plot)
dev.off()

#Plot samples as proportion or percentage of cluster
cluster_prop_bar_plot <- ggplot(combined_seurat@meta.data, aes(x = Sample_name, fill = seurat_clusters)) +
  geom_bar(position = "fill") +
  RotatedAxis() +
  xlab("Sample Name")
pdf("Plots/Quality_Control/Cluster_prop_barPlot.pdf", height = 4, width = 6)
print(cluster_prop_bar_plot)
dev.off()

#plot proportions of each condition in each cluster
cond_prop_bar_plot <- ggplot(combined_seurat@meta.data, aes(x = seurat_clusters, fill = Treatment)) +
  geom_bar(position = "fill") +
  geom_hline(yintercept = 12 / 33, linetype = "dashed", size = 1) +
  geom_hline(yintercept = 22 / 33, linetype = "dashed", size = 1) +
  RotatedAxis() +
  xlab("Sample Name")
pdf("Plots/Quality_Control/Cond_prop_barPlot.pdf", height = 4, width = 6)
print(cond_prop_bar_plot)
dev.off()

#plot proportions of each sex in each cluster
sex_prop_bar_plot <- ggplot(combined_seurat@meta.data, aes(x = seurat_clusters, fill = Sex)) +
  geom_bar(position = "fill") +
  geom_hline(yintercept = 0.5, linetype = "dashed", size = 1) +
  RotatedAxis() +
  xlab("Cluster")
pdf("Plots/Quality_Control/Sex_prop_barPlot.pdf", height = 4, width = 6)
print(sex_prop_bar_plot)
dev.off()

# Clustered QC plots ####

# Examine the clustered UMAP
Idents(combined_seurat) <- combined_seurat$seurat_clusters

clustered_umap_plot <- DimPlot(combined_seurat)
pdf("Plots/Quality_Control/Clustered_UMAP.pdf", height = 4, width = 6)
print(clustered_umap_plot)
dev.off()

# Visualize QC metrics in each cluster
percent_mito_vln_plot <- VlnPlot(combined_seurat, features = "percent_mito", pt.size = 0)
pdf("Plots/Quality_Control/PercentMito_Violin.pdf")
print(percent_mito_vln_plot)
dev.off()

nFeature_vln_plot <- VlnPlot(combined_seurat, features = "nFeature_RNA", pt.size = 0)
pdf("Plots/Quality_Control/nFeature_Violin.pdf", height = 4, width = 8)
print(nFeature_vln_plot)
dev.off()

nCount_vln_plot <- VlnPlot(combined_seurat, features = "nCount_RNA", pt.size = 0)
pdf("Plots/Quality_Control/nCount_Violin.pdf", height = 4, width = 8)
print(nCount_vln_plot)
dev.off()

doublet_score_vln_plot <- VlnPlot(combined_seurat, features = "Doublet_Score", pt.size = 0)
pdf("Plots/Quality_Control/Doublet_Score_Violin.pdf", height = 4, width = 8)
print(doublet_score_vln_plot)
dev.off()

percent_ribo_vln_plot <- VlnPlot(combined_seurat, features = "percent_ribo", pt.size = 0)
pdf("Plots/Quality_Control/percent_ribo_Violin.pdf", height = 4, width = 8)
print(percent_ribo_vln_plot)
dev.off()

# Look at the marker genes identified by Seurat ####
all_markers <- read.csv("CSV_Results/Marker_Genes_All/All_marker_genes.csv")

# Visualize the top marker gene expression in each cluster
top_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
top_marker_dotplot <- DotPlot(combined_seurat, features = unique(top_markers$gene))+RotatedAxis()
pdf("Plots/Clustering_Plots/Top_Marker_DotPlot.pdf", height = 8, width = 12)
print(top_marker_dotplot)
dev.off()

# Visualize the top 2 marker genes expression in each cluster
top_markers2 <- all_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2markers_dotplot <- DotPlot(combined_seurat, features = unique(top_markers2$gene))+RotatedAxis()
pdf("Plots/Clustering_Plots/Top2_Markers_DotPlot.pdf", height = 8, width = 12)
print(top2markers_dotplot)
dev.off()

# Visualize the top 2 marker genes expression in each cluster
top_markers2_diff <- all_markers %>% group_by(cluster) %>% top_n(n = 2, wt = PCT_Delta)
top2markers_diff_dotplot <- DotPlot(combined_seurat, features = unique(top_markers2_diff$gene))+RotatedAxis()
pdf("Plots/Clustering_Plots/Top2_Markers_diff_DotPlot.pdf", height = 8, width = 12)
print(top2markers_diff_dotplot)
dev.off()


# Marker gene plots ####
# Using markers from Dilly et al 2022 
id_features <- c("Mbp", "Mobp", "Plp1", "Gad1", "Gad2",
                 "Ndrg2", "Slc1a2", "Slc4a4",
                 "Slc17a7", "Satb1", "Neurod6", "Vcan",
                 "Pdgfra", "Pcdh15", "Csf1r", "Apbb1ip", "P2ry12",
                 "Flt1", "B2m", "Bmp4", "Cnp", "Ccdc153",
                 "Rsph1", "Tmem212", "Rbfox3")

major_cells_dotplot <- DotPlot(combined_seurat, features = id_features) + RotatedAxis()
pdf("Plots/Clustering_Plots/Major_cells_dotplot.pdf")
print(major_cells_dotplot)
dev.off()

# Relevant marker genes
# Major Cell Types ####

# Neurons
All_Neurons <- c("Rbfox3", "Syn1", "Gria1")
Inhibitory_Neurons <- c("Gad1", "Gad2")
Excitatory_Neurons <- c("Neurod6", "Satb1", "Slc17a7") # Can be divided into 2 classes
Cholinergic_Neurons <- c("Slc5a7", "Prima1")

# Glia
Astrocytes <- c("Slc1a2", "Slc1a3", "Ndrg2")
Oligodendrocytes <- c("Mog", "Mobp", "Hapln2")
Microglia <- c("Csf1r", "Apbb1ip", "Inpp5d")
OPCs <- c("Vcan", "Pdgfra", "Pcdh15")
Immature_Oligos <- c("Bmp4", "Cnp", "Tcf7l2") # - These might benefit from a more conservative name

# Other Cells
Endothelial_Cells <- c("Flt1", "B2m", "Rgs5")
Ependymal_Cells <- c("Ccdc153", "Rsph1", "Tmem212")

# Lists for programming
DS_cell_types <- c("Astrocytes", "Oligodendrocytes", "Microglia", "OPCs",
                   "Immature_Oligos", "All_Neurons", "Inhibitory_Neurons", "Excitatory_Neurons",
                   "Cholinergic_Neurons", "Endothelial_Cells", "Ependymal_Cells")

major_marker_genes <- c(Astrocytes, Oligodendrocytes, Microglia, OPCs,
                        Immature_Oligos, All_Neurons, Inhibitory_Neurons, Excitatory_Neurons,
                        Cholinergic_Neurons, Endothelial_Cells, Ependymal_Cells)

major_marker_genes_unique <- unique(major_marker_genes)

major_marker_genes_list <- list(Astrocytes, Oligodendrocytes, Microglia, OPCs,
                                Immature_Oligos, All_Neurons, Inhibitory_Neurons, Excitatory_Neurons,
                                Cholinergic_Neurons, Endothelial_Cells, Ependymal_Cells)


# DS Cell Subtypes ####
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

# Examine the resulting UMAP ####
clustered_umap_plot_lbl <- DimPlot(combined_seurat, label = TRUE)
pdf("Plots/Clustering_Plots/Clustered_UMAP.pdf", height = 6, width = 8)
print(clustered_umap_plot)
dev.off()

all_markers_list <- c(major_marker_genes_unique, DS_marker_genes_unique)

# Generate featureplots and violinplots for genes of interest
for (marker in all_markers_list){
  gene_feature_plot <- FeaturePlot(combined_seurat, features = marker)
  pdf(paste0("Plots/Clustering_Plots/Marker_Feature_Plots/", marker, "_FeaturePlot.pdf"), height = 4, width = 6)
  print(gene_feature_plot)
  dev.off()

  gene_violin_plot <- VlnPlot(combined_seurat, features = marker, pt.size = 0)
  pdf(paste0("Plots/Clustering_Plots/Marker_Violin_Plots/", marker, "_ViolinPlot.pdf"), height = 4, width = 6)
  print(gene_violin_plot)
  dev.off()
}

# Log the completion time
write(paste0("06_make_plots - Finish: ", Sys.time()), file = "snRNA_Log.txt", append = TRUE)
