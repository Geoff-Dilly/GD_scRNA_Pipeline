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

# Load custom functions
source("R/modules/plot_utils.R")

# Log the start time and a timestamped copy of the script
write(paste0("06_make_plots - Start: ", Sys.time()), file = "snRNA_Log.txt", append = TRUE)
file.copy("R/06_make_plots.R", paste0("Logs/Time_", format(Sys.time(), "%Y-%m-%d_%H-%M-%S"), "_", "06_make_plots.R"), overwrite = FALSE)

# Load the configuration file
source("sc_experiment_config.R")

# Load the clustered Seurat object
combined_seurat <- readRDS(paste0("R_Data/", scConfig.Prefix, "_combined_clustered.rds"))
DefaultAssay(combined_seurat) <- "RNA"

# Make all of the figures from previous scripts ####
# Check quality metrics for each cell
qc_vlns_plot <- VlnPlot(combined_seurat, features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "Doublet_Score"),
pt.size = 0, stack = TRUE, flip = TRUE)
save_plot_pdf(qc_vlns_plot, "Plots/Quality_Control/QC_VlnPlot.pdf", height = 6, width = 10)

# Examine the raw UMAP
raw_umap_plot <- DimPlot(combined_seurat)
save_plot_pdf(raw_umap_plot, "Plots/Clustering_Plots/Raw_UMAP.pdf", height = 4, width = 6)

# Examine QC metrics on the UMAP plot
qc_umap_plot <- FeaturePlot(combined_seurat, features = c("percent_mito", "nFeature_RNA", "Doublet_Score"))
save_plot_pdf(qc_umap_plot, "Plots/Quality_Control/QC_UMAP.pdf", height = 4, width = 10)

# Examine percent ribosomal RNA on a UMAP
ribo_umap_plot <- FeaturePlot(combined_seurat, features = "percent_ribo")
save_plot_pdf(ribo_umap_plot, "Plots/Quality_Control/Ribo_UMAP.pdf", height = 4, width = 6)

# Examine doublet score in doublets and non doublets
doublet_call_umap_plot <- FeaturePlot(combined_seurat, features = "Doublet_Score", split.by = "Doublet_Call")
save_plot_pdf(doublet_call_umap_plot, "Plots/Quality_Control/Doublet_Call_UMAP.pdf", height = 4, width = 10)

# Examine QC metrics by animal
Idents(combined_seurat) <- combined_seurat$Sample_name

# Umap
sample_umap_plot <- DimPlot(combined_seurat)
save_plot_pdf(sample_umap_plot, "Plots/Quality_Control/Sample_UMAP.pdf", height = 4, width = 6)

# Violin plots
nFeature_vln_by_animal <- VlnPlot(combined_seurat, features = "nFeature_RNA", pt.size = 0)
save_plot_pdf(nFeature_vln_by_animal, "Plots/Quality_Control/nFeature_vln_byAnimal.pdf", height = 4, width = 6)

nCount_vln_by_animal <- VlnPlot(combined_seurat, features = "nCount_RNA", pt.size = 0)
save_plot_pdf(nCount_vln_by_animal, "Plots/Quality_Control/nCount_vln_byAnimal.pdf", height = 4, width = 6)

mito_vln_by_animal <- VlnPlot(combined_seurat, features = "percent_mito", pt.size = 0)
save_plot_pdf(mito_vln_by_animal, "Plots/Quality_Control/percentMito_vln_byAnimal.pdf", height = 4, width = 6)

ribo_vln_by_animal <- VlnPlot(combined_seurat, features = "percent_ribo", pt.size = 0)
save_plot_pdf(ribo_vln_by_animal, "Plots/Quality_Control/percentRibo_vln_byAnimal.pdf", height = 4, width = 6)

# Examine proportions
#plot samples as proportion or percentage of cluster
sample_prop_bar_plot <- ggplot(combined_seurat@meta.data, aes(x = seurat_clusters, fill = Sample_name)) +
  geom_bar(position = "fill") +
  RotatedAxis() +
  xlab("Sample Name")
save_plot_pdf(sample_prop_bar_plot, "Plots/Quality_Control/Sample_prop_barPlot.pdf", height = 4, width = 6)

#Plot samples as proportion or percentage of cluster
cluster_prop_bar_plot <- ggplot(combined_seurat@meta.data, aes(x = Sample_name, fill = seurat_clusters)) +
  geom_bar(position = "fill") +
  RotatedAxis() +
  xlab("Sample Name")
save_plot_pdf(cluster_prop_bar_plot, "Plots/Quality_Control/Cluster_prop_barPlot.pdf", height = 4, width = 6)

#plot proportions of each condition in each cluster
cond_prop_bar_plot <- ggplot(combined_seurat@meta.data, aes(x = seurat_clusters, fill = Treatment)) +
  geom_bar(position = "fill") +
  geom_hline(yintercept = 12 / 33, linetype = "dashed", size = 1) +
  geom_hline(yintercept = 22 / 33, linetype = "dashed", size = 1) +
  RotatedAxis() +
  xlab("Sample Name")
save_plot_pdf(cond_prop_bar_plot, "Plots/Quality_Control/Cond_prop_barPlot.pdf", height = 4, width = 6)

#plot proportions of each sex in each cluster
sex_prop_bar_plot <- ggplot(combined_seurat@meta.data, aes(x = seurat_clusters, fill = Sex)) +
  geom_bar(position = "fill") +
  geom_hline(yintercept = 0.5, linetype = "dashed", size = 1) +
  RotatedAxis() +
  xlab("Cluster")
save_plot_pdf(sex_prop_bar_plot, "Plots/Quality_Control/Sex_prop_barPlot.pdf", height = 4, width = 6)

# Clustered QC plots ####

# Examine the clustered UMAP
Idents(combined_seurat) <- combined_seurat$seurat_clusters

clustered_umap_plot <- DimPlot(combined_seurat)
save_plot_pdf(clustered_umap_plot, "Plots/Clustering_Plots/Clustered_UMAP.pdf", height = 4, width = 6)

# Visualize QC metrics in each cluster
percent_mito_vln_plot <- VlnPlot(combined_seurat, features = "percent_mito", pt.size = 0)
save_plot_pdf(percent_mito_vln_plot, "Plots/Clustering_Plots/PercentMito_Violin.pdf", height = 4, width = 6)

nFeature_vln_plot <- VlnPlot(combined_seurat, features = "nFeature_RNA", pt.size = 0)
save_plot_pdf(nFeature_vln_plot, "Plots/Clustering_Plots/nFeature_Violin.pdf", height = 4, width = 6)

nCount_vln_plot <- VlnPlot(combined_seurat, features = "nCount_RNA", pt.size = 0)
save_plot_pdf(nCount_vln_plot, "Plots/Clustering_Plots/nCount_Violin.pdf", height = 4, width = 6)

doublet_score_vln_plot <- VlnPlot(combined_seurat, features = "Doublet_Score", pt.size = 0)
save_plot_pdf(doublet_score_vln_plot, "Plots/Quality_Control/Doublet_Violin.pdf", height = 4, width = 6)

ribo_vln_plot <- VlnPlot(combined_seurat, features = "percent_ribo", pt.size = 0)
save_plot_pdf(ribo_vln_plot, "Plots/Quality_Control/Ribo_Violin.pdf", height = 4, width = 6)

# Look at the marker genes identified by Seurat ####
all_markers <- read.csv("CSV_Results/Marker_Genes_All/All_marker_genes.csv")

# Visualize the top marker gene expression in each cluster
top_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
top_marker_dotplot <- DotPlot(combined_seurat, features = unique(top_markers$gene))+RotatedAxis()
save_plot_pdf(top_marker_dotplot, "Plots/Clustering_Plots/Top_Marker_DotPlot.pdf", height = 8, width = 12)

# Visualize the top 2 marker genes expression in each cluster
top_markers2 <- all_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2markers_dotplot <- DotPlot(combined_seurat, features = unique(top_markers2$gene))+RotatedAxis()
save_plot_pdf(top2markers_dotplot, "Plots/Clustering_Plots/Top2_Markers_DotPlot.pdf", height = 8, width = 12)

# Visualize the top 2 marker genes expression in each cluster
top_markers2_diff <- all_markers %>% group_by(cluster) %>% top_n(n = 2, wt = PCT_Delta)
top2markers_diff_dotplot <- DotPlot(combined_seurat, features = unique(top_markers2_diff$gene))+RotatedAxis()
save_plot_pdf(top2markers_diff_dotplot, "Plots/Clustering_Plots/Top2_Markers_Diff_DotPlot.pdf", height = 8, width = 12)


# Marker gene plots ####
# Using markers from Dilly et al 2022 
id_features <- c("Mbp", "Mobp", "Plp1", "Gad1", "Gad2",
                 "Ndrg2", "Slc1a2", "Slc4a4",
                 "Slc17a7", "Satb1", "Neurod6", "Vcan",
                 "Pdgfra", "Pcdh15", "Csf1r", "Apbb1ip", "P2ry12",
                 "Flt1", "B2m", "Bmp4", "Cnp", "Ccdc153",
                 "Rsph1", "Tmem212", "Rbfox3")

major_cells_dotplot <- DotPlot(combined_seurat, features = id_features) + RotatedAxis()
save_plot_pdf(major_cells_dotplot, "Plots/Clustering_Plots/Major_cells_dotplot.pdf", height = 4, width = 6)

# DS Cell Subtypes ####
# From
dSPN <- c("Drd1", "Pdyn")
iSPN <- c("Drd2", "Adora2a")
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

DS_marker_genes_unique <- unique(unlist(DS_marker_genes_list))

# Examine the resulting UMAP ####
clustered_umap_plot_lbl <- DimPlot(combined_seurat, label = TRUE)
save_plot_pdf(clustered_umap_plot_lbl, "Plots/Clustering_Plots/Clustered_UMAP_Labeled.pdf", height = 4, width = 6)

all_markers_list <- c(major_marker_genes_unique, DS_marker_genes_unique)

# Generate featureplots and violinplots for genes of interest
for (marker in all_markers_list){
  gene_feature_plot <- FeaturePlot(combined_seurat, features = marker)
  save_plot_pdf(gene_feature_plot, paste0("Plots/Clustering_Plots/Marker_Feature_Plots/", marker, "_FeaturePlot.pdf"), height = 4, width = 6)

  gene_violin_plot <- VlnPlot(combined_seurat, features = marker, pt.size = 0)
  save_plot_pdf(gene_violin_plot, paste0("Plots/Clustering_Plots/Marker_Violin_Plots/", marker, "_ViolinPlot.pdf"), height = 4, width = 6)
}

# Make a stacked violin plot for the DS marker genes
make_stacked_vln_plot(seurat_obj = combined_seurat,
                      features = DS_marker_genes_unique,
                      assay = "SCT",
                      slot = "data",
                      adjust = 2.5,
                      pt.size = 0,
                      gene_label_size = 8)
save_plot_pdf(stackers, "Plots/Clustering_Plots/marker_stacked_vln_plot.pdf", height = 12, width = 4)

# Log the completion time
write(paste0("06_make_plots - Finish: ", Sys.time()), file = "snRNA_Log.txt", append = TRUE)
