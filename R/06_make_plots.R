# 06_make_plots.R
# Purpose: Make various plots for single-cell analysis
# Author: Geoff Dilly

library(here)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

scRNA_home_dir <- here()
setwd(scRNA_home_dir)

# Setup ####
# Load custom functions
source("R/modules/plot_utils.R")
source("R/modules/log_utils.R")

# Log the start time and a timestamped copy of the script
write(paste0("06_make_plots - Start: ", Sys.time()), file = "scRNA_Log.txt", append = TRUE)
write_script_log("R/06_make_plots.R")

# Load the configuration file
source("sc_experiment_config.R")

# Load the clustered Seurat object
combined_seurat <- readRDS(paste0("R_Data/", scConfig.Prefix, "_combined_clustered.rds"))
DefaultAssay(combined_seurat) <- "SCT"

# Set the identity to label clusters
cluster_ident <- "seurat_clusters"

# Overall QC: Ident = Project_name ####
Idents(combined_seurat) <- combined_seurat$orig.ident

# Quality metric violin plots
qc_vlns_plot <- VlnPlot(
  combined_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "Doublet_Score"),
  pt.size = 0, ncol = 5
) + NoLegend() & theme(axis.text.x = element_blank())

save_plot_pdf(qc_vlns_plot, "Plots/Quality_Control/QC_VlnPlot.pdf", height = 3, width = 18)

# Examine QC metrics on the UMAP plot
qc_umap_plot <- FeaturePlot(
  combined_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "Doublet_Score"),
  ncol = 5
) & theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

save_plot_pdf(qc_umap_plot, "Plots/Quality_Control/QC_UMAP.pdf", height = 3, width = 18)

# Examine doublet score in doublets and non doublets
if (scConfig.remove_doublets == FALSE) {
  doublet_call_umap_plot <- FeaturePlot(
    combined_seurat,
    features = "Doublet_Score",
    split.by = "Doublet_Call"
  )
  save_plot_pdf(doublet_call_umap_plot, "Plots/Quality_Control/Doublet_Call_UMAP.pdf", height = 4, width = 10)
}

# Sex-level QC: Ident = Sex ####
Idents(combined_seurat) <- combined_seurat$Sex

# Plot by sex on the UMAP
sex_umap_plot <- DimPlot(combined_seurat)
save_plot_pdf(sex_umap_plot, "Plots/Quality_Control/Sex_UMAP.pdf", height = 4, width = 6)

# Quality metric violin plots
qc_bysex_vlns_plot <- VlnPlot(
  combined_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "Doublet_Score"),
  pt.size = 0, stack = TRUE, flip = TRUE
) + NoLegend() & theme(axis.text.x = element_blank(), axis.title = element_blank())

save_plot_pdf(qc_bysex_vlns_plot, "Plots/Quality_Control/QC_bySex_VlnPlot.pdf")

# Treatment-level QC: Ident = Treatment ####
Idents(combined_seurat) <- combined_seurat$Treatment

# Plot by sex on the UMAP
treatment_umap_plot <- DimPlot(combined_seurat)
save_plot_pdf(treatment_umap_plot, "Plots/Quality_Control/Treatment_UMAP.pdf", height = 4, width = 6)

# Quality metric violin plots
qc_bytreatment_vlns_plot <- VlnPlot(
  combined_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "Doublet_Score"),
  pt.size = 0, ncol = 5
) & theme(axis.title.x = element_blank())

save_plot_pdf(qc_bytreatment_vlns_plot, "Plots/Quality_Control/QC_byTreatment_VlnPlot.pdf", height = 3, width = 18)

# Sample-level QC: Ident = Sample_name ####
Idents(combined_seurat) <- combined_seurat$Sample_name

# Plot by sample on the UMAP
sample_umap_plot <- DimPlot(combined_seurat)
save_plot_pdf(sample_umap_plot, "Plots/Quality_Control/Sample_UMAP.pdf", height = 4, width = 6)

# By sample QC violin plots
qc_by_sample_vln_plot <- VlnPlot(
  combined_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "Doublet_Score"),
  pt.size = 0, stack = TRUE, flip = TRUE
) + NoLegend()

save_plot_pdf(qc_by_sample_vln_plot, "Plots/Quality_Control/QC_bySample_VlnPlot.pdf", height = 4, width = 12)

# Scatterplot: nCount_RNA vs nFeature_RNA, colored by sample
scatter_nc_vs_nf <- ggplot(
  combined_seurat@meta.data,
  aes(x = nFeature_RNA, y = nCount_RNA, color = Sample_name)
) +
  geom_point(alpha = 0.4, size = 1) +
  theme_classic() +
  xlab("Number of detected genes (nFeature_RNA)") +
  ylab("Number of UMIs (nCount_RNA)") +
  ggtitle("nCount_RNA vs nFeature_RNA by Sample") +
  theme(
    legend.title = element_blank(),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    plot.title = element_text(size = 14, face = "bold")
  ) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))

save_plot_pdf(scatter_nc_vs_nf, "Plots/Quality_Control/Scatter_nCount_vs_nFeature_bySample.pdf", height = 5, width = 7)

# Cluster-level QC: Ident = Seurat_clusters ####
Idents(combined_seurat) <- cluster_ident

# Clustered QC plots
# Plot by cluster on the UMAP
sample_umap_plot <- DimPlot(combined_seurat)
save_plot_pdf(sample_umap_plot, "Plots/Quality_Control/Cluster_UMAP.pdf", height = 4, width = 6)
save_plot_pdf(sample_umap_plot, "Plots/Clustering_Plots/Cluster_UMAP.pdf", height = 4, width = 6)

sample_umap_plot <- DimPlot(combined_seurat, label = TRUE)
save_plot_pdf(sample_umap_plot, "Plots/Quality_Control/Cluster_UMAP_labeled.pdf", height = 4, width = 6)
save_plot_pdf(sample_umap_plot, "Plots/Clustering_Plots/Cluster_UMAP.pdf", height = 4, width = 6)

# Proportion bar plots
# Cell count barplot per sample (counts as height)
cell_count_barplot <- ggplot(
  combined_seurat@meta.data,
  aes(x = Sample_name, fill = .data[[cluster_ident]])
) +
  geom_bar() +
  theme_classic() +
  xlab("Sample Name") +
  ylab("Cell Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")

save_plot_pdf(cell_count_barplot, "Plots/Quality_Control/Cell_Counts_Barplot.pdf", height = 4, width = 6)

# Plot samples as proportion or percentage of cluster
sample_prop_bar_plot <- ggplot(
  combined_seurat@meta.data,
  aes(x = .data[[cluster_ident]], fill = Sample_name)
) +
  geom_bar(position = "fill") +
  RotatedAxis() +
  xlab("Sample Name")

save_plot_pdf(sample_prop_bar_plot, "Plots/Quality_Control/Sample_prop_barPlot.pdf", height = 4, width = 6)

# Plot samples as proportion or percentage of cluster
cluster_prop_bar_plot <- ggplot(
  combined_seurat@meta.data,
  aes(x = Sample_name, fill = .data[[cluster_ident]])
) +
  geom_bar(position = "fill") +
  RotatedAxis() +
  xlab("Sample Name")

save_plot_pdf(cluster_prop_bar_plot, "Plots/Quality_Control/Cluster_prop_barPlot.pdf", height = 4, width = 6)

# Plot proportions of each condition in each cluster
cond_prop_bar_plot <- ggplot(
  combined_seurat@meta.data,
  aes(x = .data[[cluster_ident]], fill = Treatment)
) +
  geom_bar(position = "fill") +
  geom_hline(yintercept = 0.5, linetype = "dashed", size = 1) +
  RotatedAxis() +
  xlab("Sample Name")

save_plot_pdf(cond_prop_bar_plot, "Plots/Quality_Control/Cond_prop_barPlot.pdf", height = 4, width = 6)

# Plot proportions of each sex in each cluster
sex_prop_bar_plot <- ggplot(
  combined_seurat@meta.data,
  aes(x = .data[[cluster_ident]], fill = Sex)
) +
  geom_bar(position = "fill") +
  geom_hline(yintercept = 0.5, linetype = "dashed", size = 1) +
  RotatedAxis() +
  xlab("Cluster")

save_plot_pdf(sex_prop_bar_plot, "Plots/Quality_Control/Sex_prop_barPlot.pdf", height = 4, width = 6)

# By cluster QC violin plots
qc_by_cluster_vln_plot <- VlnPlot(
  combined_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "Doublet_Score"),
  pt.size = 0, stack = TRUE, flip = TRUE
) + NoLegend()

save_plot_pdf(QC_bycluster_vln_plot, "Plots/Quality_Control/QC_byCluster_VlnPlot.pdf", height = 4, width = 12)

# Marker gene plots ####
# Load the marker genes previously identified with FindAllMarkers()
all_markers <- read.csv("CSV_Results/Marker_Genes_All/All_marker_genes.csv")

# Visualize the top marker gene expression in each cluster
top_markers <- all_markers %>% group_by(cluster) %>% top_n(n = 1, wt = avg_log2FC)
top_marker_dotplot <- DotPlot(combined_seurat, features = unique(top_markers$gene)) + RotatedAxis()
save_plot_pdf(top_marker_dotplot, "Plots/Clustering_Plots/Top_Marker_DotPlot.pdf", height = 12, width = 18)

# Visualize the top 2 marker genes expression in each cluster
top_markers2 <- all_markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_log2FC)
top2markers_dotplot <- DotPlot(combined_seurat, features = unique(top_markers2$gene)) + RotatedAxis()
save_plot_pdf(top2markers_dotplot, "Plots/Clustering_Plots/Top2_Markers_DotPlot.pdf", height = 12, width = 18)

# Make a stacked violin plot for the DS marker genes
marker_stack_vln_plot <- make_stacked_vln_plot(
  seurat_obj = combined_seurat,
  features = unique(top_markers2$gene),
  assay = "SCT",
  slot = "data",
  adjust = 2.5,
  pt.size = 0,
  gene_label_size = 8
)
save_plot_pdf(marker_stack_vln_plot, "Plots/Clustering_Plots/Top_Marker_stacked_vln_plot.pdf", height = 12, width = 4)

# Load markers from Dilly et al 2022
major_marker_genes_list <- c(
  "Mbp", "Mobp", "Plp1", "Gad1", "Gad2", "Ndrg2", "Slc1a2", "Slc4a4",
  "Slc17a7", "Satb1", "Neurod6", "Vcan", "Pdgfra", "Pcdh15", "Csf1r", "Apbb1ip", "P2ry12",
  "Flt1", "B2m", "Bmp4", "Cnp", "Ccdc153", "Rsph1", "Tmem212", "Rbfox3"
)
major_marker_genes_unique <- unique(unlist(major_marker_genes_list))

major_cells_dotplot <- DotPlot(combined_seurat, features = major_marker_genes_unique) + RotatedAxis()
save_plot_pdf(major_cells_dotplot, "Plots/Clustering_Plots/Major_cells_dotplot.pdf", height = 4, width = 6)

# Make a stacked violin plot for the DS marker genes
brain_marker_stack_vln_plot <- make_stacked_vln_plot(
  seurat_obj = combined_seurat,
  features = major_marker_genes_unique,
  assay = "SCT",
  slot = "data",
  adjust = 2.5,
  pt.size = 0,
  gene_label_size = 8
)
save_plot_pdf(brain_marker_stack_vln_plot, "Plots/Clustering_Plots/Brainmarker_stacked_vln_plot.pdf", height = 12, width = 4)

# DS Cell Subtypes
# From Wildermuth et al 2025
dSPN      <- c("Drd1", "Pdyn")
iSPN      <- c("Drd2", "Adora2a")
eSPN      <- c("Otof", "Col11a1")
PVALB_IN  <- c("Kit", "Pvalb")
Chol_IN   <- c("Tacr1", "Chat")
SST_IN    <- c("Sst", "Npy")
Astro     <- c("Slc1a3", "Rorb")
Oligo     <- c("Mog", "Aspa")
OPC       <- c("Pdgfra")
Endo      <- c("Flt1", "Slco1a4")
Ependy    <- c("Dnah12", "Rsph1", "Gm973")
Mural     <- c("Pdgfrb", "Rgs5")
Microglia <- c("C1qc", "Cx3cr1")

DS_cell_types <- c(
  dSPN, iSPN, eSPN, PVALB_IN, Chol_IN, SST_IN,
  Astro, Oligo, OPC, Endo, Ependy, Mural, Microglia
)

DS_marker_genes_list <- list(
  dSPN, iSPN, eSPN, PVALB_IN, Chol_IN, SST_IN,
  Astro, Oligo, OPC, Endo, Ependy, Mural, Microglia
)
DS_marker_genes_unique <- unique(unlist(DS_marker_genes_list))

DS_cells_dotplot <- DotPlot(combined_seurat, features = DS_marker_genes_unique) + RotatedAxis()
save_plot_pdf(DS_cells_dotplot, "Plots/Clustering_Plots/DS_markers_dotplot.pdf", height = 4, width = 6)

# Make a stacked violin plot for the DS marker genes
DS_marker_stack_vln_plot <- make_stacked_vln_plot(
  seurat_obj = combined_seurat,
  features = DS_marker_genes_unique,
  assay = "SCT",
  slot = "data",
  adjust = 2.5,
  pt.size = 0,
  gene_label_size = 8
)
save_plot_pdf(DS_marker_stack_vln_plot, "Plots/Clustering_Plots/DS_marker_stacked_vln_plot.pdf", height = 12, width = 4)

# Generate featureplots and violinplots for genes of interest
all_markers_list <- c(major_marker_genes_unique, DS_marker_genes_unique)

for (marker in all_markers_list) {
  gene_feature_plot <- FeaturePlot(combined_seurat, features = marker)
  save_plot_pdf(
    gene_feature_plot,
    paste0("Plots/Clustering_Plots/Marker_Feature_Plots/", marker, "_FeaturePlot.pdf"),
    height = 4, width = 6
  )

  gene_violin_plot <- VlnPlot(combined_seurat, features = marker, pt.size = 0) + NoLegend()
  save_plot_pdf(
    gene_violin_plot,
    paste0("Plots/Clustering_Plots/Marker_Violin_Plots/", marker, "_ViolinPlot.pdf"),
    height = 4, width = 6
  )
}

# Log the completion time
write(paste0("06_make_plots - Finish: ", Sys.time()), file = "scRNA_Log.txt", append = TRUE)
