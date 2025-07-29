# 06_make_plots.R
# Purpose: Make various plots for single-cell analysis
# Author: Geoff Dilly

library(here)
library(yaml)
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# Setup ####
# Load custom functions
source(here::here("R/modules/run_utils.R"))
source(here::here("R/modules/qc_utils.R"))
source(here::here("R/modules/plot_utils"))

# Load the configuration file and metadata
scConfig <- yaml::read_yaml(here::here("sc_experiment_config.yaml"))
scConfig$Sample_metadata <- read.csv(here::here("sc_sample_metadata.csv"))

# Check for required directories
check_required_dirs()

# Function to get results directory based on run time 
output_dir <- get_results_dir(run_time = Sys.getenv("RUN_TIME"), 
                           prefix = scConfig$prefix)

# Log the start time and a timestamped copy of the script
write(paste0("06_make_plots - Start: ", Sys.time()), file = here::here("scRNA_Log.txt"), append = TRUE)
log_connection <- write_script_log(here::here("R/06_make_plots.R"), log_dir = here::here(output_dir, "logs"))

# Log all output to the end of the log file
sink(log_connection, append = TRUE)
sink(log_connection, type = "message", append = TRUE)
on.exit({
  sink(NULL)
  sink(NULL, type = "message")
})

# Load data ####
# Load the clustered Seurat object
combined_seurat <- readRDS(here::here("data", "R_Data", paste0(scConfig$prefix, "_combined_clustered.rds")))
DefaultAssay(combined_seurat) <- "SCT"

# Set the identity to label clusters
cluster_col <- scConfig$cluster_plot_ident

# Ensure cluster_col exists in meta.data, default to seurat_clusters if not
if (!cluster_col %in% colnames(combined_seurat@meta.data)) {
  message("cluster_col not found in metadata; defaulting to 'seurat_clusters'")
  cluster_col <- "seurat_clusters"
}

# Overall QC: Ident = project_name ####
Idents(combined_seurat) <- combined_seurat$orig.ident

# Quality metric violin plots
qc_vlns_plot <- VlnPlot(
  combined_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "Doublet_Score"),
  pt.size = 0, ncol = 5
) + NoLegend() & theme(axis.text.x = element_blank())

save_plot_pdf(qc_vlns_plot, here::here(output_dir, "plots", "Quality_Control", "QC_VlnPlot.pdf"), height = 3, width = 18)

# Examine QC metrics on the UMAP plot
qc_umap_plot <- FeaturePlot(
  combined_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "Doublet_Score"),
  ncol = 5
) & theme(axis.title = element_blank(), axis.text = element_blank(), axis.ticks = element_blank())

save_plot_pdf(qc_umap_plot, here::here(output_dir, "plots", "Quality_Control", "QC_UMAP.pdf"), height = 3, width = 18)

# Examine doublet score in doublets and non doublets
if (scConfig$remove_doublets == FALSE) {
  doublet_call_umap_plot <- FeaturePlot(
    combined_seurat,
    features = "Doublet_Score",
    split.by = "Doublet_Call"
  )
  save_plot_pdf(doublet_call_umap_plot, here::here(output_dir, "plots", "Quality_Control", "Doublet_Call_UMAP.pdf"), height = 4, width = 10)
}

# Sex-level QC: Ident = Sex ####
Idents(combined_seurat) <- combined_seurat$Sex

# Plot by sex on the UMAP
sex_umap_plot <- DimPlot(combined_seurat)
save_plot_pdf(sex_umap_plot, here::here(output_dir, "plots", "Quality_Control", "Sex_UMAP.pdf"), height = 4, width = 6)

# Quality metric violin plots
qc_bysex_vlns_plot <- VlnPlot(
  combined_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "Doublet_Score"),
  pt.size = 0, stack = TRUE, flip = TRUE
) + NoLegend() & theme(axis.text.x = element_blank(), axis.title = element_blank())

save_plot_pdf(qc_bysex_vlns_plot, here::here(output_dir, "plots", "Quality_Control", "QC_bySex_VlnPlot.pdf"))

# Treatment-level QC: Ident = Treatment ####
Idents(combined_seurat) <- combined_seurat$Treatment

# Plot by sex on the UMAP
treatment_umap_plot <- DimPlot(combined_seurat)
save_plot_pdf(treatment_umap_plot, here::here(output_dir, "plots", "Quality_Control", "Treatment_UMAP.pdf"), height = 4, width = 6)

# Quality metric violin plots
qc_bytreatment_vlns_plot <- VlnPlot(
  combined_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "Doublet_Score"),
  pt.size = 0, ncol = 5
) & theme(axis.title.x = element_blank())

save_plot_pdf(qc_bytreatment_vlns_plot, here::here(output_dir, "plots", "Quality_Control", "QC_byTreatment_VlnPlot.pdf"), height = 3, width = 18)

# Sample-level QC: Ident = Sample_name ####
Idents(combined_seurat) <- combined_seurat$Sample_name

# Plot by sample on the UMAP
sample_umap_plot <- DimPlot(combined_seurat)
save_plot_pdf(sample_umap_plot, here::here(output_dir, "plots", "Quality_Control", "Sample_UMAP.pdf"), height = 4, width = 6)

# By sample QC violin plots
qc_by_sample_vln_plot <- VlnPlot(
  combined_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "Doublet_Score"),
  pt.size = 0, stack = TRUE, flip = TRUE
) + NoLegend()

save_plot_pdf(qc_by_sample_vln_plot, here::here(output_dir, "plots", "Quality_Control", "QC_bySample_VlnPlot.pdf"), height = 4, width = 12)

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

save_plot_pdf(scatter_nc_vs_nf, here::here(output_dir, "plots", "Quality_Control", "Scatter_nCount_vs_nFeature_bySample.pdf"), height = 5, width = 7)

# Cluster-level QC: Ident = Seurat_clusters ####
Idents(combined_seurat) <- cluster_col

# Clustered QC plots
# Plot by cluster on the UMAP
sample_umap_plot <- DimPlot(combined_seurat)
save_plot_pdf(sample_umap_plot, here::here(output_dir, "plots", "Quality_Control", "Cluster_UMAP.pdf"), height = 4, width = 6)
save_plot_pdf(sample_umap_plot, here::here(output_dir, "plots", "Cluster_Plots", "Cluster_UMAP.pdf"), height = 4, width = 6)

sample_umap_plot_lbl <- DimPlot(combined_seurat, label = TRUE)
save_plot_pdf(sample_umap_plot_lbl, here::here(output_dir, "plots", "Quality_Control", "Cluster_UMAP_labeled.pdf"), height = 4, width = 6)
save_plot_pdf(sample_umap_plot_lbl, here::here(output_dir, "plots", "Cluster_Plots", "Cluster_UMAP.pdf"), height = 4, width = 6)

# Proportion bar plots
cell_count_barplot <- ggplot(
  combined_seurat@meta.data,
  aes(x = Sample_name, fill = .data[[cluster_col]])
) +
  geom_bar() +
  theme_classic() +
  xlab("Sample Name") +
  ylab("Cell Count") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  guides(fill = "none")

save_plot_pdf(cell_count_barplot, here::here(output_dir, "plots", "Quality_Control", "Cell_Counts_Barplot.pdf"), height = 4, width = 6)

sample_prop_bar_plot <- ggplot(
  combined_seurat@meta.data,
  aes(x = .data[[cluster_col]], fill = Sample_name)
) +
  geom_bar(position = "fill") +
  RotatedAxis() +
  xlab("Sample Name")

save_plot_pdf(sample_prop_bar_plot, here::here(output_dir, "plots", "Quality_Control", "Sample_prop_barPlot.pdf"), height = 4, width = 6)

cluster_prop_bar_plot <- ggplot(
  combined_seurat@meta.data,
  aes(x = Sample_name, fill = .data[[cluster_col]])
) +
  geom_bar(position = "fill") +
  RotatedAxis() +
  xlab("Sample Name")

save_plot_pdf(cluster_prop_bar_plot, here::here(output_dir, "plots", "Quality_Control", "Cluster_prop_barPlot.pdf"), height = 4, width = 6)

cond_prop_bar_plot <- ggplot(
  combined_seurat@meta.data,
  aes(x = .data[[cluster_col]], fill = Treatment)
) +
  geom_bar(position = "fill") +
  geom_hline(yintercept = 0.5, linetype = "dashed", size = 1) +
  RotatedAxis() +
  xlab("Sample Name")

save_plot_pdf(cond_prop_bar_plot, here::here(output_dir, "plots", "Quality_Control", "Cond_prop_barPlot.pdf"), height = 4, width = 6)

sex_prop_bar_plot <- ggplot(
  combined_seurat@meta.data,
  aes(x = .data[[cluster_col]], fill = Sex)
) +
  geom_bar(position = "fill") +
  geom_hline(yintercept = 0.5, linetype = "dashed", size = 1) +
  RotatedAxis() +
  xlab("Cluster")

save_plot_pdf(sex_prop_bar_plot, here::here(output_dir, "plots", "Quality_Control", "Sex_prop_barPlot.pdf"), height = 4, width = 6)

# By cluster QC violin plots
qc_by_cluster_vln_plot <- VlnPlot(
  combined_seurat,
  features = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_ribo", "Doublet_Score"),
  pt.size = 0, stack = TRUE, flip = TRUE
) + NoLegend()

save_plot_pdf(qc_by_cluster_vln_plot, here::here(output_dir, "plots", "Quality_Control", "QC_byCluster_VlnPlot.pdf"), height = 4, width = 12)

# Marker gene plots ####
all_markers <- read.csv(here::here(output_dir, "csv_results", "Marker_Genes_All", "All_marker_genes.csv"))

top_markers <- all_markers %>% group_by(cluster) %>% slice_max(n = 1, order_by = avg_log2FC)
top_marker_dotplot <- DotPlot(combined_seurat, features = unique(top_markers$gene)) + RotatedAxis()
save_plot_pdf(top_marker_dotplot, here::here(output_dir, "plots", "Cluster_Plots", "Top_Marker_DotPlot.pdf"), height = 12, width = 18)

top_markers2 <- all_markers %>% group_by(cluster) %>% slice_max(n = 2, order_by = avg_log2FC)
top2markers_dotplot <- DotPlot(combined_seurat, features = unique(top_markers2$gene)) + RotatedAxis()
save_plot_pdf(top2markers_dotplot, here::here(output_dir, "plots", "Cluster_Plots", "Top2_Markers_DotPlot.pdf"), height = 12, width = 18)

marker_stack_vln_plot <- make_stacked_vln_plot(
  seurat_obj = combined_seurat,
  features = unique(top_markers2$gene),
  assay = "SCT",
  slot = "data",
  adjust = 2.5,
  pt.size = 0,
  gene_label_size = 8
)
save_plot_pdf(marker_stack_vln_plot, here::here(output_dir, "plots", "Cluster_Plots", "Top_Marker_stacked_vln_plot.pdf"), height = 12, width = 4)

marker_tbl <- read.csv(here::here("reference", "marker_gene_db.csv"), stringsAsFactors = FALSE)
major_marker_genes <- marker_tbl %>% filter(reference == scConfig$marker_gene_reference) %>% pull(gene)

major_marker_genes_unique <- unique(unlist(major_marker_genes))
major_cells_dotplot <- DotPlot(combined_seurat, features = major_marker_genes_unique) + RotatedAxis()
save_plot_pdf(major_cells_dotplot, here::here(output_dir, "plots", "Cluster_Plots", "Major_cells_dotplot.pdf"), height = 10, width = 12)

brain_marker_stack_vln_plot <- make_stacked_vln_plot(
  seurat_obj = combined_seurat,
  features = major_marker_genes_unique,
  assay = "SCT",
  slot = "data",
  adjust = 2.5,
  pt.size = 0,
  gene_label_size = 8
)
save_plot_pdf(brain_marker_stack_vln_plot, here::here(output_dir, "plots", "Cluster_Plots", "Brainmarker_stacked_vln_plot.pdf"), height = 12, width = 4)

all_markers_list <- c(major_marker_genes_unique)

for (marker in all_markers_list) {
  gene_feature_plot <- FeaturePlot(combined_seurat, features = marker)
  save_plot_pdf(
    gene_feature_plot,
    here::here(output_dir, "plots", "Cluster_Plots", "Marker_Feature_Plots", paste0(marker, "_FeaturePlot.pdf")),
    height = 6, width = 8
  )

  gene_violin_plot <- VlnPlot(combined_seurat, features = marker, pt.size = 0) + NoLegend()
  save_plot_pdf(
    gene_violin_plot,
    here::here(output_dir, "plots", "Cluster_Plots", "Marker_Violin_Plots", paste0(marker, "_ViolinPlot.pdf")),
    height = 6, width = 8
  )
}

cell_types_dotgrid <- dotplot_by_marker_group(
  seurat_obj = combined_seurat,
  marker_csv = here::here("reference", "marker_gene_db.csv"),
  reference_name = "Dilly_et_al_2022",
  group_by = "seurat_clusters",
  gene_group_col = "cell_type",
  add_separators = TRUE
)
save_plot_pdf(cell_types_dotgrid, here::here(output_dir, "plots", "Cluster_Plots", "Brainmarker_dotgrid_plot.pdf"), height = 12, width = 10)

# Log the completion time
write(paste0("06_make_plots - Finish: ", Sys.time()), file = here::here("scRNA_Log.txt"), append = TRUE)
