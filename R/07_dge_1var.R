# 07_dge_1var.R
# Purpose: Perform DESEQ2 analysis on pseudobulked Seurat objects and make plots
# Author: Geoff Dilly

library(here)
library(yaml)
library(foreach)
library(doParallel)
library(Seurat)
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(pheatmap)

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
output_dir <- Get_results_dir(run_time = Sys.getenv("RUN_TIME"), 
                            prefix = scConfig$prefix)

# Log the start time and a time stamped copy of the script
write(paste0("07_dge_1var - Start: ", Sys.time()), file = here::here("scRNA_Log.txt"), append = TRUE)
log_connection <- write_script_log(here::here("R/07_dge_1var.R"), log_dir = here::here(output_dir, "Logs"))

# Setup parallel backend
n_cores <- max(1, parallel::detectCores() - 1)
cl <- makeCluster(n_cores)
registerDoParallel(cl)

# Log all output to the end of the log file
sink(log_connection, append = TRUE)
sink(log_connection, type = "message", append = TRUE)
on.exit({
  sink(NULL)
  sink(NULL, type = "message")
  stopCluster(cl)
})

# Pseudobulking and data processing ####
# Load the clustered Seurat file
combined_seurat <- readRDS(here::here("Data", "R_Data", paste0(scConfig$prefix, "_combined_clustered.rds")))
print(unique(combined_seurat$Sample_name))

# Define the cluster identity for plotting and DGE analysis
cluster_col <- scConfig$cluster_plot_ident

# Ensure cluster_col exists in meta.data, default to seurat_clusters if not
if (!cluster_col %in% colnames(combined_seurat@meta.data)) {
  message("cluster_col not found in metadata; defaulting to 'seurat_clusters'")
  cluster_col <- "seurat_clusters"
}

# Select the assay for DESeq2 analysis
if (scConfig$soupx_adjust) {
  # If SoupX was used, set the default assay to SoupX
  dge_assay <- "SoupX"
} else {
  # If SoupX was not used, use RNA
  dge_assay <- "RNA"
}

# Create pseudobulked seurat object
pseudobulked_seurat <- AggregateExpression(combined_seurat, assays = dge_assay, return.seurat = TRUE,
                                           group.by = c("Treatment", "Sample_name", cluster_col))

pseudobulked_seurat$celltype.treatment <- paste(pseudobulked_seurat$Treatment, pseudobulked_seurat[[cluster_col]][, 1], sep = "_")

Idents(pseudobulked_seurat) <- "celltype.treatment"

colnames_vec <- colnames(pseudobulked_seurat[[dge_assay]])

# Get clusters and conditions
clusters <- unique(pseudobulked_seurat[[cluster_col]][, 1])
conditions <- unique(pseudobulked_seurat$Treatment)

# DGE analysis and plotting ####
summary_list <- foreach(cluster = clusters,
                        .packages = c("Seurat", "DESeq2", "dplyr", "pheatmap", "ggrepel",
                                      "tibble", "stringr", "ggplot2")) %dopar% {
  # Find columns that end with the cluster ID
  cols <- grep(paste0("_", cluster, "$"), colnames_vec, value = TRUE)

  # Subset matrix for this cluster
  mat <- GetAssayData(pseudobulked_seurat, slot = "counts")[, cols, drop = FALSE]

  # Save to CSV
  write.csv(mat, file = here::here(output_dir, "CSV_Results", "DEGs_All", paste0("Pseudobulk_counts_cluster_", cluster, ".csv")))

  # Build colData data.frame from the matrix column names
  sample_names <- colnames(mat)
  Treatment <- sub("_.*", "", sample_names)  # everything before first underscore
  col_data <- data.frame(Treatment = Treatment,
                         row.names = sample_names,
                         Sample_name = sample_names)

  # Create DESeqDataSet and run DESeq2
  dds <- DESeqDataSetFromMatrix(countData = mat, colData = col_data, design = ~Treatment)
  dds <- DESeq(dds)

  # Transform counts for data visualization
  rld <- rlog(dds, blind = TRUE)

  # Plot PCAs by condition and sample
  pca_plot_by_condition <- DESeq2::plotPCA(rld, intgroup = "Treatment", ntop = 50)
  save_plot_pdf(pca_plot_by_condition, here::here(output_dir, "Plots", "DESEQ_Plots", "PCAs", paste0(str_replace(cluster, "_", " "), "_PCA.pdf")),
                height = 6, width = 8)

  pca_plot_by_sample <- DESeq2::plotPCA(rld, intgroup = "Sample_name", ntop = 50)
  save_plot_pdf(pca_plot_by_sample, here::here(output_dir, "Plots", "DESEQ_Plots", "PCAs", paste0(str_replace(cluster, "_", " "), "_sample_PCA.pdf")),
                height = 6, width = 8)

  # Extract the rlog matrix from the object and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)

  # Plot heatmap
  cluster_heatmap <-  pheatmap(rld_cor, annotation = col_data[, c("Treatment"), drop = FALSE])
  save_plot_pdf(cluster_heatmap, here::here(output_dir, "Plots", "DESEQ_Plots", "Heatmaps", paste0(str_replace(cluster, "_", " "), "_heatmap_plot.pdf")),
                height = 6, width = 8)

  # Plot dispersion estimates
  dispersion_plot <- plotDispEsts(dds)
  save_plot_pdf(dispersion_plot, here::here(output_dir, "Plots", "DESEQ_Plots", "Dispersion_Plots", paste0(str_replace(cluster, "_", " "), "_dispersion_plot.pdf")),
                height = 6, width = 8)

  # Set up the contrast for DESeq2, run DGE, and apply LFC shrinkage
  contrast <- c("Treatment", levels(as.factor(col_data$Treatment))[2], levels(as.factor(col_data$Treatment))[1])
  res <- results(dds, contrast = contrast, alpha = 0.05)
  res <- lfcShrink(dds, contrast =  contrast, res = res, type = "normal")

  # Turn the results object into a tibble for use with tidyverse functions
  res_tbl <- res %>%
    data.frame() %>%
    rownames_to_column(var = "gene") %>%
    as_tibble()

  res_tbl %>%
    arrange(desc(padj))

  # Write all results to file
  write.csv(res_tbl,
            here::here(output_dir, "CSV_Results", "DEGs_All", paste0("DEG_Results_All_Genes_", conditions[1], "_vs_", conditions[2], "_Cluster_", cluster, ".csv")),
            quote = FALSE,
            row.names = FALSE)

  # Make a CSV with normalized counts
  normalized_counts <- data.frame(round(counts(dds, normalized = TRUE), 1))
  colnames(normalized_counts) <- paste(colnames(normalized_counts), "norm", sep = "_")
  write.csv(normalized_counts,
            here::here(output_dir, "CSV_Results", "DEGs_All", paste0("Normalized_counts_cluster", cluster, "_all_genes.csv")),
            quote = FALSE,
            row.names = FALSE)

  # Set thresholds
  padj_cutoff <- 0.05
  lfc_cutoff <- 0.25

  # Subset the significant results
  sig_res <- dplyr::filter(res_tbl, padj < padj_cutoff & abs(log2FoldChange) >= lfc_cutoff) %>%
    dplyr::arrange(padj)

  # DEG summary counts
  n_de <- nrow(sig_res)
  n_up <- sum(sig_res$log2FoldChange > 0)
  n_down <- sum(sig_res$log2FoldChange < 0)

  # Make summary
  summary <- data.frame(
    Cluster = cluster,
    Num_DE_Genes = n_de,
    Num_Upregulated = n_up,
    Num_Downregulated = n_down
  )

  # Write significant results to file
  write.csv(sig_res,
            here::here(output_dir, "CSV_Results", "DEGs_All", paste0("DEG_Results_Significant_Genes_", conditions[1], "_vs_", conditions[2], "_Cluster_", cluster, ".csv")),
            quote = FALSE,
            row.names = FALSE)

  # Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 0.25 in either direction
  res_table_thres <- res_tbl %>%
    mutate(threshold = padj < padj_cutoff & abs(log2FoldChange) >= lfc_cutoff)

  # Volcano plot
  volc_plot <- ggplot(res_table_thres, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(aes(colour = threshold)) +
    ggtitle(cluster) +
    xlab("log2 fold change") +
    ylab("-log10 adj p-value") +
    #!scale_y_continuous(limits = c(0,20)) +
    #!scale_x_continuous(limits = c(-2,2)) +
    scale_color_manual(values = c("black", "red")) +
    ggrepel::geom_text_repel(data = head(sig_res, 10), aes(x = log2FoldChange, y = -log10(padj), label = gene)) +
    theme(legend.position = "none",
          plot.title = element_text(size = rel(1.5), hjust = 0.5),
          axis.title = element_text(size = rel(1.25)))

  # Save the volcano plot
  save_plot_pdf(volc_plot, here::here(output_dir, "Plots", "DESEQ_Plots", "Volcano_Plots", paste0(str_replace(cluster, "_", " "), "_volcano_plot.pdf")),
                height = 6, width = 8)

  # Make an MA plot for quality control
  ma_plot <- plotMA(res, ylim = c(-3,3), alpha = 0.05, main = (paste0("Cluster: ", as.character(cluster))))
  save_plot_pdf(ma_plot, here::here(output_dir, "Plots", "DESEQ_Plots", "MA_Plots", paste0(str_replace(cluster, "_", " "), "_MA_plot.pdf")),
                height = 6, width = 8)

  # Return the summary
  summary
}

# Combine and write summary
summary_df <- dplyr::bind_rows(summary_list)
write.csv(summary_df, here::here(output_dir, "CSV_Results", "DGE_summary.csv"), row.names = FALSE)

write(paste0("07_dge_1var - Finish: ", Sys.time()), file = here::here("scRNA_Log.txt"), append = TRUE)
