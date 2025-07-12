# 07_dge_1var.R
# Purpose: Perform DESEQ2 analysis on pseudobulked Seurat objects and make plots
# Author: Geoff Dilly

library(here)
library(foreach)
library(doParallel)
library(Seurat)
library(tidyverse)
library(DESeq2)
library(ggrepel)
library(pheatmap)
scRNA_home_dir <- here()
setwd(scRNA_home_dir)

# Setup ####
# Load custom functions
source("R/modules/plot_utils.R")
source("R/modules/log_utils.R")
write_script_log("R/01b_soupx.R")

# Log the start time and a time stamped copy of the script
write(paste0("07_dge_1var - Start: ", Sys.time()), file = "scRNA_Log.txt", append = TRUE)
write_script_log("R/07_dge_1var.R")

# Load the configuration file and metadata
source("sc_experiment_config.R")
scConfig.Sample_metadata <- read.csv("sc_sample_metadata.csv")

# Setup parallel backend
n_cores <- parallel::detectCores() - 1
cl <- makeCluster(n_cores/2)
registerDoParallel(cl)

# Pseudobulking and data processing ####
# Load the clustered Seurat file
combined_seurat <- readRDS(paste0("R_Data/", scConfig.Prefix, "_combined_clustered.rds"))
print(unique(combined_seurat$Sample_name))

# Select the assay for DESeq2 analysis
if (scConfig.soupx_adjust == TRUE) {
  # If SoupX was used, set the default assay to SoupX
  dge_assay <- "SoupX"
} else {
  # If SoupX was not used, use RNA
  dge_assay <- "RNA"
}

# Create pseudobulked seurat object
pseudobulked_seurat <- AggregateExpression(combined_seurat, assays = dge_assay, return.seurat = TRUE,
                                           group.by = c("Treatment", "Sample_name", "seurat_clusters"))

pseudobulked_seurat$celltype.treatment <- paste(pseudobulked_seurat$Treatment, pseudobulked_seurat$seurat_clusters, sep = "_")

Idents(pseudobulked_seurat) <- "celltype.treatment"

colnames_vec <- colnames(pseudobulked_seurat$RNA)

# Get clusters and conditions
clusters <- unique(pseudobulked_seurat$seurat_clusters)
conditions <- unique(pseudobulked_seurat$Treatment)

# DGE Analysis and Plotting ####
summary_list <- foreach(cluster = clusters,
                        .packages = c("Seurat", "DESeq2", "dplyr", "pheatmap", "ggrepel",
                                      "tibble", "stringr", "ggplot2")) %dopar% {
  # Find columns that end with the cluster ID
  cols <- grep(paste0("_", cluster, "$"), colnames_vec, value = TRUE)

  # Subset matrix for this cluster
  mat <- GetAssayData(pseudobulked_seurat, slot = "counts")[, cols, drop = FALSE]

  # Save to CSV (change file path as needed)
  write.csv(mat, file = paste0("CSV_Results/DEGs_All/", "Pseudobulk_counts_cluster_", cluster, ".csv"))

  # Build colData data.frame from the matrix column names
  sample_names <- colnames(mat)

  # Extract treatment from the names (assumes format: Treatment_Sample-#_Cluster#)
  Treatment <- sub("_.*", "", sample_names)  # everything before first underscore # nolint
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
  save_plot_pdf(pca_plot_by_condition, paste0("Plots/DESEQ_Plots/PCAs/", str_replace(cluster, "_", " "), "_PCA.pdf"),
                height = 6, width = 8)

  pca_plot_by_sample <- DESeq2::plotPCA(rld, intgroup = "Sample_name", ntop = 50)
  save_plot_pdf(pca_plot_by_sample, paste0("Plots/DESEQ_Plots/PCAs/", str_replace(cluster, "_", " "), "_sample_PCA.pdf"),
                height = 6, width = 8)

  # Extract the rlog matrix from the object and compute pairwise correlation values
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)

  # Plot heatmap
  cluster_heatmap <-  pheatmap(rld_cor, annotation = col_data[, c("Treatment"), drop = FALSE])
  save_plot_pdf(cluster_heatmap, paste0("Plots/DESEQ_Plots/Heatmaps/", str_replace(cluster, "_", " "), "_heatmap_plot.pdf"), 
                height = 6, width = 8)

  # Plot dispersion estimates
  dispersion_plot <- plotDispEsts(dds)
  save_plot_pdf(dispersion_plot, paste0("Plots/DESEQ_Plots/Dispersion_Plots/", str_replace(cluster, "_", " "), "_dispersion_plot.pdf"), 
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
            paste0("CSV_Results/DEGs_All/DEG_Results_All_Genes_", conditions[1], "_vs_", conditions[2], "_Cluster_", cluster, ".csv"),
            quote = FALSE,
            row.names = FALSE)

  # Make a CSV with normalized counts
  normalized_counts <- data.frame(round(counts(dds, normalized = TRUE), 1))
  colnames(normalized_counts) <- paste(colnames(normalized_counts), "norm", sep = "_")
  write.csv(normalized_counts,
            paste0("CSV_Results/DEGs_All/", "Normalized_counts_cluster", cluster, "_all_genes.csv"),
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
            paste0("CSV_Results/DEGs_All/", "DEG_Results_Significant_Genes_", conditions[1], "_vs_", conditions[2], "_Cluster_", cluster, ".csv"),
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
  save_plot_pdf(volc_plot, paste0("Plots/DESEQ_Plots/Volcano_Plots/", str_replace(cluster, "_", " "), "_volcano_plot.pdf"), 
                height = 6, width = 8)

  # Make an MA plot for quality control
  ma_plot <- plotMA(res, ylim = c(-3,3), alpha = 0.05, main = (paste0("Cluster: ", as.character(cluster)))) # nolint
  save_plot_pdf(ma_plot, paste0("Plots/DESEQ_Plots/MA_Plots/", str_replace(cluster, "_", " "), "_MA_plot.pdf"), 
                height = 6, width = 8)

  # Return the summary
  summary
}

# Stop the parallel cluster
stopCluster(cl)

# Combine and write summary
summary_df <- dplyr::bind_rows(summary_list)
write.csv(summary_df, "CSV_Results/DEGs_All/DGE_summary.csv", row.names = FALSE)

write(paste0("07_dge_1var - Finish: ", Sys.time()), file = "scRNA_Log.txt", append = TRUE)
