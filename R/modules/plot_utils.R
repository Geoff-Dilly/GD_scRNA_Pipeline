# plot_utils.R

library(Seurat)
library(patchwork)
library(ggplot2)
library(purrr)

#' @title Save ggplot to PDF
#' @description Saves a ggplot object to a PDF file with specified dimensions.
save_plot_pdf <- function(plot, filename, height = 4, width = 6) {
  pdf(filename, height = height, width = width)
  print(plot)
  dev.off()
}

#' @title Save ggplot to PNG
#' @description Saves a ggplot object to a PNG file with specified dimensions.
save_plot_png <- function(plot, filename, height = 4, width = 6, dpi = 300) {
  ggsave(filename, plot, height = height, width = width, dpi = dpi)
}

#' @title Make stacked violin plot
#' @description Create a stacked violin plot for a set of features (ultrastack style)
#'
#' @param seurat_obj      Seurat object
#' @param features        Character vector of features (genes) to plot
#' @param assay           Assay to use (default "SCT")
#' @param layer           Layer to use (default "data")
#' @param adjust          Density adjust parameter (default 2)
#' @param pt.size         Point size (default 0)
#' @param gene_label_size Strip text (gene) label size (default 8)
#' @param ...             Additional arguments for VlnPlot
#' @return                A ggplot/patchwork object
make_stacked_vln_plot <- function(
  seurat_obj,
  features,
  assay = "SCT",
  layer = "data",
  adjust = 2,
  pt.size = 0,
  gene_label_size = 8,
  ...
) {
  # Calculate global max for y-axis
  df <- FetchData(seurat_obj, vars = features, layer = layer, assay = assay)
  global_max <- ceiling(max(df, na.rm = TRUE))

  # Make the plot
  p <- VlnPlot(
    seurat_obj,
    features = features,
    pt.size = pt.size,
    stack = TRUE,
    flip = TRUE,
    assay = assay,
    layer = layer,
    adjust = adjust,
    ...
  ) +
    scale_y_continuous(breaks = global_max) +
    theme(
      legend.position = "none",
      axis.ticks.y = element_blank(),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 8),
      strip.text.y = element_text(size = gene_label_size)
    )
  
  return(p)
}

#' @title Plot aligned dotplots by cell type from marker gene CSV (grid output)
#' @description Returns a single patchwork grid of all cell type dotplots, well aligned.
plot_dotplots_by_celltype_grid <- function(
    seurat_obj,
    marker_csv = "reference/marker_genes_db.csv",
    reference_name,
    group_by = "seurat_clusters",
    min_genes = 1,
    ncol = 3, # Number of columns in grid (change as needed)
    height = 6, width = 8 # Default plot size if saving
) {
  require(ggplot2)
  require(dplyr)
  require(readr)
  require(patchwork)
  require(Seurat)

  # Read marker table
  marker_tbl <- read_csv(marker_csv)

  # Filter by reference
  filtered_tbl <- marker_tbl %>% filter(reference == reference_name)

  # Get all genes and cell_types (in order)
  all_genes <- unique(filtered_tbl$gene)
  all_cell_types <- unique(filtered_tbl$cell_type)

  plot_list <- list()
  for (ct in all_cell_types) {
    genes <- filtered_tbl %>% filter(cell_type == ct) %>% pull(gene)
    genes <- genes[!is.na(genes)]
    if (length(genes) >= min_genes) {
      ordered_genes <- all_genes[all_genes %in% genes]
      plot_genes <- all_genes

      p <- DotPlot(seurat_obj, features = plot_genes, group.by = group_by) +
        RotatedAxis() +
        ggtitle(paste("Cell Type:", ct, "| Reference:", reference_name)) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
          axis.text.y = element_text(size = 8),
          plot.title = element_text(face = "bold", size = 12)
        ) +
        scale_y_discrete(limits = rev(plot_genes)) +
        geom_point(
          data = function(data) {
            data$highlight <- data$features.plot %in% ordered_genes
            data
          },
          aes(alpha = highlight)
        ) +
        scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.2), guide = "none")
      plot_list[[ct]] <- p
    }
  }
  # Combine with patchwork
  grid_plot <- wrap_plots(plot_list, ncol = ncol)
  return(grid_plot)
}