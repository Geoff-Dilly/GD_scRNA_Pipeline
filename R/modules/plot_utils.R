# plot_utils.R

library(Seurat)
library(patchwork)
library(ggplot2)
library(purrr)
library(dplyr)
library(readr)

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

dotplot_by_marker_group <- function(
  seurat_obj,
  marker_csv = "reference/marker_gene_db.csv",
  reference_name = "Dilly_et_al_2022",
  group_by = "seurat_clusters",
  gene_group_col = "cell_type",   # or "cell_class", "function", etc.
  add_separators = TRUE,
  separator_color = "grey70",
  separator_linetype = "dashed",
  separator_size = 0.7
) {

  # Read marker gene CSV
  marker_tbl <- read_csv(marker_csv, show_col_types = FALSE) %>%
    filter(reference == reference_name)

  # Prepare ordered gene list, grouped by category
  marker_tbl <- marker_tbl %>%
    mutate(!!gene_group_col := factor(.data[[gene_group_col]], levels = unique(.data[[gene_group_col]]))) %>%
    arrange(.data[[gene_group_col]], gene)
  gene_ordered <- marker_tbl$gene %>% unique()

  # Make DotPlot
  p <- DotPlot(seurat_obj, features = gene_ordered, group.by = group_by) +
    RotatedAxis() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = element_text(size = 9),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(face = "bold", size = 14)
    ) +
    ggtitle(paste0("DotPlot: Marker Genes Grouped by ", gene_group_col))

  # Add vertical separators between groups (optional)
  if (add_separators) {
    group_lengths <- table(marker_tbl[[gene_group_col]])
    breaks <- cumsum(group_lengths)
    # Remove last break to avoid line at the end
    if (length(breaks) > 1) {
      for (b in breaks[-length(breaks)]) {
        p <- p + geom_vline(xintercept = b + 0.5, linetype = separator_linetype,
                            color = separator_color, size = separator_size)
      }
    }
  }

  return(p)
}
