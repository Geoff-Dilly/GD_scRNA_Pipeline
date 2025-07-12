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
#' @param slot            Slot to use (default "data")
#' @param adjust          Density adjust parameter (default 2)
#' @param pt.size         Point size (default 0)
#' @param gene_label_size Strip text (gene) label size (default 8)
#' @param ...             Additional arguments for VlnPlot
#' @return                A ggplot/patchwork object
make_stacked_vln_plot <- function(
  seurat_obj,
  features,
  assay = "SCT",
  slot = "data",
  adjust = 2,
  pt.size = 0,
  gene_label_size = 8,
  ...
) {
  # Calculate global max for y-axis
  df <- FetchData(seurat_obj, vars = features, slot = slot, assay = assay)
  global_max <- ceiling(max(df, na.rm = TRUE))

  # Make the plot
  p <- VlnPlot(
    seurat_obj,
    features = features,
    pt.size = pt.size,
    stack = TRUE,
    flip = TRUE,
    assay = assay,
    slot = slot,
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