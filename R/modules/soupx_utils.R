# soupx_utils.R

#' @title Run SoupX Correction
#' @description Run SoupX on an unfiltered Seurat object and add as assay
run_soupx_correction <- function(sample, seurat_obj) {
  require(Seurat)
  require(SoupX)

  # Load the unfiltered matrix for SoupX
  filt_matrix <- Read10X(file.path(sample$Raw_data_dir, "filtered_feature_bc_matrix"))
  raw_matrix <- Read10X(file.path(sample$Raw_data_dir, "raw_feature_bc_matrix"))

  # Minimally cluster the seurat seurat_object
  seurat_obj <- SCTransform(seurat_obj, verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj, verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
  seurat_obj <- FindClusters(seurat_obj)

  # Assign clusters to SoupX
  cluster_labels <- seurat_obj$seurat_clusters
  cell_names     <- rownames(seurat_obj@meta.data)

  # Run SoupX
  soup_channel <- SoupChannel(tod = raw_matrix, toc = filt_matrix)
  soup_channel <- setClusters(soup_channel, setNames(cluster_labels, cell_names))
  soup_channel <- autoEstCont(soup_channel)
  adj_matrix   <- adjustCounts(soup_channel, roundToInt = TRUE)

  # In your SoupX correction step (after autoEstCont and adjustCounts)
  saveRDS(soup_channel, file = file.path("R_Data", paste0(sample$Sample_name, "_SoupChannel.rds")))

  # Remake a clean Seurat object and add SoupX results
  out_obj <- CreateSeuratObject(counts = filt_matrix)
  out_obj[["SoupX"]] <- CreateAssayObject(counts = adj_matrix)

  # Save the top ambient genes
  top_ambient <- head(soup_channel$soupProfile[order(soup_channel$soupProfile$est, decreasing = TRUE), ], n = 30)
  top_ambient <- tibble::rownames_to_column(top_ambient, "gene")
  top_ambient$Sample <- sample$Sample_name

  # Return as a list for easy unpacking
  return(list(seurat_obj = out_obj, top_ambient = top_ambient))
}