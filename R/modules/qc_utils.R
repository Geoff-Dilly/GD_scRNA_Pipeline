# qc_utils.R
# Author: Geoff Dilly

library(here)
library(Seurat)

#' @title QC Filtering
#' @description Flters a seurat object base on scConfig
qc_filter_seurat <- function(seurat_obj, config = scConfig) {
  # Set soupx to default assay if needed
  if (config$soupx_adjust == TRUE) {
    DefaultAssay(seurat_obj) <- "SoupX"
  }

  # Remove cells that don't meet QC cutoffs
  required_cols <- c("nFeature_RNA", "percent_mito", "percent_ribo")
  missing_cols <- setdiff(required_cols, colnames(seurat_obj@meta.data))
  if (length(missing_cols) > 0) {
    stop("Missing columns in meta.data: ", paste(missing_cols, collapse = ", "))
  }
  seurat_obj <- subset(
    seurat_obj,
    subset = nFeature_RNA > config$nFeature_RNA_cutoff &
      percent_mito < config$percent_mito_cutoff &
      percent_ribo < config$percent_ribo_cutoff
  )

  # Remove the top quartile high UMI cells if specified
  if (config$remove_top_nUMIs == TRUE) {
    umi_threshold <- quantile(seurat_obj$nCount_RNA, 0.75)
    seurat_obj <- subset(seurat_obj, subset = nCount_RNA < umi_threshold)
  }

  # Remove called doublets if specified
  if (config$remove_doublets == TRUE) {
    seurat_obj <- subset(seurat_obj, subset = Doublet_Call == "Singlet")
  }

  # Remove mitochondrial genes if specified
  if (config$remove_mito_genes == TRUE) {
    mito_genes <- grep(config$mito_pattern, rownames(seurat_obj), value = TRUE)
    seurat_obj <- subset(
      seurat_obj,
      features = setdiff(rownames(seurat_obj), mito_genes)
    )
  }

  # Remove ribosomal genes if specified
  if (config$remove_ribo_genes == TRUE) {
    ribo_genes <- grep(config$ribo_pattern, rownames(seurat_obj), value = TRUE)
    seurat_obj <- subset(
      seurat_obj,
      features = setdiff(rownames(seurat_obj), ribo_genes)
    )
  }

  return(seurat_obj)
}

#' @title Strip Seurat Object
#' @description Removes unused data from a Seurat object
strip_seurat_obj <- function(seurat_obj) {
  # Strip reductions and related data
  seurat_obj@reductions <- list()
  VariableFeatures(seurat_obj) <- character(0)
  seurat_obj@graphs <- list()
  seurat_obj@neighbors <- list()
  seurat_obj@commands <- list()
  seurat_obj@misc <- list()

  # Strip Normalized and scaled data
  seurat_obj[["RNA"]]$scale.data <- NULL
  seurat_obj[["RNA"]]$data <- NULL
  seurat_obj[["SCT"]] <- NULL

  return(seurat_obj)
}

#' @title Run SoupX Correction
#' @description Run SoupX on an unfiltered Seurat object and add as assay
run_soupx_correction <- function(sample, seurat_obj, output_dir = "R_Data") {
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
  saveRDS(soup_channel, file = file.path(output_dir, paste0(sample$Sample_name, "_SoupChannel.rds")))

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
#' @title Run DoubletFinder
#' @description Run DoubletFinder to detect doublets in a Seurat object
run_doubletfinder <- function(sample, seurat_obj, config = scConfig) {
  # Check if DoubletFinder is installedd
  if (!requireNamespace("DoubletFinder", quietly = TRUE)) stop("DoubletFinder package is required.")

  # Normalize and run UMAP
  seurat_obj <- Seurat::SCTransform(seurat_obj, verbose = FALSE)
  seurat_obj <- Seurat::RunPCA(seurat_obj)
  seurat_obj <- Seurat::RunUMAP(seurat_obj, dims = 1:10)

  # Sweep to find the optimal pK value for DoubletFinder
  sweep_res_list <- DoubletFinder::paramSweep(seurat_obj, PCs = 1:10, sct = FALSE)
  sweep_stats <- DoubletFinder::summarizeSweep(sweep_res_list, GT = FALSE)
  bcmvn <- DoubletFinder::find.pK(sweep_stats)
  pK_value <- as.numeric(as.vector(bcmvn$pK[which.max(bcmvn$BCmetric), drop = TRUE]))

  homotypic_prop <- DoubletFinder::modelHomotypic(seurat_obj@meta.data$seurat_clusters)
  nExp_poi <- round((config$expct_doublet_pct/100) * nrow(seurat_obj@meta.data))
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic_prop))

  # Run DoubletFinder
  seurat_obj <- DoubletFinder::doubletFinder(seurat_obj,
                                             PCs = 1:10,
                                             pN = 0.25,
                                             pK = pK_value,
                                             nExp = nExp_poi.adj,
                                             sct = TRUE)

  # Standardize and save DoubletFinder results
  meta_cols <- colnames(seurat_obj@meta.data)
  score <- stringr::str_subset(meta_cols, "^pANN")
  call <- stringr::str_subset(meta_cols, "^DF.cl")
  seurat_obj$Doublet_Score <- seurat_obj[[score]]
  seurat_obj$Doublet_Call <- seurat_obj[[call]]

  # Strip unnecessary data from the seurat object
  seurat_obj <- strip_seurat_obj(seurat_obj)

  saveRDS(
    seurat_obj,
    file = here::here("R_Data", paste0(sample, "_seurat_Doublets.rds"))
  )
  return(list(sample, score, call))
}