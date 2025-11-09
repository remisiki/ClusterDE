#' Find Differentially Expressed Genes (DEGs) in Synthetic Null Data
#'
#' Runs a Seurat-based null differential expression analysis for provided synthetic null datasets.
#'
#' @param synthetic_null A matrix (genes x cells), or a list of such matrices representing synthetic null data. If not a list, it will be coerced to a list.
#' @param normalize Logical; if \code{TRUE} (default), apply Seurat's NormalizeData to each dataset. If \code{FALSE}, sets assay data directly. Set to false if null data is generated using PCA approximation. It is recommended to provide hvg if normalize is false.
#' @param hvg A list of feaures; will be applied if normalization is set to false.
#' @param seed Numeric; random seed for clustering reproducibility (default: 123).
#' @param nCores An integer. The number of cores to use for Parallel processing.
#'
#' @return A vector (if single input) or list (if multiple) of named p-values for each gene, corresponding to the null DEGs found by Seurat.
#'
#' @export
findNullDeg <- function(synthetic_null, normalize = T, hvg = NULL, seed = 123, nCores = 1) {
  if (!is.list(synthetic_null)) {
    synthetic_null <- list(synthetic_null)
  }
  null_pval_list <- bettermc::mclapply(synthetic_null, function(count_mat) {
    set.seed(seed)
    if (normalize) {
      data <- Seurat::CreateSeuratObject(count_mat)
      data <- Seurat::NormalizeData(data)
      data <- Seurat::FindVariableFeatures(data)
    } else {
      zero_mat <- Matrix::Matrix(0, nrow(count_mat), ncol(count_mat), sparse = T)
      rownames(zero_mat) <- rownames(count_mat)
      colnames(zero_mat) <- colnames(count_mat)
      data <- Seurat::CreateSeuratObject(zero_mat)
      data <- Seurat::SetAssayData(data, layer = "data", new.data = count_mat)
      Seurat::VariableFeatures(data) <- if (!is.null(hvg)) {
        hvg
      } else {
        warning("hvg not given, use all features")
        rownames(data)
      }
    }
    data <- Seurat::ScaleData(data)
    data <- Seurat::RunPCA(data)
    data <- Seurat::FindNeighbors(data)

    #### find two clusters in the null data ####
    right <- 0.3
    data <- Seurat::FindClusters(data, resolution = right)
    number_of_clusters <- length(unique(data$seurat_clusters))
    while (number_of_clusters < 2) {
      right <- right * 2
      data <- Seurat::FindClusters(data, resolution = right)
      number_of_clusters <- length(unique(data$seurat_clusters))
    }
    left <- 0
    while (number_of_clusters != 2) {
      mid <- (left + right) / 2
      data <- Seurat::FindClusters(data, resolution = mid)
      number_of_clusters <- length(unique(data$seurat_clusters))
      if (number_of_clusters < 2) {
        left <- mid
      } else {
        right <- mid
      }
    }
    data <- Seurat::FindClusters(data, resolution = right)
    deg <- Seurat::FindMarkers(
      data,
      ident.1 = "0",
      ident.2 = "1",
      min.pct = 0,
      logfc.threshold = 0
    )
    null_pval <- deg$p_val
    names(null_pval) <- rownames(deg)
    null_pval
  }, mc.cores = nCores, mc.retry = 5)
  if (length(null_pval_list) == 1) {
    null_pval_list[[1]]
  } else {
    null_pval_list
  }
}
