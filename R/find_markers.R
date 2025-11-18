#' Find differential markers using ClusterDE
#'
#' @param obj Seurat object containing expression data.
#' @param ident.1 Group identity to test as group 1.
#' @param ident.2 Group identity to test as group 2.
#' @param nCores Number of CPU cores used for parallelization.
#' @param nRep Number of null replicates.
#' @param flavour Either `"classic"` or `"pca"`; selects null-construction method. **classic**: uses default count matrix as simulation input. **pca**: uses PCA approximation as simulation input, faster and more scalable when large number of replicates is needed.
#'
#' @return A table of DEG results.
#'
#' @examples
#' dummy_obj <- DUMMY
#' findMarkers(dummy_obj, ident.1 = "cell_type_1", ident.2 = "cell_type_2")
#' @export
findMarkers <- function(
  obj,
  ident.1,
  ident.2,
  nCores = 1,
  nRep = 1,
  flavour = "classic"
) {
  original_markers <- Seurat::FindMarkers(
    obj,
    ident.1 = ident.1,
    ident.2 = ident.2,
    min.pct = 0,
    logfc.threshold = 0
  )
  original_pval <- original_markers$p_val
  names(original_pval) <- rownames(original_markers)
  if (flavour == "classic") {
    null_data <- constructNull(obj, nCores = nCores, nRep = nRep, fastVersion = T)
    null_pval <- calcNullPval(null_data)
    callDE(original_pval, null_pval)
  } else if (flavour == "pca") {
    null_data <- constructNull(obj, usePca = T, nCores = nCores, nRep = nRep)
    null_pval <- calcNullPval(null_data, normalize = F, hvg = Seurat::VariableFeatures(obj), nCores = nCores)
    callDE(original_pval, null_pval, nCores = nCores)
  }
}
