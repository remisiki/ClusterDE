#' Plot reduced-dimensional embeddings for reference and null data
#'
#' @description
#' Generates a UMAP (or other reduced dimension) plot comparing a reference
#' dataset and one or multiple null datasets. Input objects can be either
#' \code{Seurat} or \code{SingleCellExperiment} objects. All objects are
#' internally converted to \code{SingleCellExperiment}.
#'
#' @param ref_obj A \code{Seurat} or \code{SingleCellExperiment} object
#'   representing the reference dataset.
#' @param null_obj A single object or a list of objects, each being either
#'   \code{Seurat} or \code{SingleCellExperiment}, representing the null datasets.
#' @param color_by A character string specifying the cell metadata column used
#'   to color the embedding. Default: \code{"seurat_clusters"}.
#' @param seed Integer random seed for reproducibility. Default: \code{0}.
#' @param nrow Number of rows for faceting multiple datasets. Default: \code{NULL}.
#' @param ncol Number of columns for faceting multiple datasets. Default: \code{NULL}.
#'
#' @return A ggplot object showing the embedded reference and null datasets.
#' @export
dimPlot <- function(
  ref_obj,
  null_obj,
  color_by = "seurat_clusters",
  seed = 0,
  nrow = NULL,
  ncol = NULL
) {
  if (inherits(ref_obj, "Seurat")) {
    ref_obj <- Seurat::as.SingleCellExperiment(ref_obj)
  } else if (!inherits(ref_obj, "SingleCellExperiment")) {
    stop("Non supported ref_obj class")
  }
  if (!is.list(null_obj)) {
    null_obj <- list(null_obj)
  }
  null_obj <- lapply(null_obj, function(obj) {
    if (inherits(obj, "Seurat")) {
      Seurat::as.SingleCellExperiment(obj)
    } else if (!inherits(obj, "SinglCellExperiment")) {
      stop("Non supported null_obj class")
    } else {
      obj
    }
  })
  set.seed(seed)
  null_names <- paste0("Null", seq_along(null_obj))
  fig <- scDesign3::plot_reduceddim(ref_obj, null_obj, c("Ref", null_names), color_by = color_by)$p_umap +
    ggplot2::facet_wrap(~ Method, nrow = nrow, ncol = ncol) +
    ggplot2::geom_point(alpha = 0.1) +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "bottom", aspect.ratio = 1)
}