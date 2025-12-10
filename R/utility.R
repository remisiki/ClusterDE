#' @useDynLib ClusterDE, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#' Log-transformed Counts per 10,000 (logCP10K)
#'
#' Computes the per-cell log-transformed counts per 10,000 for a matrix of raw counts.
#'
#' @param x A matrix or \code{Matrix::dgCMatrix} of raw counts with genes as rows and cells as columns.
#' @param scale.factor A numeric value representing the scaling factor for total counts per cell (default: 1e4).
#'
#' @return A matrix of the same dimensions as \code{x}, containing log-transformed normalized counts per cell.
logcp10k <- function(x, scale.factor = 1e4) {
  # L_j
  lib <- Matrix::colSums(x)
  # per-cell scaling; 0 for empty columns
  sf <- ifelse(lib > 0, scale.factor / lib, 0)
  # cp10k = C * (scale.factor / L_j)
  X <- t(t(x) * sf)
  # log(1 + cp10k)
  log1p(X)
}
