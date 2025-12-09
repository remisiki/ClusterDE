# Calculate pval for each gene in synthetic null data

Runs a Seurat-based null differential expression analysis for provided
synthetic null datasets.

## Usage

``` r
calcNullPval(
  synthetic_null,
  spatial_coords = NULL,
  normalize = T,
  hvg = NULL,
  seed = 123,
  nCores = 1
)
```

## Arguments

  - synthetic\_null:
    
    A matrix (genes x cells), or a list of such matrices representing
    synthetic null data. If not a list, it will be coerced to a list.

  - spatial\_coords:
    
    A data frame, should contain two columns representing X and Y
    coordinates if using spatial data. Default is NULL.

  - normalize:
    
    Logical; if `TRUE` (default), apply Seurat's NormalizeData to each
    dataset. If `FALSE`, sets assay data directly. Set to false if null
    data is generated using PCA approximation. It is recommended to
    provide hvg if normalize is false.

  - hvg:
    
    A list of feaures; will be applied if normalization is set to false.

  - seed:
    
    Numeric; random seed for clustering reproducibility (default: 123).

  - nCores:
    
    An integer. The number of cores to use for parallel processing.

## Value

A named list with preprocess data and gene p-values
