# Plot reduced-dimensional embeddings for reference and null data

Generates a UMAP (or other reduced dimension) plot comparing a reference
dataset and one or multiple null datasets. Input objects can be either
`Seurat` or `SingleCellExperiment` objects. All objects are internally
converted to `SingleCellExperiment`.

## Usage

``` r
dimPlot(
  ref_obj,
  null_obj,
  color_by = "seurat_clusters",
  seed = 0,
  nrow = NULL,
  ncol = NULL
)
```

## Arguments

  - ref\_obj:
    
    A `Seurat` or `SingleCellExperiment` object representing the
    reference dataset.

  - null\_obj:
    
    A single object or a list of objects, each being either `Seurat` or
    `SingleCellExperiment`, representing the null datasets.

  - color\_by:
    
    A character string specifying the cell metadata column used to color
    the embedding. Default: `"seurat_clusters"`.

  - seed:
    
    Integer random seed for reproducibility. Default: `0`.

  - nrow:
    
    Number of rows for faceting multiple datasets. Default: `NULL`.

  - ncol:
    
    Number of columns for faceting multiple datasets. Default: `NULL`.

## Value

A ggplot object showing the embedded reference and null datasets.
