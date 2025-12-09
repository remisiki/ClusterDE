# Find differential markers using ClusterDE

Find differential markers using ClusterDE

## Usage

``` r
findMarkers(
  obj,
  ident.1,
  ident.2,
  nCores = 1,
  nRep = 1,
  flavour = "classic",
  spatial = NULL
)
```

## Arguments

  - obj:
    
    Seurat object containing expression data.

  - ident.1:
    
    Group identity to test as group 1.

  - ident.2:
    
    Group identity to test as group 2.

  - nCores:
    
    Number of CPU cores used for parallelization.

  - nRep:
    
    Number of null replicates.

  - flavour:
    
    Either `"classic"` or `"pca"`; selects null-construction method.
    **classic**: uses default count matrix as simulation input. **pca**:
    uses PCA approximation as simulation input, faster and more scalable
    when large number of replicates is needed.

  - spatial:
    
    A vector of 2 strings, the meta data column name representing X and
    Y coordinates if using spatial data. Default is NULL.

## Value

A table of DEG results.

## Examples

``` r
dummy_obj <- DUMMY
#> Error: object 'DUMMY' not found
findMarkers(dummy_obj, ident.1 = "cell_type_1", ident.2 = "cell_type_2")
#> Error: object 'dummy_obj' not found
```
