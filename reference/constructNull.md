# Construct the synthetic null data

`constructNull` takes the target data as the input and returns the
corresponding synthetic null data.

## Usage

``` r
constructNull(
  obj,
  family = "nb",
  formula = NULL,
  extraInfo = NULL,
  nCores = 1,
  nRep = 1,
  parallelization = "mcmapply",
  fastVersion = TRUE,
  ifSparse = FALSE,
  corrCut = 0.1,
  BPPARAM = NULL,
  approximation = FALSE,
  usePca = F,
  nPcs = 200
)
```

## Arguments

  - obj:
    
    A Seurat object. The reference data.

  - family:
    
    A string or a vector of strings of the distribution of your data.
    Must be one of 'nb', 'binomial', 'poisson', 'zip', 'zinb' or
    'gaussian', which represent 'poisson distribution', 'negative
    binomial distribution', 'zero-inflated poisson distribution',
    'zero-inflated negative binomail distribution', and 'gaussian
    distribution' respectively. For UMI-counts data, we usually use
    'nb'. Default is 'nb'.

  - formula:
    
    A string of the mu parameter formula. It defines the relationship
    between gene expression in synthetic null data and the extra
    covariates. Default is NULL (cell type case). For example, if your
    input data is a spatial data with X, Y coordinates, the formula can
    be 's(X, Y, bs = 'gp', k = 4)'.

  - extraInfo:
    
    A data frame of the extra covariates used in `formula`. For example,
    the 2D spatial coordinates. Default is NULL.

  - nCores:
    
    An integer. The number of cores to use for Parallel processing.

  - nRep:
    
    An integer. The number of sampled synthetic null datasets. Default
    value is 1.

  - parallelization:
    
    A string indicating the specific parallelization function to use.
    Must be one of 'mcmapply', 'bpmapply', or 'pbmcmapply', which
    corresponds to the parallelization function in the package
    `parallel`,`BiocParallel`, and `pbmcapply` respectively. The default
    value is 'pbmcmapply'.

  - fastVersion:
    
    A logic value. If TRUE, the fast approximation is used. Default is
    FALSE.

  - ifSparse:
    
    A logic value. For high-dimensional data (gene number is much larger
    than cell number), if a sparse correlation estimation will be used.
    Default is FALSE.

  - corrCut:
    
    A numeric value. The cutoff for non-zero proportions in genes used
    in modelling correlation.

  - BPPARAM:
    
    A `MulticoreParam` object or NULL. When the parameter
    parallelization = 'mcmapply' or 'pbmcmapply', this parameter must be
    NULL. When the parameter parallelization = 'bpmapply', this
    parameter must be one of the

  - usePca:
    
    A logic value. Whether to use PCA approximation. Default is FALSE.

  - nPcs:
    
    A numeric value. Number of PCs to use when usePca=T. Default is 200.
    `MulticoreParam` object offered by the package 'BiocParallel. The
    default value is NULL.

  - Approximation:
    
    A logic value. For a high-latitude counting matrix, Approximation
    can increase the speed of data generation while ensuring accuracy.
    Note that it only takes effect if "fastVersion=TRUE,
    Approximation=TRUE". Default is FALSE.

## Value

The expression matrix of the synthetic null data.

## Details

This function constructs the synthetic null data based on the target
data (real data). The input is a expression matrix (gene by cell); the
user should specify a distribution, which is usually Negative Binomial
for count matrix.

## Examples

``` r
data(exampleCounts)
nullData <- constructNull(mat = exampleCounts)
#> Error in constructNull(mat = exampleCounts): unused argument (mat = exampleCounts)
```
