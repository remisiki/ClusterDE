# Log-transformed Counts per 10,000 (logCP10K)

Computes the per-cell log-transformed counts per 10,000 for a matrix of
raw counts.

## Usage

``` r
logcp10k(x, scale.factor = 10000)
```

## Arguments

  - x:
    
    A matrix or `Matrix::dgCMatrix` of raw counts with genes as rows and
    cells as columns.

  - scale.factor:
    
    A numeric value representing the scaling factor for total counts per
    cell (default: 1e4).

## Value

A matrix of the same dimensions as `x`, containing log-transformed
normalized counts per cell.
