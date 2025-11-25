# Call differentially expressed (DE) genes by Clipper algorithm

`callDE` takes two vectors representing each gene's significance. They
are usually p-values from the target data and synthetic null data.

## Usage

``` r
callDE(
  targetScores,
  nullScores,
  nlogTrans = TRUE,
  FDR = 0.05,
  contrastScore = "diff",
  correct = FALSE,
  threshold = "BC",
  ordering = TRUE,
  nCores = 1
)
```

## Arguments

  - targetScores:
    
    A named numeric vector of the DE scores from the target data, e.g.,
    the p-values between two clusters from the real data.

  - nullScores:
    
    A named numeric vector of the DE scores from the synthetic null
    data, e.g., the p-values between two clusters from the null data.

  - nlogTrans:
    
    A logical value. If the input scores are p-values, take the `-log10`
    transformation since Clipper require larger scores represent more
    significant DE. Default is TRUE.

  - FDR:
    
    A numeric value of the target False Discovery Rate (FDR). Must be
    'diff' or 'max'.

  - contrastScore:
    
    A string value of the way to construct contrast scores. The choice
    can be

  - correct:
    
    A logical value. If TRUE, perform the correction to make the
    distribution of contrast scores approximately symmetric. Default is
    FALSE.

  - threshold:
    
    A string value of the threshold method. Must be 'BC' or 'DS'.

  - ordering:
    
    A logic value. If TRUE, order the genes in the returned table by
    their significance. Default is TRUE.

## Value

A list of target FDR, DE genes, and the detailed summary table.

## Details

This function constructs the contrast scores by taking the difference
between target DE scores and null DE scores.

## Examples

``` r
targetScores <- runif(10000)
nullScores <- runif(10000)
names(targetScores) <- names(nullScores) <- paste0("Gene", 1:10000)
res <- callDE(targetScores, nullScores, correct = FALSE)
```
