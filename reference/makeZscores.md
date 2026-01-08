# Calculates Z-scores per chromosome

Calculates Z-scores per chromosome

## Usage

``` r
makeZscores(indat, refType = "auto", min_reads = 30000, max_reads = 1e+09)
```

## Arguments

- indat:

  full input tibble

- refType:

  "auto"="autosomal"; "sca"="sex chromosomal"; "diagnostic"

- min_reads:

  min number of reads per sample

- max_reads:

  max number of reads per sample

## Value

Z-scores
