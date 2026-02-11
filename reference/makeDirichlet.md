# function to make a Dirichlet prior

function to make a Dirichlet prior

## Usage

``` r
makeDirichlet(
  indat,
  refType,
  min_reads = 30000,
  max_reads = 1e+09,
  show_plot = TRUE
)
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

- show_plot:

  boolean to show plot

## Value

Dirichlet prior
