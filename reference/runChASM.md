# Title

Title

## Usage

``` r
runChASM(
  rawReadCountsIn,
  minSamplesPerProtocol = 30,
  min_reads = 60000,
  max_reads = 1e+09,
  p_contamination = 0.01,
  show_plot = TRUE,
  printMissingIDs = FALSE
)
```

## Arguments

- rawReadCountsIn:

  the reads counts for each chromosome

- minSamplesPerProtocol:

  minimum number of reads per protocol for parameter estimation

- min_reads:

  the minimum number of reads for Dirichlet parameter estimation

- max_reads:

  the maximum number of reads for Dirichlet parameter estimation

- p_contamination:

  the probability of a sample yielding significant contamination

- show_plot:

  show the clustering plot for sex chromosomal aneuploidy Dirichlet
  estimation?

- printMissingIDs:

  when combining karyotype calls, return names that are missing (or just
  the number of missing IDs)?

## Value

summary
