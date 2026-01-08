# Plots diagnostic plots

Plots diagnostic plots

## Usage

``` r
plot_diagnostic(
  IDs,
  calls.auto,
  calls.sca,
  z.scores,
  inDirichlet.auto,
  inDirichlet.sca,
  min_reads = 30000,
  max_reads = 1e+09,
  addLabels = FALSE
)
```

## Arguments

- IDs:

  a list or predefined piped string of sample IDs to look for

- calls.auto:

  the autosomal karyotype calls for the samples

- calls.sca:

  the sex chromosomal karyotype calls for the samples

- z.scores:

  the autosomal z-scores for the samples

- inDirichlet.auto:

  the dirichlet parameters for the autosomes

- inDirichlet.sca:

  the dirichlet parameters for the sex chromosomes

- min_reads:

  min number of reads per sample

- max_reads:

  max number of reads per sample

- addLabels:

  add "A", "B" and "C" to the plot panels?

## Value

plots
