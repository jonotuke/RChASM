# summary_calls

Combine both autosomal and sex chromosomal analyses into one output

## Usage

``` r
summary_calls(
  inChASM,
  minTotal = 60000,
  minPosterior = 0.95,
  ignoreUnusual = FALSE,
  printProtocol = FALSE
)
```

## Arguments

- inChASM:

  the result of combining both autosomal and sca analyses

- minTotal:

  minimum total for autosomal or sca analyses

- minPosterior:

  minimum value maxP can take (i.e. minimum posterior probability
  accepted)

- ignoreUnusual:

  filter our unusual observations?

- printProtocol:

  print protocol column?

## Value

combined calls
