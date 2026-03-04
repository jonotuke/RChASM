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

## Examples

``` r
example_calls <- runChASM(rawReadCountsIn = example_data)

summary_calls(inChASM = example_calls, minTotal = 6e4, minPosterior = 0.95)
#> # A tibble: 3 × 10
#>   sample  unusual flags autosomal_call sca_call C_call autosomal_total sca_total
#>   <chr>   <lgl>   <int> <chr>          <chr>    <chr>            <dbl>     <dbl>
#> 1 Ind_71… TRUE        4 No Aneuploidy  XXX      No Si…           70688     76444
#> 2 Ind_18… FALSE       1 No Aneuploidy  XXX      No Si…          209933    225713
#> 3 Ind_25… FALSE       0 Trisomy 21     XX       No Si…          254544    268132
#> # ℹ 2 more variables: automsomal_maxP <dbl>, sca_maxP <dbl>
```
