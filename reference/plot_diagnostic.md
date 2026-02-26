# Plots diagnostic plots

Plots diagnostic plots

## Usage

``` r
plot_diagnostic(IDs, inChASM, addLabels = FALSE)
```

## Arguments

- IDs:

  a list or predefined piped string of sample IDs to look for

- inChASM:

  parameter object

- addLabels:

  add "A", "B" and "C" to the plot panels?

## Value

plots

## Examples

``` r
example_calls <- runChASM(rawReadCountsIn = example_data)

plot_diagnostic(IDs = 'Ind_255_1', inChASM = example_calls, addLabels = TRUE)
#> Warning: Removed 1 row containing missing values or values outside the scale range
#> (`geom_text()`).

```
