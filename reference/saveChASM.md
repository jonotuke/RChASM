# A function for saving the results of a ChASM analysis as a tsv

A function for saving the results of a ChASM analysis as a tsv

## Usage

``` r
saveChASM(inChASM, file, sort_by_samplename = FALSE)
```

## Arguments

- inChASM:

  the result of a full ChASM analysis (output from runChASM)

- file:

  the full path and file name to place output

- sort_by_samplename:

  reorder alphabetically by sample names?

## Value

saves results

## Examples

``` r
example_calls <- runChASM(rawReadCountsIn = example_data)

if (FALSE) { # \dontrun{
saveChASM(inChASM = example_calls)
} # }
```
