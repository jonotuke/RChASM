# combine data

Combine both autosomal and sex chromosomal analyses into one output

## Usage

``` r
combineData(calls.auto, calls.sca, z.scores, printMissingIDs = FALSE)
```

## Arguments

- calls.auto:

  a list or predefined piped string of sample IDs to look for

- calls.sca:

  the autosomal karyotype calls for the samples

- z.scores:

  the sex chromosomal karyotype calls for the samples

- printMissingIDs:

  print the IDs of individuals not in all 3 input files (or just
  numbers)?

## Value

combined data
