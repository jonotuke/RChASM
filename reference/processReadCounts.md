# processReadCounts

Task: takes read counts and makes useful file

## Usage

``` r
processReadCounts(
  rawReadCountsIn,
  refType,
  minTotal = 1000,
  minSamplesPerProtocol = 30
)
```

## Arguments

- rawReadCountsIn:

  the reads counts via KP

- refType:

  for capture, use only on-target (F)?

- minTotal:

  minimum number of protocol counts to be included

- minSamplesPerProtocol:

  "auto"="autosomal"; "sca"="sex chromosomal"

## Value

outreadcounts
