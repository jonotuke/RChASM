# calculate Dirichlet probabilities with filtering

calculate Dirichlet probabilities with filtering

## Usage

``` r
ddirichletultinomial(nvector, avector, cvector, log = TRUE, correction = 1e-06)
```

## Arguments

- nvector:

  observed read counts

- avector:

  alpha vector to change

- cvector:

  change vector for new karyotype

- log:

  return value in log space

- correction:

  error rate for Y mapping for non-Y karyotypes

## Value

Dirichlet probabilities
