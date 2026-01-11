
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RChASM

<!-- badges: start -->

[![R-CMD-check](https://github.com/jonotuke/RChASM/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/jonotuke/RChASM/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

RChASM: A Statistically Rigorous Method for the Detection of Chromosomal
Aneuploidies in Ancient DNA Studies.

RChASM is an R implementation of ChASM (Chromosomal Aneuploidy Screening
Methodology): a statistically rigorous Bayesian approach for screening
data sets for autosomal and sex chromosomal aneuploidies. RChASM takes
as input the number of (deduplicated) reads mapping to chromosomes 1-22
and the X and Y chromosomes, and models these using a
Dirichlet-multinomial distribution. From this, RChASM returns posterior
probabilities of sex chromosomal karyotypes (XX, XY, XXY, XYY, XXX and
X) and full autosomal aneuploidies (trisomy 13, trisomy 18 and trisomy
21). RChASM also returns two diagnostic statistics: (i) a posterior
probability addressing whether contamination between XX and XY may
explain the observed sex chromosomal aneuploidy, and (ii) a chi-squared
statistic measuring whether the observed read counts are too divergent
from the underlying distribution (and may represent abnormal
sequencing/quality issues).

## Installation

You can install the development version of RChASM from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("jonotuke/RChASM")
```

## Example

Please see `vignette("example_analysis")` for a tutorial
