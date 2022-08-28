
<!-- README.md is generated from README.Rmd. Please edit that file -->

# stcp: Sequential Test and Change-Point detection algorithms based on E-values / E-detectors

<!-- badges: start -->

[![R-CMD-check](https://github.com/shinjaehyeok/stcp/workflows/R-CMD-check/badge.svg)](https://github.com/shinjaehyeok/stcp/actions)
<!-- badges: end -->

stcp is an R package built to run nonparametric sequential tests and
online change point detection algorithms in \[CITE SRR 21’ and SRR
22’\]. This package supports APIs of nonparametric sequential test and
online change-point detection for streams of univariate sub-Gaussian,
binary, and bounded random variables. This package also provides APIs of
general E-value based sequential test and E-detectors based change-point
detection frameworks for advanced users.

## Installation

You can install the development version of stcp from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("shinjaehyeok/stcp")
```

<!-- You'll still need to render `README.Rmd` regularly, to keep `README.md` up-to-date. `devtools::build_readme()` is handy for this.  -->
<!-- TODO 1. Make test also work for list; 2. To make CI works for variance process-->
