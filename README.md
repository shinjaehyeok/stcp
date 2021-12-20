
<!-- README.md is generated from README.Rmd. Please edit that file -->

# EDCP: E-Detector-based online Change-Point detection algorithms

<!-- badges: start -->

[![R-CMD-check](https://github.com/shinjaehyeok/EDCP/workflows/R-CMD-check/badge.svg)](https://github.com/shinjaehyeok/EDCP/actions)
<!-- badges: end -->

EDCP is a R package built to run e-detector-based nonparametric online
change-point detection algorithms in \[CITE PAPER\].

## Installation

You can install the development version of EDCP from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("shinjaehyeok/EDCP")
```

## Example

Suppose we have a stream of observations
*X*<sub>1</sub>, *X*<sub>2</sub>, … ∈ \[0,1\]. Before an unknown
change-point *ν*, each pre-change observation has mean less than or
equal to 0.3. But after the change-point *ν*, each post-change
observation has mean larger than 0.3.

``` r
set.seed(1)
library(EDCP)
# Generate a stream of observations
# For simplicity, we use beta distributions to generate samples
max_sample <- 1000L
nu <- 500L
x_vec <- numeric(max_sample)
x_vec[1:nu] <- rbeta(nu, 1.5, 3.5) # Pre-change dist: B(1.5, 3.5)
x_vec[seq(nu+1, max_sample)] <- rbeta((max_sample - nu), 2, 3) # Post-change dist: B(2,3)
m_pre <- 0.3
m_post <- 0.4
plot(1:max_sample, x_vec, pch=20, 
     xlab = "n", ylab = "X_n", main = "Simulated Data")
```

<img src="man/figures/README-example-1.png" width="100%" />

Below, we compute a mixture e-detectors by which we can infer when the
change of the mean happened.

``` r
alpha <- 1e-3 # Inverse of target ARL
m_pre <- 0.3 # Upper bound of mean of pre-change observations
delta_lower <- 0.01  # Guess on the minimum gap between pre- and post-change means
baseline_bounded <- compute_baseline_bounded(alpha,
                                             m_pre,
                                             delta_lower)

# Compute mixture of SR-type e-detectors.
mix_SR_bounded <- edcp(x_vec, baseline_bounded)
mix_SR_stop <- mix_SR_bounded$stopped_ind
threshold <- mix_SR_bounded$threshold # = log(1/alpha)

# Plot 
plot(1:max_sample, mix_SR_bounded$log_mix_e_val, type = "l",
     xlab = "n", ylab = "log e-detectors", main = paste0("v = ", nu, " v_hat = ",mix_SR_stop))
abline(h = threshold, col = 2)
abline(v = nu, col = 1, lty = 2)
abline(v = mix_SR_stop, col = 2, lty = 2)
```

<img src="man/figures/README-edcp-1.png" width="100%" />

``` r
plot(1:max_sample, x_vec, pch=20, 
     xlab = "n", ylab = "X_n", main = "Simulated Data with the detected CP")
lines(x = c(0,nu), y = c(0.3,0.3), col = 2, lwd = 2)
lines(x = c(nu,max_sample), y = c(0.4,0.4), col = 2, lwd =2 )
abline(v = nu, col = 1, lty = 2)
abline(v = mix_SR_stop, col = 2, lty = 2, lwd = 2)
```

<img src="man/figures/README-example2-1.png" width="100%" /> You’ll
still need to render `README.Rmd` regularly, to keep `README.md`
up-to-date. `devtools::build_readme()` is handy for this.
