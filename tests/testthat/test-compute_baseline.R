test_that("Gaussian case works", {
# Gaussian case
# Pre-change : N(0,1)
# Post-change: N(mu, sig)
# Change-point: v = 0 (immediate chane), 200, 500, 1000 (no change)

max_sample <- 1000L
mu <- 1
sig <- sqrt(2)
v <- max_sample
alpha <- 1/max_sample

# Compute optimal delta star
delta_star <- mu / sig^2

# Generate sample
  generator <- function(v = 100){
    x_vec <- numeric(max_sample)
    if (v > 0){
      x_vec[1:v] <- rnorm(v)
    }
    if (v < max_sample){
      x_vec[seq(v+1, max_sample)] <- rnorm(max_sample - v, mu, sig)
    }
    return(x_vec)
  }
  x_vec <- generator(v)
  # plot(1:max_sample, x_vec)

  # When delta_lower = delta_upper = delta_star
  base_param <- compute_baseline(
    alpha,
    delta_lower = delta_star,
    delta_upper = delta_star
    )

  # It should return a single lambda with trivial weight and threshold
  expect_equal(length(base_param$lambda), 1)
  expect_equal(base_param$omega, 1)
  expect_true(abs(base_param$g_alpha - log(1/alpha)) < 1e-8)

  # Compute e-detectors
  compute_log_baseline <- function(x){
    base_param$lambda * x - base_param$lambda^2 / 2
  }

  # e-value for testing (it is not the scope of this package but good for debugging)
  e_val_vec <- cumsum(compute_log_baseline(x_vec))
  # plot(1:max_sample, e_val_vec, type = "l")
  # abline(h = base_param$g_alpha, col = 2)

  # e-detector 1. SR type
  current_ind <- 1
  prev_log_e <- 0
  e_detect_val <- numeric(max_sample)
  updater <- function(){
    prev_log_e <<- update_log_e_detector(
      x_current = x_vec[current_ind],
      prev_log_e = prev_log_e,
      compute_log_baseline = compute_log_baseline,
      is_SR_type = FALSE
    )
    e_detect_val[current_ind] <<- prev_log_e
    current_ind <<- current_ind + 1
  }
  while(current_ind <= max_sample){
    updater()
  }
  # plot(1:max_sample, e_val_vec, type = "l")
  # abline(h = base_param$g_alpha, col = 2)
  # lines(1:max_sample, e_detect_val_SR, type = "l", col = 3)
  # lines(1:max_sample, e_detect_val, type = "l", col = 4)
  # plot(1:max_sample, e_detect_val, type = "l")
  # abline(h = base_param$g_alpha, col = 2)

})

