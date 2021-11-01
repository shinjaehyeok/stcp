# Gaussian case
# Pre-change : N(0,1)
# Post-change: N(mu, sig)
# Change-point: v = 0 (immediate chane), 200, 500, 1000 (no change)

max_sample <- 1000L
mu <- 1
sig <- sqrt(2)
v <- 200
alpha <- 1/max_sample

# Compute optimal delta star
delta_star <- mu / sig^2
# delta_star <- 2
delta_upper <- 10
delta_lower <- 0.001


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

# e-detector 1. SR-type
current_ind <- 1
prev_log_e <- 0
e_detect_val_SR <- numeric(max_sample)
updater <- function(){
  prev_log_e <<- update_log_e_detector(
    x_current = x_vec[current_ind],
    prev_log_e = prev_log_e,
    compute_log_baseline = compute_log_baseline,
    is_SR_type = TRUE
  )
  e_detect_val_SR[current_ind] <<- prev_log_e
  current_ind <<- current_ind + 1
}
while(current_ind <= max_sample){
  updater()
}

# e-dector 2. CUSUM-type
current_ind <- 1
prev_log_e <- 0
e_detect_val_CS <- numeric(max_sample)
updater <- function(){
  prev_log_e <<- update_log_e_detector(
    x_current = x_vec[current_ind],
    prev_log_e = prev_log_e,
    compute_log_baseline = compute_log_baseline,
    is_SR_type = FALSE
  )
  e_detect_val_CS[current_ind] <<- prev_log_e
  current_ind <<- current_ind + 1
}
while(current_ind <= max_sample){
  updater()
}

plot(1:max_sample, e_val_vec, type = "l",
     xlab = "n", ylab = "log e-val", main = paste0("v = ", v))
abline(h = base_param$g_alpha, col = 2)
abline(v = v, col = 1, lty = 2)
lines(1:max_sample, e_detect_val_SR, type = "l", col = 3)
lines(1:max_sample, e_detect_val_CS, type = "l", col = 4)

g_alpha_single <- base_param$g_alpha

# When delta_lower < delta_star < delta_upper
base_param <- compute_baseline(
  alpha,
  delta_lower = delta_lower,
  delta_upper = delta_upper
)

# Compute e-detectors
generate_log_base_fn <- function(lambda){
  compute_log_baseline <- function(x){
    lambda * x - lambda^2 / 2
  }
  return(compute_log_baseline)
}
compute_log_baseline_fn_list <- lapply(base_param$lambda, generate_log_base_fn)


# e-value for testing (it is not the scope of this package but good for debugging)
e_val_mat <- sapply(compute_log_baseline_fn_list, function(f) f(x_vec))
log_omega <- log(base_param$omega)
mix_e_val_vec <- apply(e_val_mat, 1, function(e_vec) matrixStats::logSumExp(e_vec + log_omega)) %>% cumsum()
# plot(1:max_sample, mix_e_val_vec, type = "l")
# abline(h = base_param$g_alpha, col = 2)
# lines(1:max_sample, e_val_vec, col = 3)

# e-detector 1. SR-type
current_ind <- 1
prev_log_e_vec <- numeric(length(base_param$lambda))
mix_e_detect_val_SR <- numeric(max_sample)

updater <- function(is_SR_type = TRUE){
  update_log_e_for_mix_ind <- function(mix_ind){
    update_log_e_detector(
      x_current = x_vec[current_ind],
      prev_log_e = prev_log_e_vec[mix_ind],
      compute_log_baseline = compute_log_baseline_fn_list[[mix_ind]],
      is_SR_type = is_SR_type
    )
  }
  prev_log_e_vec <<- sapply(seq_along(base_param$lambda), update_log_e_for_mix_ind)
  mix_e_detect_val_SR[current_ind] <<- matrixStats::logSumExp(prev_log_e_vec + log_omega)
  current_ind <<- current_ind + 1
}
while(current_ind <= max_sample){
  updater(is_SR_type = TRUE)
}

# e-dector 2. CUSUM-type
current_ind <- 1
prev_log_e_vec <- numeric(length(base_param$lambda))
mix_e_detect_val_CS <- numeric(max_sample)

updater <- function(is_SR_type = TRUE){
  update_log_e_for_mix_ind <- function(mix_ind){
    update_log_e_detector(
      x_current = x_vec[current_ind],
      prev_log_e = prev_log_e_vec[mix_ind],
      compute_log_baseline = compute_log_baseline_fn_list[[mix_ind]],
      is_SR_type = is_SR_type
    )
  }
  prev_log_e_vec <<- sapply(seq_along(base_param$lambda), update_log_e_for_mix_ind)
  mix_e_detect_val_CS[current_ind] <<- matrixStats::logSumExp(prev_log_e_vec + log_omega)
  current_ind <<- current_ind + 1
}
while(current_ind <= max_sample){
  updater(is_SR_type = FALSE)
}

SR_stop <- which(e_detect_val_SR > g_alpha_single) %>% min()
CS_stop <- which(e_detect_val_CS > g_alpha_single) %>% min()
mix_SR_stop <- which(mix_e_detect_val_SR > base_param$g_alpha) %>% min()
mix_CS_stop <- which(mix_e_detect_val_SR > base_param$g_alpha) %>% min()

plot(1:max_sample, mix_e_val_vec, type = "l",
     xlab = "n", ylab = "log e-val", main = paste0("v = ", v),
     xlim = c(0, max_sample),
     ylim = c(min(e_val_vec), max(mix_e_detect_val_SR)))
# plot(1:max_sample, mix_e_val_vec, type = "l",
#      xlab = "n", ylab = "log e-val", main = paste0("v = ", v),
#      xlim = c(v, v+20))
abline(h = base_param$g_alpha, col = 2)
abline(h = g_alpha_single, col = 2, lty = 2)
abline(v = v, col = 1, lty = 2)
lines(1:max_sample, e_val_vec, type = "l", lty = 2)
lines(1:max_sample, mix_e_detect_val_SR, type = "l", col = 3)
abline(v = mix_SR_stop, col = 3)
lines(1:max_sample, e_detect_val_SR, type = "l", col = 3, lty =2)
abline(v = SR_stop, col = 3, lty = 2)
lines(1:max_sample, mix_e_detect_val_CS, type = "l", col = 4)
abline(v = mix_CS_stop, col = 4)
lines(1:max_sample, e_detect_val_CS, type = "l", col = 4, lty = 2)
abline(v = CS_stop, col = 4, lty = 2)

print(SR_stop)
print(CS_stop)
print(mix_SR_stop)
print(mix_CS_stop)

