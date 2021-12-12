# If EDCP is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/EDCP")

library(EDCP)

# Bernoulli case ----
# 5 * Beta(2, 3) - 1 (m = 1, sig^2 = 5^2 * 6 / (5^2 * 6) = 1)

max_sample <- 1000L
a <- 2
b <- 3
m <- 1
sig <- 1
bound_lower <- -1
bound_upper <- 5
alpha <- 0.025

psi_fn_list <- generate_sub_E_fn()
v_min <- 0
k_max <- 1e+3

check_bf_test <- TRUE

# Compute target interval
n_lower <- max_sample / 10
n_upper <- max_sample * 2
delta_mix <- convert_time_to_delta_bound(alpha,
                                         n_upper,
                                         n_lower,
                                         psi_fn_list,
                                         v_min,
                                         k_max)

# Generate data
x_vec <-  5 * rbeta(max_sample, a, b) - 1
x_bar <- cumsum(x_vec) / seq_along(x_vec)

# Build CI
ci_fixed <-
  sig * qnorm(alpha, lower.tail = FALSE) / sqrt(1:max_sample)

# For bounded case, ci_helper based method does not yield a valid CI but it can be used as a proxy.
baseline_mix <- compute_baseline_bounded(
  alpha,
  (bound_lower + bound_upper) / 2, # Dummy value
  delta_lower = 1, # Dummy value
  delta_upper = 1, # Dummy value
  bound_lower = bound_lower,
  bound_upper = bound_upper,
  k_max = k_max,
  delta_lower_sub_E = delta_mix$delta_lower,
  delta_upper_sub_E = delta_mix$delta_upper
)

# For bounded case, it return ci ratio
ci_mix <- compute_ci_width(x_vec,
                           baseline_mix$ci_helper,
                           ci_upper_const = 1000)
ci_mix_lower <- ci_mix$x_bar / ci_mix$ci_width

# For the bounded case, we can only use BF method to get a valid CI.
bruth_force_ci_helper <- function(m, x_vec) {
  baseline_m <- compute_baseline_bounded(
    alpha,
    max(m, bound_lower + 1e-6),
    delta_lower = 1, # Dummy value
    delta_upper = 1, # Dummy value
    bound_lower,
    bound_upper,
    k_max,
    delta_lower_sub_E = delta_mix$delta_lower,
    delta_upper_sub_E = delta_mix$delta_upper
  )
  e_val <- edcp(x_vec,
                baseline_m,
                is_test = TRUE)
  return(last(e_val$log_mix_e_val) - e_val$threshold)
}

# Warning::This code is very slow O(n^2)
n_bf_vec <- c(seq(1, length(x_vec), by = 10L), length(x_vec))
ci_bf <-
  sapply(n_bf_vec, function(n) {
    brute_force_lower_ci(x_vec[1:n],
                         bruth_force_ci_helper,
                         ci_upper_const = 100)
  })

# Plot CI

plot(
  1:max_sample,
  x_bar,
  type = "l",
  xlab = "n",
  ylab = "X_bar",
  main = "Running average",
  ylim = c(-0.5, 0.5) + m,
)
graphics::abline(h = m,
                 col = 2,
                 lwd = 2)

graphics::lines(1:max_sample,
                x_bar - ci_fixed,
                lty = 2,
                col = 1)
graphics::lines(1:max_sample, ci_mix_lower, lty = 2, col = 2)
graphics::lines(n_bf_vec, ci_bf, lty = 2, col = 3)


# Compare width
plot(n_bf_vec,
     (x_bar[n_bf_vec] - ci_bf) / ci_fixed[n_bf_vec],
     type = "l")
