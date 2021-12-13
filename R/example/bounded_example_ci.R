# If EDCP is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/EDCP")

library(EDCP)
set.seed(1)
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
n_lower <- max_sample / 5
n_upper <- max_sample * 2

# Generate data
x_vec <-  5 * rbeta(max_sample, a, b) - 1
x_bar <- cumsum(x_vec) / seq_along(x_vec)

# Build CI
ci_fixed <-
  sig * qnorm(alpha, lower.tail = FALSE) / sqrt(1:max_sample)

base_param <- compute_baseline_for_sample_size(alpha,
                                               n_upper,
                                               n_lower,
                                               psi_fn_list,
                                               v_min,
                                               k_max)


brute_force_ci_helper <- function(m, x_vec, tol = 1e-6) {
  log_base_fn_list <- sapply(
    base_param$lambda,
    generate_log_base_fn,
    psi_fn = base_param$psi_fn_list$psi,
    s_fn = function(x) {
      x - m
    }
  )

  baseline_m <- list(
    omega = base_param$omega,
    log_base_fn_list =  log_base_fn_list,
    alpha = alpha,
    lambda = base_param$lambda,
    g_alpha = base_param$g_alpha
  )

  e_val <- edcp(x_vec,
                baseline_m,
                is_test = TRUE)
  return(e_val$log_mix_e_val[length(e_val$log_mix_e_val)] - e_val$threshold)
}
# Warning::This code is very slow O(n^2)
n_bf_vec <- c(seq(1, length(x_vec), by = 10L), length(x_vec))
ci_bf <-
  sapply(n_bf_vec, function(n) {
    compute_brute_force_ci(x_vec[1:n],
                           brute_force_ci_helper,
                           width_upper = 100)
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
graphics::lines(n_bf_vec, ci_bf, lty = 2, col = 2)


# Compare width
plot(n_bf_vec,
     (x_bar[n_bf_vec] - ci_bf) / ci_fixed[n_bf_vec],
     type = "l",
     ylim = c(1, 2))
# graphics::lines(n_bf_vec, (x_bar[n_bf_vec] - ci_bf_mix) / ci_fixed[n_bf_vec], lty = 2, col = 2)
