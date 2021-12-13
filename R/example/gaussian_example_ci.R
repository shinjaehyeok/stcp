# If EDCP is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/EDCP")

library(EDCP)

# Gaussian case ----
# N(mu, 1)
max_sample <- 1000L
mu <- -1
sig <- 1
alpha <- 0.025

psi_fn_list <- generate_sub_G_fn()
v_min <- 1
k_max <- 1e+3

check_bf_test <- FALSE

# Generate data
x_vec <- rnorm(max_sample, mu)
x_bar <- cumsum(x_vec) / seq_along(x_vec)

# Build CI
# When delta_lower = delta_upper = delta_star
# Compute optimal delta star
n_target <- round(max_sample / 2)
baseline_ci_star <- compute_baseline_for_sample_size(alpha,
                                                     n_target,
                                                     n_target,
                                                     psi_fn_list,
                                                     v_min,
                                                     k_max)
ci_helper_star <- generate_ci_helper(
  baseline_ci_star$alpha,
  baseline_ci_star$omega,
  baseline_ci_star$lambda,
  generate_psi_fn_list = generate_sub_G_fn,
  is_psi_depend_on_m = FALSE
)
ci_star <- compute_ci(x_vec, ci_helper_star)

# When delta_lower < delta_star < delta_upper
# Compute target interval
n_lower <- max_sample / 5
n_upper <- max_sample * 2
baseline_ci_mix <- compute_baseline_for_sample_size(alpha,
                                                    n_upper,
                                                    n_lower,
                                                    psi_fn_list,
                                                    v_min,
                                                    k_max)

ci_helper_mix <- generate_ci_helper(
  baseline_ci_mix$alpha,
  baseline_ci_mix$omega,
  baseline_ci_mix$lambda,
  generate_psi_fn_list = generate_sub_G_fn,
  is_psi_depend_on_m = FALSE
)
ci_mix <- compute_ci(x_vec, ci_helper_mix)

ci_fixed <-
  x_bar - qnorm(alpha, lower.tail = FALSE) / sqrt(1:max_sample)

# Plot CI

plot(
  1:max_sample,
  x_bar,
  type = "l",
  xlab = "n",
  ylab = "X_bar",
  main = "Running average",
  ylim = c(-0.5, 0.5) + mu,
)
graphics::abline(h = mu,
                 col = 2,
                 lwd = 2)

graphics::lines(1:max_sample,
                ci_fixed,
                lty = 2,
                col = 1)
graphics::lines(1:max_sample,
                ci_star$ci_lower,
                lty = 2,
                col = 2)
graphics::lines(1:max_sample,
                ci_mix$ci_lower,
                lty = 2,
                col = 3)

# Compare width
plot(
  1:max_sample,
  (x_bar - ci_star$ci_lower) / (x_bar - ci_fixed),
  type = "l",
  ylim = c(1, 2)
)
graphics::lines(1:max_sample, (x_bar - ci_mix$ci_lower) / (x_bar - ci_fixed), col = 2)

# Check correctness via brute-force method
if (check_bf_test) {
  bruth_force_ci_helper <- function(m, x_vec) {
    baseline_m <- compute_baseline_simple_exp(
      alpha,
      m,
      baseline_ci_mix$delta_lower,
      baseline_ci_mix$delta_upper,
      baseline_ci_mix$psi_fn_list,
      s_fn = function(x) {
        x - m
      },
      v_min = v_min,
      k_max = k_max
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
                             bruth_force_ci_helper,
                             width_upper = 100)
    })

  ratio_vec <- ci_mix$ci_lower[n_bf_vec] / ci_bf
  if (abs(mean(ratio_vec) - 1) > 1e-8)
    stop("BF test failed")
  if (var(ratio_vec) > 1e-8)
    stop("BF test failed")
}
