# If EDCP is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/EDCP")

library(EDCP)
set.seed(1)
# Bernoulli case ----
# Ber(p)
max_sample <- 1000L
pre_sample <- 10L
p <- 0.3
alpha <- 0.025
v_min <- 1
k_max <- 1e+3

check_bf_test <- FALSE # Use it only for debugging purpose.

# Compute target interval
n_target <- round(max_sample / 2)
n_lower <- max_sample / 10
n_upper <- max_sample * 2

# Generate data
x_vec <-  rbinom(max_sample, 1, p)
x_bar <- cumsum(x_vec) / seq_along(x_vec)

# Build CI
# When delta_lower = delta_upper = delta_star
# We can use sub_B fn of different p. Here, we use 0.5. Ideally, it should be similar to the true underlying p.
# This make CI computation efficient as we use fixed parameters but it makes potentially wider CI compared to the brute-force.
baseline_ci_star <- compute_baseline_for_sample_size(alpha,
                                                     n_target,
                                                     n_target,
                                                     generate_sub_B_fn(0.5),
                                                     v_min,
                                                     k_max)

ci_helper_star <- generate_ci_helper(
  baseline_ci_star$alpha,
  baseline_ci_star$omega,
  baseline_ci_star$lambda,
  generate_psi_fn_list = generate_sub_B_fn,
  is_psi_depend_on_m = TRUE
)
ci_star <- compute_ci(x_vec, ci_helper_star, ci_lower_trivial = 0)

# When delta_lower < delta_star < delta_upper
baseline_ci_mix <- compute_baseline_for_sample_size(alpha,
                                                    n_upper,
                                                    n_lower,
                                                    generate_sub_B_fn(0.5),
                                                    v_min,
                                                    k_max)

ci_helper_mix <- generate_ci_helper(
  baseline_ci_mix$alpha,
  baseline_ci_mix$omega,
  baseline_ci_mix$lambda,
  generate_psi_fn_list = generate_sub_B_fn,
  is_psi_depend_on_m = TRUE
)
ci_mix <- compute_ci(x_vec, ci_helper_mix, ci_lower_trivial = 0)

# If we want to treat this example as an additive sub-psi case
# then set is_psi_depend_on_m = FALSE and ci_lower_trivial = -Inf
# It will return  more conservative CI but the width now does not depend on x_bar (so you can pre-compute).
# The additional conservatives depends on the gap between true p and the parameter used in CI construction.
# This is because by treating it as an additive sub-psi, we no longer adaptively track the variance term.
# Hence, use it only if you have a good prior knowledge on it.
ci_helper_mix2 <- generate_ci_helper(
  baseline_ci_mix$alpha,
  baseline_ci_mix$omega,
  baseline_ci_mix$lambda,
  generate_psi_fn_list = generate_sub_B_fn,
  is_psi_depend_on_m = FALSE
)
ci_mix2 <- compute_ci(x_vec, ci_helper_mix2)


ci_fixed <-
  x_bar - sqrt(p * (1 - p)) * qnorm(alpha, lower.tail = FALSE) / sqrt(1:max_sample)

# Plot CI
plot(
  1:max_sample,
  x_bar,
  type = "l",
  xlab = "n",
  ylab = "X_bar",
  main = "Running average",
  ylim = c(0, 1),
)
graphics::abline(h = p,
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
graphics::lines(1:max_sample,
                ci_mix2$ci_lower,
                lty = 2,
                col = 4)


# Compare width
plot(
  1:max_sample,
  (x_bar - ci_star$ci_lower) / (x_bar - ci_fixed),
  type = "l",
  ylim = c(1, 2),
  col = 2
)
graphics::lines(1:max_sample, (x_bar - ci_mix$ci_lower) / (x_bar - ci_fixed), col = 3)
graphics::lines(1:max_sample, (x_bar - ci_mix2$ci_lower) / (x_bar - ci_fixed), col = 4)

# Check correctness via brute-force method
# For the Bernoulli case,
if (check_bf_test) {
  brute_force_ci_helper <- function(m, x_vec, tol = 1e-6) {
    if (m <= 0)
      m <- tol
    if (m + baseline_ci_mix$delta_upper >= 1)
      m <- 1 - baseline_ci_mix$delta_upper - tol

    base_param <- compute_baseline_for_sample_size(alpha,
                                                   n_upper,
                                                   n_lower,
                                                   generate_sub_B_fn(m),
                                                   v_min,
                                                   k_max)

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
      compute_brute_force_ci(
        x_vec[1:n],
        brute_force_ci_helper,
        width_upper = 100,
        ci_lower_trivial = 0
      )
    })

  ratio_vec <-
    (x_bar[n_bf_vec] - ci_mix$ci_lower[n_bf_vec]) / (x_bar[n_bf_vec] - ci_bf)
  print(ratio_vec)
  # Note in general bf method yields tighter bound by trading off increased computations.
}
