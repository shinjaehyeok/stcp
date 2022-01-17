# If STCP is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/STCP")

library(STCP)
set.seed(1)
# Bernoulli case ----
# Ber(p)
max_sample <- 1000L
pre_sample <- 10L
p <- 0.3
p_guess <- 0.5
alpha <- 0.025
v_min <- 1
k_max <- 1e+3
# Make a wrapper to change default value.
psi_fn_list_generator <- function(p = p_guess) {
  generate_sub_B_fn(p)
}

check_bf_test <- TRUE # Use it only for debugging purpose.

# Compute target interval
n_target <- round(max_sample / 2)
n_lower <- max_sample / 5
n_upper <- max_sample * 2

# Generate data
x_vec <-  rbinom(max_sample, 1, p)
x_bar <- cumsum(x_vec) / seq_along(x_vec)

# Build CI
# When delta_lower = delta_upper = delta_star
# We can use sub_B fn of different p. Here, we use p_guess. Ideally, it should be similar to the true underlying p.
# This make CI computation efficient as we use fixed parameters but it makes potentially wider CI compared to the brute-force.
ci_model_star <- build_ci_exp(alpha,
                              n_target,
                              n_target,
                              psi_fn_list_generator,
                              v_min = v_min,
                              k_max = k_max)

ci_star <- compute_ci(x_vec, ci_model_star, ci_lower_trivial = 0)



# When delta_lower < delta_star < delta_upper
ci_model_mix <- build_ci_exp(alpha,
                             n_upper,
                             n_lower,
                             psi_fn_list_generator,
                             v_min = v_min,
                             k_max = k_max)

ci_mix <- compute_ci(x_vec, ci_model_mix, ci_lower_trivial = 0)



# If we want to treat this example as an additive sub-psi case
# then set is_psi_depend_on_m = FALSE and ci_lower_trivial = -Inf
# It will return  more conservative CI but the width now does not depend on x_bar (so you can pre-compute it).
# The degree of conservative depends on the gap between true p and the parameter used in CI construction.
# This is because by treating it as an additive sub-psi, we no longer adaptively track the variance term.
# Hence, use it only if you have a good prior knowledge on it.

# This procedure is not officially supported by this package.
ci_helper_mix2 <- generate_ci_helper(
  ci_model_mix$baseline_obj$alpha,
  ci_model_mix$baseline_obj$omega,
  ci_model_mix$baseline_obj$lambda,
  psi_fn_list_generator = psi_fn_list_generator,
  is_psi_depend_on_m = FALSE
)
tmp <- list(ci_helper = ci_helper_mix2)
class(tmp) <- "STCP_CI"
ci_mix2 <- compute_ci(x_vec, tmp)



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
graphics::lines(1:max_sample,
                (x_bar - ci_mix2$ci_lower) / (x_bar - ci_fixed),
                col = 4)

# Check correctness via brute-force method
# For the Bernoulli case,
if (check_bf_test) {
  bf_ci_model <- build_bf_ci_exp(alpha,
                                 n_upper,
                                 n_lower,
                                 psi_fn_list_generator,
                                 v_min = v_min,
                                 k_max = k_max)

  bf_ci <- compute_bf_ci(x_vec,
                         bf_ci_model,
                         ci_lower_trivial = 0,
                         max_num_ci = 100)

  ratio_vec <-
    (x_bar[bf_ci$n] - ci_mix$ci_lower[bf_ci$n]) / ( bf_ci$x_bar -  bf_ci$ci_lower)
  print(ratio_vec)

  graphics::lines(bf_ci$n,
                  (bf_ci$x_bar -  bf_ci$ci_lower) / (x_bar[bf_ci$n] - ci_fixed[bf_ci$n]),
                  col = 5)

}
