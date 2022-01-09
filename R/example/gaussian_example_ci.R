# If EDCP is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/EDCP")

library(EDCP)

# Gaussian case ----
# N(mu, 1)
max_sample <- 1000L
mu <- -1
sig <- 5
alpha <- 0.025
v_min <- 1
k_max <- 1e+3
psi_fn_list_generator <- function() {
  generate_sub_G_fn(sig)
}

check_bf_test <- FALSE

# Generate data
x_vec <- rnorm(max_sample, mu, sig)
x_bar <- cumsum(x_vec) / seq_along(x_vec)

# Build CI
# When delta_lower = delta_upper = delta_star
# Compute optimal delta star
n_target <- round(max_sample / 2)
ci_model_star <- build_ci_exp(alpha,
                              n_target,
                              n_target,
                              psi_fn_list_generator,
                              v_min = v_min,
                              k_max = k_max)

ci_star <- compute_ci(x_vec, ci_model_star, width_upper = sig * 100)


# When delta_lower < delta_star < delta_upper
# Compute target interval
n_lower <- max_sample / 5
n_upper <- max_sample * 2
ci_model_mix <- build_ci_exp(alpha,
                             n_upper,
                             n_lower,
                             psi_fn_list_generator,
                             v_min = v_min,
                             k_max = k_max)

ci_mix <- compute_ci(x_vec, ci_model_mix, width_upper = sig * 100)

ci_fixed <-
  x_bar - sig * qnorm(alpha, lower.tail = FALSE) / sqrt(1:max_sample)

# Plot CI

plot(
  1:max_sample,
  x_bar,
  type = "l",
  xlab = "n",
  ylab = "X_bar",
  main = "Running average",
  ylim = c(-0.5, 0.5) * sig + mu,
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
  ci_model <- ci_model_mix
  bruth_force_ci_helper <- function(m, x_vec) {
    edcp_model <- build_edcp_exp(
      alpha,
      m,
      ci_model$baseline_obj$delta_lower,
      ci_model$baseline_obj$delta_upper,
      is_test = TRUE,
      ci_model$baseline_obj$psi_fn_list,
      s_fn = function(x) {
        x - m
      },
      v_min = v_min,
      k_max = k_max
    )
    e_val <- run_edcp(x_vec, edcp_model)
    return(e_val$log_mix_e_vec[length(e_val$log_mix_e_vec)] - e_val$edcp_obj$log_one_over_alpha)
  }

  # Warning::This code is very slow O(n^2)
  n_bf_vec <- c(seq(1, length(x_vec), by = 10L), length(x_vec))
  ci_bf <-
    sapply(n_bf_vec, function(n) {
      compute_bf_ci_single(x_vec[1:n],
                             bruth_force_ci_helper,
                             width_upper = 100)
    })

  ratio_vec <- ci_mix$ci_lower[n_bf_vec] / ci_bf
  if (abs(mean(ratio_vec) - 1) > 1e-8)
    stop("BF test failed")
  if (var(ratio_vec) > 1e-8)
    stop("BF test failed")
}
