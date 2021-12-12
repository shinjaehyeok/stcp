# If EDCP is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/EDCP")

library(EDCP)

# Bernoulli case ----
# Ber(p)

max_sample <- 1000L
pre_sample <- 10L
p <- 0.1
alpha <- 0.025

# For simplicity, we use some historical data to build psi function
# We may build more accurate CS by using the method in SGLRT paper.
p_hat_pre <- mean(rbinom(pre_sample, 1, p))
p_hat_pre <- min(0.999, max(p_hat_pre, 0.001))

psi_fn_list <- generate_sub_B_fn(p = p)
v_min <- 1
k_max <- 1e+3

check_bf_test <- TRUE

# Compute optimal delta star
n_target <- round(max_sample / 2)
delta_star <- convert_time_to_delta_bound(
  alpha,
  n_target,
  n_target,
  psi_fn_list,
  v_min,
  k_max
)

# Compute target interval
n_lower <- max_sample / 10
n_upper <- max_sample * 2
delta_mix <- convert_time_to_delta_bound(
  alpha,
  n_upper,
  n_lower,
  psi_fn_list,
  v_min,
  k_max
)

# Generate data
x_vec <-  rbinom(max_sample, 1, p)
x_bar <- cumsum(x_vec) / seq_along(x_vec)

plot(
  1:max_sample,
  x_bar,
  type = "l",
  xlab = "n",
  ylab = "X_bar",
  main = "Running average",
  ylim = c(0,1),
)
graphics::abline(h = p,
                 col = 2,
                 lwd = 2)

# Build CI
# When delta_lower = delta_upper = delta_star
baseline_star <- compute_baseline_simple_exp(alpha,
                                             0,
                                             delta_star$delta_lower,
                                             delta_star$delta_upper,
                                             psi_fn_list,
                                             v_min = v_min,
                                             k_max = k_max)

ci_star <- compute_ci_width(x_vec,
                            baseline_star$ci_helper)

# When delta_lower < delta_star < delta_upper
baseline_mix <- compute_baseline_simple_exp(alpha,
                                            0,
                                            delta_mix$delta_lower,
                                            delta_mix$delta_upper,
                                            psi_fn_list,
                                            v_min = v_min,
                                            k_max = k_max)

ci_mix <- compute_ci_width(x_vec,
                           baseline_mix$ci_helper)

ci_fixed <- sqrt(p * (1-p)) * qnorm(alpha, lower.tail = FALSE) / sqrt(1:max_sample)

# Plot CI
graphics::lines(1:max_sample,
                ci_star$x_bar - ci_fixed,
                lty = 2,
                col = 1)
graphics::lines(1:max_sample,
                ci_star$x_bar - ci_star$ci_width,
                lty = 2,
                col = 2)
graphics::lines(1:max_sample,
                ci_mix$x_bar - ci_mix$ci_width,
                lty = 2,
                col = 3)

# Compare width
plot(1:max_sample,
     ci_star$ci_width / ci_fixed,
     type = "l",
     ylim = c(1, 2))
graphics::lines(1:max_sample, ci_mix$ci_width / ci_fixed, col = 2)

# Check correctness via brute-force method
if (check_bf_test) {
  bruth_force_ci_helper <- function(m, x_vec) {
    baseline_m <- compute_baseline_simple_exp(
      alpha,
      m,
      delta_mix$delta_lower,
      delta_mix$delta_upper,
      psi_fn_list,
      s_fn = function(x) {
        x - m
      },
      v_min = v_min,
      k_max = k_max
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

  ratio_vec <- (ci_mix$x_bar - ci_mix$ci_width)[n_bf_vec] / ci_bf
  if (abs(mean(ratio_vec) - 1) > 1e-8) stop("BF test failed")
  if (  var(ratio_vec) > 1e-8) stop("BF test failed")
}
