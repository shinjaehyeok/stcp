test_that("CI and BF-CI are consistent for sub-G", {
  # N(mu, 1)
  max_sample <- 100L
  mu <- -1
  sig <- 5
  alpha <- 0.025
  v_min <- 1
  k_max <- 1e+3
  psi_fn_list_generator <- function() {
    generate_sub_G_fn(sig)
  }

  # Generate data
  x_vec <- rnorm(max_sample, mu, sig)
  x_bar <- cumsum(x_vec) / seq_along(x_vec)

  # Build CI
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


  bf_ci_model <- build_bf_ci_exp(alpha,
                                 n_upper,
                                 n_lower,
                                 psi_fn_list_generator,
                                 v_min = v_min,
                                 k_max = k_max)

  bf_ci <- compute_bf_ci(x_vec,
                         bf_ci_model,
                         max_num_ci = 10)


  ratio_vec <- ci_mix$ci_lower[bf_ci$n] / bf_ci$ci_lower

  expect_true(abs(mean(ratio_vec) - 1) < 1e-6)
  expect_true(var(ratio_vec) < 1e-6)
})
