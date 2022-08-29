compute_SR <- function(x){
  out <- numeric(length(x))
  sr <- 0
  for (i in seq_along(out)) {
    out[i] <- log((sr + 1)) + x[i]
    sr <- exp(out[i])
  }
  return(out)
}

compute_CS <- function(x){
  out <- numeric(length(x))
  sr <- 0
  for (i in seq_along(out)) {
    out[i] <- log(max(sr,1)) + x[i]
    sr <- exp(out[i])
  }
  return(out)
}

test_that("Consistent with the parametric cases 1. Gaussian", {
  # sub G vs Gaussian
  m_post <- 1
  m_pre <- -1
  delta_star <- m_post - m_pre

  lrt_G <- function(x){
    log(dnorm(x, mean = m_post) / dnorm(x, mean = m_pre))
  }

  sub_G_test <- build_stcp_exp(alpha = 0.1,
                               m_pre = m_pre,
                               delta_lower = delta_star,
                               delta_upper = delta_star,
                               is_test = TRUE)

  sub_G_cp <- build_stcp_exp(alpha = 0.1,
                             m_pre = m_pre,
                             delta_lower = delta_star,
                             delta_upper = delta_star,
                             is_test = FALSE)


  x_vec <- rnorm(100)
  lrt_G_out_test <- cumsum(lrt_G(x_vec))
  lrt_G_out_sr <- compute_SR(lrt_G(x_vec))
  lrt_G_out_cs <- compute_CS(lrt_G(x_vec))
  sub_G_out_test <- run_stcp(x_vec, sub_G_test)$log_mix_e_vec
  sub_G_out_sr <- run_stcp(x_vec, sub_G_cp)$log_mix_e_vec
  sub_G_out_cs <- run_stcp(x_vec, sub_G_cp, is_SR_type = FALSE)$log_mix_e_vec

  expect_true(
    mean(abs(lrt_G_out_test - sub_G_out_test)) < 1e-8
  )

  expect_true(
    mean(abs(lrt_G_out_sr - sub_G_out_sr)) < 1e-8
  )

  expect_true(
    mean(abs(lrt_G_out_cs - sub_G_out_cs)) < 1e-8
  )

})

test_that("Consistent with the parametric cases 2. Bernoulli", {
  # sub B vs Bernoulli
  p_post <- 0.5
  p_pre <- 0.2
  delta_star <- p_post - p_pre

  lrt_B <- function(x){
    log(dbinom(x, 1, p_post) / dbinom(x, 1, p_pre))
  }

  sub_B_test <- build_stcp_exp(alpha = 0.1,
                               m_pre = p_pre,
                               delta_lower = delta_star,
                               delta_upper = delta_star,
                               is_test = TRUE,
                               psi_fn_list = generate_sub_B_fn(p_pre))

  sub_B_cp <- build_stcp_exp(alpha = 0.1,
                             m_pre = p_pre,
                             delta_lower = delta_star,
                             delta_upper = delta_star,
                             is_test = FALSE,
                             psi_fn_list = generate_sub_B_fn(p_pre))

  x_vec <- rbinom(100, 1, 0.5)
  lrt_B_out_test <- cumsum(lrt_B(x_vec))
  lrt_B_out_sr <- compute_SR(lrt_B(x_vec))
  lrt_B_out_cs <- compute_CS(lrt_B(x_vec))
  sub_B_out_test <- run_stcp(x_vec, sub_B_test)$log_mix_e_vec
  sub_B_out_sr <- run_stcp(x_vec, sub_B_cp)$log_mix_e_vec
  sub_B_out_cs <- run_stcp(x_vec, sub_B_cp, is_SR_type = FALSE)$log_mix_e_vec

  expect_true(
    mean(abs(lrt_B_out_test - sub_B_out_test)) < 1e-8
  )

  expect_true(
    mean(abs(lrt_B_out_sr - sub_B_out_sr)) < 1e-8
  )

  expect_true(
    mean(abs(lrt_B_out_cs - sub_B_out_cs)) < 1e-8
  )


})


test_that("Type 1 error control", {
  n_rep <- 1000
  x_len <- 100
  alpha <- 0.1

  # Bounded
  # Get proper delta for test of length x_len
  # This step is not for actual test but just for a tighter simulation
  baseline_obj <- compute_baseline_for_sample_size(alpha,
                                                   n_upper = x_len,
                                                   n_lower = x_len / 10,
                                                   psi_fn_list = generate_sub_E_fn(),
                                                   v_min = 0)

  m_pre <- 0.2 # Upper bound of mean of the null distribution.
  delta_lower_sub_E <- baseline_obj$delta_lower
  delta_upper_sub_E <- baseline_obj$delta_upper
  bounded_test_model <- build_stcp_bounded(alpha,
                                           m_pre,
                                           delta_lower = 0.1,
                                           is_test = TRUE,
                                           delta_lower_sub_E = delta_lower_sub_E,
                                           delta_upper_sub_E = delta_upper_sub_E)

  is_rejected <- function(){
    x_vec <- rbeta(x_len, 2, 8) # H_0 dist: B(2, 8)
    out <- run_stcp(x_vec, bounded_test_model)
    return(ifelse(is.infinite(out$stopped_ind), FALSE, TRUE))
  }

  simul <- replicate(n_rep, is_rejected())
  expect_true(mean(simul) <= alpha)

  # Check Bounded delta lower and upper were recovered properly
  bounded_test_model2 <- build_stcp_bounded(
    alpha,
    m_pre,
    delta_lower = bounded_test_model$delta_lower,
    delta_upper = bounded_test_model$delta_upper,
    is_test = TRUE,
    bound_lower = bounded_test_model$bound_lower,
    bound_upper = bounded_test_model$bound_upper,
    var_lower = bounded_test_model$var_lower,
    var_upper = bounded_test_model$var_upper)

  expect_true(
    mean(abs(bounded_test_model$lambda - bounded_test_model2$lambda)) < 1e-8
  )

  # sub_G
  # Get proper delta for test of length x_len
  # This step is not for actual test but just for a tighter simulation
  baseline_obj <- compute_baseline_for_sample_size(alpha,
                                                   n_upper = x_len * 2,
                                                   n_lower = x_len / 10,
                                                   psi_fn_list = generate_sub_G_fn())

  m_pre <- 0 # Upper bound of mean of the null distribution.
  delta_lower <- baseline_obj$delta_lower
  delta_upper <- baseline_obj$delta_upper
  sub_G_test_model <- build_stcp_exp(alpha, m_pre,
                                     delta_lower, delta_upper,
                                     is_test = TRUE)

  is_rejected <- function(){
    x_vec <- rnorm(x_len) # H_0 dist: B(2, 8)
    out <- run_stcp(x_vec, sub_G_test_model)
    return(ifelse(is.infinite(out$stopped_ind), FALSE, TRUE))
  }

  simul <- replicate(n_rep, is_rejected())
  expect_true(mean(simul) <= alpha)

  # sub_B
  m_pre <- 0.2 # Upper bound of mean of the null distribution.

  # Get proper delta for test of length x_len
  # This step is not for actual test but just for a tighter simulation
  baseline_obj <- compute_baseline_for_sample_size(alpha,
                                                   n_upper = x_len,
                                                   n_lower = x_len / 10,
                                                   psi_fn_list = generate_sub_B_fn(m_pre))


  delta_lower <- baseline_obj$delta_lower
  delta_upper <- baseline_obj$delta_upper
  sub_B_test_model <- build_stcp_exp(alpha, m_pre,
                                     delta_lower, delta_upper,
                                     is_test = TRUE,
                                     psi_fn_list = generate_sub_B_fn(m_pre))

  is_rejected <- function(){
    x_vec <- rbinom(x_len,1,m_pre)
    out <- run_stcp(x_vec, sub_B_test_model)
    return(ifelse(is.infinite(out$stopped_ind), FALSE, TRUE))
  }

  simul <- replicate(n_rep, is_rejected())
  expect_true(mean(simul) <= alpha)


})

test_that("is_flipped automatically converts the bounded input", {
  alpha <- 1e-2 # Level 0.01 test
  m_pre <- 0.2 # Upper bound of mean of the null distribution.
  delta_lower <- 0.01  # Guess on the minimum gap between the null and alternative means

  max_sample <- 1000L
  x_vec <- rbeta(max_sample, 2, 8) # H_0 dist: B(2, 8)


  # Run manually flipped model (Y = 1-X, theta = 1-m)
  theta_pre <- 1-m_pre # Upper bound of mean of the null distribution.
  manual_flip_model <- build_stcp_bounded(alpha, theta_pre, delta_lower, is_test = TRUE)
  manual_flip_run <- run_stcp(1-x_vec, manual_flip_model) # Note input is 1-x_vec

  # Run automatically flipped model by enabling is_flipped = TRUE
  flip_model <- build_stcp_bounded(alpha, m_pre, delta_lower, is_test = TRUE, is_flipped = TRUE)
  flip_run <- run_stcp(x_vec, flip_model)

  # 1. Models are the same
  expect_true(all.equal(manual_flip_model$omega, flip_model$omega))
  expect_true(all.equal(manual_flip_model$lambda, flip_model$lambda))

  # 2. Runs are the same
  expect_true(all.equal(manual_flip_run$log_mix_e_vec, flip_run$log_mix_e_vec))
})

test_that("Model combinataion works", {
  # 1. Bounded case

  max_sample <- 1000L
  x_vec <- 3*rbeta(max_sample, 2, 8) + 1 # H_0 dist: B(2, 8)
  bound_lower <- 1
  bound_upper <- 4
  alpha <- 1e-2 # Level 0.01 test
  m_pre <- 0.2 * 3 + 1 # Upper bound of mean of the null distribution.
  delta_lower <- 0.01  # Guess on the minimum gap between the null and alternative means

  normal_model <- build_stcp_bounded(alpha, m_pre, delta_lower, is_test = TRUE, bound_lower = bound_lower, bound_upper = bound_upper)
  flip_model <- build_stcp_bounded(alpha, m_pre, delta_lower, is_test = TRUE, is_flipped = TRUE, bound_lower = bound_lower, bound_upper = bound_upper)

  combined_model <- combine_stcp(normal_model, flip_model, w = 0.5)
  combined_model2 <- combine_stcp(normal_model, flip_model, w = 0.9)

  normal_run <- run_stcp(x_vec, normal_model)
  flip_run <- run_stcp(x_vec, flip_model)

  combine_run <- combine_stcp_run(normal_run, flip_run, w = 0.5)
  combine_run2 <- combine_stcp_run(normal_run, flip_run, w = 0.9)

  combined_model_run <- run_stcp(x_vec, combined_model)
  combined_model_run2 <- run_stcp(x_vec, combined_model2)

  expect_true(all.equal(combine_run$log_mix_e_vec, combined_model_run$log_mix_e_vec))
  expect_true(all.equal(combine_run2$log_mix_e_vec, combined_model_run2$log_mix_e_vec))
  expect_false(isTRUE(all.equal(combine_run$log_mix_e_vec, combined_model_run2$log_mix_e_vec)))

  # 2. Gaussian case (WIP)

  # 3. Bernoulli case (WIP)
})
