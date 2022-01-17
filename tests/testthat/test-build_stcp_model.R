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
  alpha <- 0.2

  # Bounded
  m_pre <- 0.2 # Upper bound of mean of the null distribution.
  delta_lower <- 0.1  # Guess on the minimum gap between the null and alternative means
  delta_upper <- 0.2
  bounded_test_model <- build_stcp_bounded(alpha, m_pre, delta_lower, delta_upper, is_test = TRUE)

  is_rejected <- function(){
    x_vec <- rbeta(100, 2, 8) # H_0 dist: B(2, 8)
    out <- run_stcp(x_vec, bounded_test_model)
    return(ifelse(is.infinite(out$stopped_ind), FALSE, TRUE))
  }

  simul <- replicate(n_rep, is_rejected())
  mean(simul)
  expect_true(mean(simul) <= alpha)

  # sub_G
  m_pre <- 0 # Upper bound of mean of the null distribution.
  delta_lower <- 1  # Guess on the minimum gap between the null and alternative means
  delta_upper <- 2
  sub_G_test_model <- build_stcp_exp(alpha, m_pre,
                                       delta_lower, delta_upper,
                                       is_test = TRUE)

  is_rejected <- function(){
    x_vec <- rnorm(100) # H_0 dist: B(2, 8)
    out <- run_stcp(x_vec, sub_G_test_model)
    return(ifelse(is.infinite(out$stopped_ind), FALSE, TRUE))
  }

  simul <- replicate(n_rep, is_rejected())
  mean(simul)
  expect_true(mean(simul) <= alpha)

  # sub_B
  m_pre <- 0.2 # Upper bound of mean of the null distribution.
  delta_lower <- 0.1  # Guess on the minimum gap between the null and alternative means
  delta_upper <- 0.2
  sub_G_test_model <- build_stcp_exp(alpha, m_pre,
                                     delta_lower, delta_upper,
                                     is_test = TRUE,
                                     psi_fn_list = generate_sub_B_fn(m_pre))

  is_rejected <- function(){
    x_vec <- rbinom(100,1,m_pre) # H_0 dist: B(2, 8)
    out <- run_stcp(x_vec, sub_G_test_model)
    return(ifelse(is.infinite(out$stopped_ind), FALSE, TRUE))
  }

  simul <- replicate(n_rep, is_rejected())
  mean(simul)
  expect_true(mean(simul) <= alpha)


})

