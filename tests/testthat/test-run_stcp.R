test_that("Sequential update works for SimpleExp test", {
  alpha <- 0.1
  p_pre <- 0.3
  p_post <- sample(c(0.3, 0.5, 0.6), 1)
  delta_lower <- 0.01
  delta_upper <- 0.5

  # SimpleExp
  x <-  rbinom(100, 1, p_post)
  model_sub_B_test <- build_stcp_exp(
    alpha,
    p_pre,
    delta_lower,
    delta_upper,
    is_test = TRUE,
    psi_fn_list = generate_sub_B_fn(p_pre)
  )

  # Batch run
  B_1_100 <- run_stcp(x, model_sub_B_test)

  # Sequential runs (1-50, 51, 52-100)
  B_1_50 <- run_stcp(x[1:50], model_sub_B_test)
  B_51 <- run_stcp(x[51], B_1_50$stcp_obj)
  B_52_100 <- run_stcp(x[52:100], B_51$stcp_obj)

  # log e-values are matched
  check_two_num_vec <- function(x, y, tol = 1e-8){
    expect_true(length(x) == length(y))
    expect_true(mean(abs(x-y)) <= tol)
  }
  check_two_num_vec(
    B_1_100$log_mix_e_vec,
    c(B_1_50$log_mix_e_vec,
      B_51$log_mix_e_vec,
      B_52_100$log_mix_e_vec)
    )

  # Stopped index matched
  if (is.infinite(B_1_100$stopped_ind)){
    expect_true(is.infinite(B_1_50$stopped_ind))
    expect_true(is.infinite(B_51$stopped_ind))
    expect_true(is.infinite(B_52_100$stopped_ind))
  } else if (B_1_100$stopped_ind <= 50) {
    expect_true(B_1_50$stopped_ind == B_1_100$stopped_ind)
  } else if (B_1_100$stopped_ind == 51) {
    expect_true(is.infinite(B_1_50$stopped_ind))
    expect_true(B_51$stopped_ind == 1)
  } else {
    expect_true(is.infinite(B_1_50$stopped_ind))
    expect_true(is.infinite(B_51$stopped_ind))
    expect_true(B_1_100$stopped_ind == (B_52_100$stopped_ind +
                                          B_52_100$stcp_obj$n -
                                          length(B_52_100$log_mix_e_vec)))
  }

})

test_that("Sequential update works for Bounded test", {
  alpha <- 0.1
  p_pre <- 0.3
  p_post <- sample(c(0.3, 0.5, 0.6), 1)
  delta_lower <- 0.01
  delta_upper <- 0.5

  # SimpleExp
  x <-  rbinom(100, 1, p_post)
  model_bounded_test <- build_stcp_bounded(
    alpha,
    p_pre,
    delta_lower,
    delta_upper,
    is_test = TRUE
  )

  # Batch run
  bound_1_100 <- run_stcp(x, model_bounded_test)

  # Sequential runs (1-50, 51, 52-100)
  bound_1_50 <- run_stcp(x[1:50], model_bounded_test)
  bound_51 <- run_stcp(x[51], bound_1_50$stcp_obj)
  bound_52_100 <- run_stcp(x[52:100], bound_51$stcp_obj)

  # log e-values are matched
  check_two_num_vec <- function(x, y, tol = 1e-8){
    expect_true(length(x) == length(y))
    expect_true(mean(abs(x-y)) <= tol)
  }
  check_two_num_vec(
    bound_1_100$log_mix_e_vec,
    c(bound_1_50$log_mix_e_vec,
      bound_51$log_mix_e_vec,
      bound_52_100$log_mix_e_vec)
  )

  # Stopped index matched
  if (is.infinite(bound_1_100$stopped_ind)){
    expect_true(is.infinite(bound_1_50$stopped_ind))
    expect_true(is.infinite(bound_51$stopped_ind))
    expect_true(is.infinite(bound_52_100$stopped_ind))
  } else if (bound_1_100$stopped_ind <= 50) {
    expect_true(bound_1_50$stopped_ind == bound_1_100$stopped_ind)
  } else if (bound_1_100$stopped_ind == 51) {
    expect_true(is.infinite(bound_1_50$stopped_ind))
    expect_true(bound_51$stopped_ind == 1)
  } else {
    expect_true(is.infinite(bound_1_50$stopped_ind))
    expect_true(is.infinite(bound_51$stopped_ind))
    expect_true(bound_1_100$stopped_ind == (bound_52_100$stopped_ind +
                                          bound_52_100$stcp_obj$n -
                                          length(bound_52_100$log_mix_e_vec)))
  }

})


test_that("Sequential update works for SimpleExp CP", {
  alpha <- 0.1
  p_pre <- 0.3
  p_post <- sample(c(0.3, 0.5, 0.6), 1)
  delta_lower <- 0.01
  delta_upper <- 0.5

  # SimpleExp
  x <-  rbinom(100, 1, p_post)
  model_sub_B_test <- build_stcp_exp(
    alpha,
    p_pre,
    delta_lower,
    delta_upper,
    is_test = FALSE,
    psi_fn_list = generate_sub_B_fn(p_pre)
  )

  # Batch run
  B_1_100 <- run_stcp(x, model_sub_B_test)

  # Sequential runs (1-50, 51, 52-100)
  B_1_50 <- run_stcp(x[1:50], model_sub_B_test)
  B_51 <- run_stcp(x[51], B_1_50$stcp_obj)
  B_52_100 <- run_stcp(x[52:100], B_51$stcp_obj)

  # log e-values are matched
  check_two_num_vec <- function(x, y, tol = 1e-8){
    expect_true(length(x) == length(y))
    expect_true(mean(abs(x-y)) <= tol)
  }
  check_two_num_vec(
    B_1_100$log_mix_e_vec,
    c(B_1_50$log_mix_e_vec,
      B_51$log_mix_e_vec,
      B_52_100$log_mix_e_vec)
  )

  # Stopped index matched
  if (is.infinite(B_1_100$stopped_ind)){
    expect_true(is.infinite(B_1_50$stopped_ind))
    expect_true(is.infinite(B_51$stopped_ind))
    expect_true(is.infinite(B_52_100$stopped_ind))
  } else if (B_1_100$stopped_ind <= 50) {
    expect_true(B_1_50$stopped_ind == B_1_100$stopped_ind)
  } else if (B_1_100$stopped_ind == 51) {
    expect_true(is.infinite(B_1_50$stopped_ind))
    expect_true(B_51$stopped_ind == 1)
  } else {
    expect_true(is.infinite(B_1_50$stopped_ind))
    expect_true(is.infinite(B_51$stopped_ind))
    expect_true(B_1_100$stopped_ind == (B_52_100$stopped_ind +
                                          B_52_100$stcp_obj$n -
                                          length(B_52_100$log_mix_e_vec)))
  }

})

test_that("Sequential update works for Bounded CP", {
  alpha <- 0.1
  p_pre <- 0.3
  p_post <- sample(c(0.3, 0.5, 0.6), 1)
  delta_lower <- 0.01
  delta_upper <- 0.5

  # SimpleExp
  x <-  rbinom(100, 1, p_post)
  model_bounded_test <- build_stcp_bounded(
    alpha,
    p_pre,
    delta_lower,
    delta_upper,
    is_test = FALSE
  )

  # Batch run
  bound_1_100 <- run_stcp(x, model_bounded_test)

  # Sequential runs (1-50, 51, 52-100)
  bound_1_50 <- run_stcp(x[1:50], model_bounded_test)
  bound_51 <- run_stcp(x[51], bound_1_50$stcp_obj)
  bound_52_100 <- run_stcp(x[52:100], bound_51$stcp_obj)

  # log e-values are matched
  check_two_num_vec <- function(x, y, tol = 1e-8){
    expect_true(length(x) == length(y))
    expect_true(mean(abs(x-y)) <= tol)
  }
  check_two_num_vec(
    bound_1_100$log_mix_e_vec,
    c(bound_1_50$log_mix_e_vec,
      bound_51$log_mix_e_vec,
      bound_52_100$log_mix_e_vec)
  )

  # Stopped index matched
  if (is.infinite(bound_1_100$stopped_ind)){
    expect_true(is.infinite(bound_1_50$stopped_ind))
    expect_true(is.infinite(bound_51$stopped_ind))
    expect_true(is.infinite(bound_52_100$stopped_ind))
  } else if (bound_1_100$stopped_ind <= 50) {
    expect_true(bound_1_50$stopped_ind == bound_1_100$stopped_ind)
  } else if (bound_1_100$stopped_ind == 51) {
    expect_true(is.infinite(bound_1_50$stopped_ind))
    expect_true(bound_51$stopped_ind == 1)
  } else {
    expect_true(is.infinite(bound_1_50$stopped_ind))
    expect_true(is.infinite(bound_51$stopped_ind))
    expect_true(bound_1_100$stopped_ind == (bound_52_100$stopped_ind +
                                              bound_52_100$stcp_obj$n -
                                              length(bound_52_100$log_mix_e_vec)))
  }

})


