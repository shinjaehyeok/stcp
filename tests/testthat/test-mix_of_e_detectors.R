test_that("update_log_mix_e and update_log_mix_e_values are matched for the univariate test", {
  alpha <- 0.1
  p_pre <- 0.3
  p_post <- sample(c(0.3, 0.5, 0.6), 1)
  delta_lower <- 0.01
  delta_upper <- 0.5

  new_x <-  rbinom(100, 1, p_post)

  # stcp_exp
  model_exp <- build_stcp_exp(
    alpha,
    p_pre,
    delta_lower,
    delta_upper,
    is_test = TRUE
  )

  normal_update <- update_log_mix_e(
    new_x,
    model_exp$omega,
    model_exp$log_base_fn_list,
    model_exp$log_e_vec,
    model_exp$is_test
  )

  faster_update <- update_log_mix_e_values(new_x,
                                           model_exp$omega,
                                           model_exp$log_base_fn_list,
                                           model_exp$log_e_vec)

  expect_equal(normal_update, faster_update)

  # stcp_bounded
  model_bounded <- build_stcp_bounded(
    alpha,
    p_pre,
    delta_lower,
    delta_upper,
    is_test = TRUE
  )

  normal_update_bounded <- update_log_mix_e(
    new_x,
    model_bounded$omega,
    model_bounded$log_base_fn_list,
    model_bounded$log_e_vec,
    model_bounded$is_test
  )

  faster_update_bounded <- update_log_mix_e_values(new_x,
                                           model_bounded$omega,
                                           model_bounded$log_base_fn_list,
                                           model_bounded$log_e_vec)

  expect_equal(normal_update_bounded, faster_update_bounded)
})
