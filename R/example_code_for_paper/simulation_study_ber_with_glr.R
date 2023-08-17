# If stcp is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/stcp")

library(stcp)
# Bernoulli case ----
p_pre <- 0.5
p_post <- 0.6
max_sample <- 1000L
v <- 0 # Worst-case for SR type detector
ARL_target <- max_sample / 2
alpha <- 1 / ARL_target
num_repeat <- 2000L

# CUSUM functions
run_ber_cusum <- function(x_vec, p_pre, p_post, thres) {
  n_max <- length(x_vec)
  n_star <- Inf
  m <- 0
  for (i in 1:n_max) {
    x <- x_vec[i]
    # Update the change statistic
    m <- max(m, 0) + ifelse(x == 1, log(p_post/p_pre), log((1-p_post)/(1-p_pre)))
    # Check whether the stopping time happens
    if (m > thres) {
      n_star <- i
      break
    }
  }
  return(n_star)
}

run_ber_cusum_with_memory <- function(x_vec, p_pre, p_post, thres) {
  n_max <- length(x_vec)
  n_star <- Inf
  m <- 0
  m_vec <- length(n_max)
  for (i in 1:n_max) {
    x <- x_vec[i]
    # Update the change statistic
    m <- max(m, 0) + ifelse(x == 1, log(p_post/p_pre), log((1-p_post)/(1-p_pre)))
    m_vec[i] <- m
    # Check whether the stopping time happens
    if (m > thres) {
      n_star <- i
    }
  }
  return(list(n_star = n_star, m_vec = m_vec))
}


run_ber_glr_cusum <- function(x_vec, p_pre, p_post, thres, window_size = 20L) {
  n_max <- length(x_vec)
  n_star <- Inf
  m <- 0
  p_hat <- 0
  x_history <- rep(0.5, window_size)
  window_index <- -1
  for (i in 1:n_max) {
    x <- x_vec[i]
    m <- max(m, 0) + ifelse(x == 1, log(p_hat/p_pre), log((1-p_hat)/(1-p_pre)))
    # Check whether the stopping time happens
    window_index <- (window_index + 1) %% 10
    x_history[window_index + 1] <- x
    # Update the change statistic
    p_hat <- mean(x_history)
    if (m > thres) {
      n_star <- i
      break
    }
  }
  return(n_star)
}

run_ber_glr_cusum_with_memory <- function(x_vec, p_pre, p_post, thres, window_size = 20L) {
  n_max <- length(x_vec)
  n_star <- Inf
  m <- 0
  p_hat <- p_post
  x_history <- rep(p_post, window_size)
  window_index <- -1
  m_vec <- length(n_max)
  p_hat_vec <- length(n_max)
  for (i in 1:n_max) {
    x <- x_vec[i]
    m <- max(m, 0) + ifelse(x == 1, log(p_hat/p_pre), log((1-p_hat)/(1-p_pre)))
    m_vec[i] <- m
    window_index <- (window_index + 1) %% 10
    x_history[window_index + 1] <- x
    # Update the change statistic
    p_hat <- mean(x_history)
    p_hat_vec[i] <- p_hat
    # Check whether the stopping time happens
    if (m > thres) {
      n_star <- i
    }
  }
  return(list(n_star = n_star, m_vec = m_vec, p_hat_vec = p_hat_vec))
}

# E-detectors
# Compute optimal delta star
delta_star <- (p_post - p_pre)
delta_upper <- 0.99 - p_pre
delta_lower <- 0.01

# Other settings
psi_fn_list <- generate_sub_B_fn(p = p_pre)
v_min <- 1
k_max <- 1e+3

# Build CP detectors
# When delta_lower = delta_upper = delta_star
stcp_star <- build_stcp_exp(
  alpha,
  p_pre,
  delta_star,
  delta_star,
  is_test = FALSE,
  psi_fn_list,
  s_fn = function(x) {
    x - p_pre
  },
  v_fn = function(x) {
    1
  },
  v_min,
  k_max,
  tol = 1e-6
)

# When delta_lower < delta_star < delta_upper
stcp_mix <- build_stcp_exp(
  alpha,
  p_pre,
  delta_lower,
  delta_upper,
  is_test = FALSE,
  psi_fn_list,
  s_fn = function(x) {
    x - p_pre
  },
  v_fn = function(x) {
    1
  },
  v_min,
  k_max,
  tol = 1e-6
)

fn_factor_run_ber_e_detector <- function(stcp_obj, is_SR_type) {
  run_ber_e_detector <- function(x_vec, p_pre, p_post, thres) {
    n_max <- length(x_vec)
    n_star <- Inf
    m <- 0
    for (i in 1:n_max) {
      x <- x_vec[i]
      # Update the change statistic
      m <- update_log_mix_e(
        x,
        stcp_obj$omega,
        stcp_obj$log_base_fn_list,
        stcp_obj$log_e_vec,
        stcp_obj$is_test,
        is_SR_type
      )
      stcp_obj$log_e_vec <- m$last_log_e_vec
      stcp_obj$n <- stcp_obj$n + 1
      # Check whether the stopping time happens
      if (m$log_mix_e_vec > stcp_obj$log_one_over_alpha) {
        n_star <- i
        break
      }
    }
    return(n_star)
  }
  return(run_ber_e_detector)
}

run_ber_eSR_star <- fn_factor_run_ber_e_detector(stcp_star, is_SR_type = TRUE)
run_ber_eCS_star <- fn_factor_run_ber_e_detector(stcp_star, is_SR_type = FALSE)
run_ber_eSR_mix <- fn_factor_run_ber_e_detector(stcp_mix, is_SR_type = TRUE)
run_ber_eCS_mix <- fn_factor_run_ber_e_detector(stcp_mix, is_SR_type = FALSE)


# Simulation functions
WAD_ber_fn <- function(p_pre, p_post, thres, v, max_sample, num_repeat,
                          run_fn = run_ber_cusum) {
  single_run <- function() {
    x_vec <- c(rbinom(v, 1, p_pre), rbinom(max_sample - v, 1, p_post))
    stopped_time <- run_fn(x_vec, p_pre, p_post, thres)
    run_length <- ifelse(is.infinite(stopped_time), max_sample, stopped_time)
    return(run_length)
  }
  simulation_with_fixed_thres <- replicate(num_repeat, single_run())
  return(mean(simulation_with_fixed_thres[simulation_with_fixed_thres > v]) - v)
}

ARL_ber_fn <- function(p_pre, p_post, thres, max_sample, num_repeat,
                          run_fn = run_ber_cusum) {
  single_run <- function() {
    x_vec <- rbinom(max_sample, 1, p_pre)
    stopped_time <- run_fn(x_vec, p_pre, p_post, thres)
    run_length <- ifelse(is.infinite(stopped_time), max_sample, stopped_time)
    return(run_length)
  }
  simulation_with_fixed_thres <- replicate(num_repeat, single_run())
  return(mean(simulation_with_fixed_thres))
}

find_exact_thres <- function(p_pre, p_post, ARL_target, max_sample, num_repeat,
                             run_fn = run_ber_cusum) {
  f <- function(thres) {
    arl <- ARL_ber_fn(p_pre, p_post, thres, max_sample, num_repeat, run_fn)
    return(arl - ARL_target * 1.01) # 1.01 factor is a butter to ensure the final ALR > target.
  }
  exact_thres <- stats::uniroot(f,
                             c(1, log(ARL_target) * 10), tol = 1e-6)
  return(exact_thres$root)
}


# 1. CUSUM with fixed threshold log(ARL target)
set.seed(100)
ARL_fixed_cusum <- ARL_ber_fn(p_pre, p_post,
                                 thres = log(ARL_target),
                                 max_sample = max_sample,
                                 num_repeat = num_repeat)
WAD_fixed_cusum <- WAD_ber_fn(p_pre, p_post, v,
                                 thres = log(ARL_target),
                                 max_sample = max_sample,
                                 num_repeat = num_repeat)
# 2. Exact CUSUM.
set.seed(100)
exact_cusum_thres <- find_exact_thres(p_pre, p_post,
                                      ARL_target, max_sample, num_repeat)
ARL_exact_cusum <- ARL_ber_fn(p_pre, p_post,
                                 thres = exact_cusum_thres,
                                 max_sample = max_sample,
                                 num_repeat = num_repeat)
WAD_exact_cusum <- WAD_ber_fn(p_pre, p_post, v,
                                 thres = exact_cusum_thres,
                                 max_sample = max_sample,
                                 num_repeat = num_repeat)

# 3. Window-limited GLR CUSUM with fixed threshold log(ARL target)
set.seed(100)
ARL_glr_cusum <- ARL_ber_fn(p_pre, p_post,
                              thres = log(ARL_target),
                              max_sample = max_sample,
                              num_repeat = num_repeat,
                              run_fn = run_ber_glr_cusum)
WAD_glr_cusum <- WAD_ber_fn(p_pre, p_post, v,
                              thres = log(ARL_target),
                              max_sample = max_sample,
                              num_repeat = num_repeat,
                              run_fn = run_ber_glr_cusum)

# 4. Window-limited GLR CUSUM with optimal threshold
set.seed(100)
exact_glr_cusum_thres <- find_exact_thres(p_pre, p_post,
                                          ARL_target, max_sample, num_repeat,
                                          run_fn = run_ber_glr_cusum)
ARL_exact_glr_cusum <- ARL_ber_fn(p_pre, p_post,
                              thres = exact_glr_cusum_thres,
                              max_sample = max_sample,
                              num_repeat = num_repeat,
                              run_fn = run_ber_glr_cusum)
WAD_exact_glr_cusum <- WAD_ber_fn(p_pre, p_post, v,
                              thres = exact_glr_cusum_thres,
                              max_sample = max_sample,
                              num_repeat = num_repeat,
                              run_fn = run_ber_glr_cusum)

# 5. e-SR with known p_post
set.seed(100)
ARL_eSR_star <- ARL_ber_fn(p_pre, p_post,
                            thres = log(ARL_target),
                            max_sample = max_sample,
                            num_repeat = num_repeat,
                            run_fn = run_ber_eSR_star)
WAD_eSR_star <- WAD_ber_fn(p_pre, p_post, v,
                            thres = log(ARL_target),
                            max_sample = max_sample,
                            num_repeat = num_repeat,
                            run_fn = run_ber_eSR_star)

# 5. e-SR with unknown p_post
set.seed(100)
ARL_eSR_mix <- ARL_ber_fn(p_pre, p_post,
                           thres = log(ARL_target),
                           max_sample = max_sample,
                           num_repeat = num_repeat,
                           run_fn = run_ber_eSR_mix)
WAD_eSR_mix <- WAD_ber_fn(p_pre, p_post, v,
                           thres = log(ARL_target),
                           max_sample = max_sample,
                           num_repeat = num_repeat,
                           run_fn = run_ber_eSR_mix)

#6. e-CUSUM with known p_post
set.seed(100)
ARL_eCS_star <- ARL_ber_fn(p_pre, p_post,
                           thres = log(ARL_target),
                           max_sample = max_sample,
                           num_repeat = num_repeat,
                           run_fn = run_ber_eCS_star)
WAD_eCS_star <- WAD_ber_fn(p_pre, p_post, v,
                           thres = log(ARL_target),
                           max_sample = max_sample,
                           num_repeat = num_repeat,
                           run_fn = run_ber_eCS_star)




# Single run visualization
# No change
x_vec <- rbinom(max_sample, 1, p_pre)
exact_cusum <- run_ber_cusum_with_memory(x_vec, p_pre, p_post,
                                         thres = exact_cusum_thres)
exact_glr_cusum <- run_ber_glr_cusum_with_memory(x_vec, p_pre, p_post,
                                         thres = exact_glr_cusum_thres,
                                         window_size = 30L)
eSR_star <- run_stcp(x_vec, stcp_star)
eSR_mix <- run_stcp(x_vec, stcp_mix)

plot(seq_along(x_vec), exact_cusum$m_vec, type = "l", ylim = c(0, 3 * log(ARL_target)))
abline(h=exact_cusum_thres)
lines(seq_along(x_vec), exact_glr_cusum$m_vec, col = 2)
abline(h=exact_glr_cusum_thres, col = 2)
lines(seq_along(x_vec), eSR_star$log_mix_e_vec, col = 3)
abline(h=eSR_star$stcp_obj$log_one_over_alpha, col = 3)
lines(seq_along(x_vec), eSR_mix$log_mix_e_vec, col = 4)
abline(h=eSR_mix$stcp_obj$log_one_over_alpha, col = 4)


plot(seq_along(x_vec), eSR_star$log_mix_e_vec - exact_cusum$m_vec, type = "l", ylim = c(0, 10))
abline(h = eSR_mix$stcp_obj$log_one_over_alpha - exact_cusum_thres)


plot(seq_along(x_vec), eSR_mix$log_mix_e_vec - exact_cusum$m_vec, type = "l", ylim = c(0, 10))
abline(h = eSR_mix$stcp_obj$log_one_over_alpha - exact_cusum_thres)

plot(seq_along(x_vec), eSR_mix$log_mix_e_vec - exact_glr_cusum$m_vec, type = "l", ylim = c(0, 10))
abline(h = eSR_mix$stcp_obj$log_one_over_alpha - exact_glr_cusum_thres)

# plot(seq_along(x_vec), exact_glr_cusum$p_hat, type = "l")

# Change at v = 0 (Worst case)
v <- 0
x_vec <- c(rbinom(v, 1, p_pre), rbinom(max_sample - v, 1, p_post))
exact_cusum <- run_ber_cusum_with_memory(x_vec, p_pre, p_post,
                                         thres = exact_cusum_thres)
exact_glr_cusum <- run_ber_glr_cusum_with_memory(x_vec, p_pre, p_post,
                                                 thres = exact_glr_cusum_thres,
                                                 window_size = 30L)
eSR_star <- run_stcp(x_vec, stcp_star)
eSR_mix <- run_stcp(x_vec, stcp_mix)

plot(seq_along(x_vec), exact_cusum$m_vec, type = "l", ylim = c(0, 3 * log(ARL_target)))
abline(h=exact_cusum_thres)
lines(seq_along(x_vec), exact_glr_cusum$m_vec, col = 2)
abline(h=exact_glr_cusum_thres, col = 2)
lines(seq_along(x_vec), eSR_star$log_mix_e_vec, col = 3)
abline(h=eSR_star$stcp_obj$log_one_over_alpha, col = 3)
lines(seq_along(x_vec), eSR_mix$log_mix_e_vec, col = 4)
abline(h=eSR_mix$stcp_obj$log_one_over_alpha, col = 4)

plot(seq_along(x_vec), eSR_star$log_mix_e_vec - exact_cusum$m_vec, type = "l", ylim = c(0, 5))
abline(h = eSR_mix$stcp_obj$log_one_over_alpha - exact_cusum_thres)

plot(seq_along(x_vec), eSR_mix$log_mix_e_vec - exact_cusum$m_vec, type = "l", ylim = c(0, 5))
abline(h = eSR_mix$stcp_obj$log_one_over_alpha - exact_cusum_thres)

plot(seq_along(x_vec), eSR_mix$log_mix_e_vec - exact_glr_cusum$m_vec, type = "l", ylim = c(0, 5))
abline(h = eSR_mix$stcp_obj$log_one_over_alpha - exact_glr_cusum_thres)

# Change at v = 200 Favorable for SR-type detector
v <- 200
x_vec <- c(rbinom(v, 1, p_pre), rbinom(max_sample - v, 1, p_post))
exact_cusum <- run_ber_cusum_with_memory(x_vec, p_pre, p_post,
                                         thres = exact_cusum_thres)
exact_glr_cusum <- run_ber_glr_cusum_with_memory(x_vec, p_pre, p_post,
                                                 thres = exact_glr_cusum_thres,
                                                 window_size = 30L)
eSR_star <- run_stcp(x_vec, stcp_star)
eSR_mix <- run_stcp(x_vec, stcp_mix)

plot(seq_along(x_vec), exact_cusum$m_vec, type = "l", ylim = c(0, 3 * log(ARL_target)))
abline(h=exact_cusum_thres)
lines(seq_along(x_vec), exact_glr_cusum$m_vec, col = 2)
abline(h=exact_glr_cusum_thres, col = 2)
lines(seq_along(x_vec), eSR_star$log_mix_e_vec, col = 3)
abline(h=eSR_star$stcp_obj$log_one_over_alpha, col = 3)
lines(seq_along(x_vec), eSR_mix$log_mix_e_vec, col = 4)
abline(h=eSR_mix$stcp_obj$log_one_over_alpha, col = 4)


plot(seq_along(x_vec), eSR_star$log_mix_e_vec - exact_cusum$m_vec, type = "l", ylim = c(0, 5))
abline(h = eSR_mix$stcp_obj$log_one_over_alpha - exact_cusum_thres)


plot(seq_along(x_vec), eSR_mix$log_mix_e_vec - exact_cusum$m_vec, type = "l", ylim = c(0, 5))
abline(h = eSR_mix$stcp_obj$log_one_over_alpha - exact_cusum_thres)

plot(seq_along(x_vec), eSR_mix$log_mix_e_vec - exact_glr_cusum$m_vec, type = "l", ylim = c(0, 5))
abline(h = eSR_mix$stcp_obj$log_one_over_alpha - exact_glr_cusum_thres)

# plot(seq_along(x_vec), exact_glr_cusum$p_hat, type = "l")
