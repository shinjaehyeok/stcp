# If stcp is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/stcp")

library(stcp)

# Bounded case ----
# Beta distribution is used for this simulation
# Pre-change : Beta(1, 4) (m = 0.2, sig^2 = 4 / (5^2 * 6) = 2/75)
# Post-change: Beta(2, 3) (m = 0.4, sig^2 = 6 / (5^2 * 6) = 2/50)

bounded_simulator <- function(a_pre = 1,
                              b_pre = 4,
                              a_post = 2,
                              b_post = 3,
                              max_sample = 200L,
                              v = max_sample / 10,
                              m_bound = a_pre / (a_pre + b_pre),
                              alpha = 10 / max_sample,
                              range_factor = 10,
                              num_repeat = 100L) {
  # Check whether CP is for ARL or WAD simulation
  if (v >= max_sample) {
    is_ARL <- TRUE
  } else {
    is_ARL <- FALSE
  }

  # compute mean and variance
  m_pre <- a_pre / (a_pre + b_pre)
  m_post <- a_post / (a_post + b_post)
  v_pre <- a_pre * b_pre / (a_pre + b_pre)^2 / (a_pre + b_pre + 1)
  v_post <- a_post * b_post / (a_post + b_post)^2 / (a_post + b_post + 1)

  # Compute optimal delta star
  delta_star <- (m_post - m_pre)
  var_star <- v_post

  # delta range = delta_star * / range_factor
  delta_upper <- delta_star * range_factor
  delta_upper <- min(c(delta_upper, 0.99 - m_pre)) # To make delta_upper + m_pre < 1
  delta_lower <- delta_star / range_factor

  # Build CP detectors
  # When delta_lower = delta_upper = delta_star
  stcp_star <- build_stcp_bounded(
    alpha,
    m_bound,
    delta_star,
    delta_star,
    is_test = FALSE,
    bound_lower = 0,
    bound_upper = 1,
    k_max = 1000,
    tol = 1e-6,
    var_lower = var_star,
    var_upper = var_star
  )

  # When delta_lower < delta_star < delta_upper
  stcp_mix <- build_stcp_bounded(
    alpha,
    m_bound,
    delta_lower,
    delta_upper,
    is_test = FALSE,
    bound_lower = 0,
    bound_upper = 1,
    k_max = 1000,
    tol = 1e-6
  )

  compute_single_run_stat <- function() {
    # Generate sample
    x_vec <- generator(
      max_sample,
      v,
      pre_sampler = function(n) {
        rbeta(n, a_pre, b_pre)
      },
      post_sampler = function(n) {
        rbeta(n, a_post, b_post)
      }
    )

    # Compute mixtures of SR and CUSUM e-detectors.
    star_SR_run <- run_stcp(x_vec, stcp_star)
    star_CS_run <- run_stcp(x_vec, stcp_star, is_SR_type = FALSE)

    mix_SR_run <- run_stcp(x_vec, stcp_mix)
    mix_CS_run <- run_stcp(x_vec, stcp_mix, is_SR_type = FALSE)


    # Stopping times of mixtures of e-SR and e-CUSUM procedures
    stopped_time <- c(
      star_SR = star_SR_run$stopped_ind,
      mix_SR = mix_SR_run$stopped_ind,
      star_CS = star_CS_run$stopped_ind,
      mix_CS = mix_CS_run$stopped_ind
    )

    ARL <- ifelse(is.infinite(stopped_time), max_sample, stopped_time)
    if (!is_ARL) {
      WAD <- ifelse(ARL > v, ARL - v, NA)
      return(WAD)
    } else {
      return(ARL)
    }
  }

  f <- function(iter) {
    result <- compute_single_run_stat()
    return(c(iter = iter, result))
  }
  multiple_run <- lapply(1:num_repeat, f)
  multiple_run <- do.call('rbind', multiple_run)
  summ <- colMeans(multiple_run, na.rm = TRUE)
  sd_summ <- apply(multiple_run, 2, sd, na.rm = TRUE)

  return(
    list(is_ARL = is_ARL,
         num_repeat = num_repeat,
         run_avg = summ[-1],
         run_sd = sd_summ[-1],
         run_raw = multiple_run,
         one_over_alpha = 1/alpha,
         m_bound = m_bound,
         max_sample = max_sample,
         v = v,
         m_pre = m_pre,
         m_post = m_post,
         a_pre = a_pre,
         b_pre = b_pre,
         a_post = a_post,
         b_post = b_post
    )
  )
}


# ARL study

m_pre <- 1 / 5
m_bound_vec <- m_pre + seq(-0.1, 0.1, by = 0.05)
max_sample <- 200L

ARL_fn <- function(m_bound){
  out <- bounded_simulator(a_pre = 1,
                           b_pre = 4,
                           a_post = 2,
                           b_post = 3,
                           max_sample = max_sample,
                           v = max_sample, # ARL study
                           m_bound = m_bound,
                           alpha = 10 / max_sample,
                           range_factor = 10,
                           num_repeat = 100L)

  return(c(
    m_bound = m_bound,
    ARL_level = out$one_over_alpha,
    ARL = out$run_avg,
    sd = out$run_sd
  ))
}

ARL_result <- lapply(m_bound_vec, ARL_fn)
ARL_result <- do.call('rbind', ARL_result)

# Worst average delay study
a_pre <- 1
b_pre <- 4
m_pre <- a_pre / (a_pre + b_pre)
m_post_vec <- m_pre + seq(0.02, 0.1, by = 0.02)
b_post <- 3

max_sample <- 200L

WAD_fn <- function(m_post){
  a_post <- m_post * b_post / (1-m_post)
  out <- bounded_simulator(a_pre = a_pre,
                       b_pre = b_pre,
                       a_post = a_post,
                       b_post = b_post,
                       max_sample = max_sample,
                       v = max_sample / 10, # WAD study
                       m_bound = a_pre / (a_pre + b_pre), # Correctly specified boundary
                       alpha = 1 / max_sample, # Set ARL >= max_sample
                       range_factor = 10,
                       num_repeat = 100L)

  return(c(
    m_pre = out$m_pre,
    m_post = out$m_post,
    WAD = out$run_avg,
    sd = out$run_sd
  ))
}


WAD_result <- lapply(m_post_vec, WAD_fn)
WAD_result <- do.call('rbind', WAD_result)
