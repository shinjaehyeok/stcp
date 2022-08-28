# If stcp is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/stcp")

library(stcp)

# Gaussian case ----

gaussian_simulator <- function(m_pre = 0 ,
                               sig_pre = 1,
                               m_post = 0.5,
                               sig_post = 1,
                               max_sample = 200L,
                               v = max_sample / 10,
                               m_bound = 0,
                               alpha = 10 / max_sample,
                               range_factor = 10,
                               num_repeat = 100L) {
  # Check whether CP is for ARL or WAD simulation
  if (v >= max_sample) {
    is_ARL <- TRUE
  } else {
    is_ARL <- FALSE
  }


  # Compute optimal delta star
  delta_star <- (m_post - m_pre) / sig_post^2

  # delta range = delta_star * / range_factor
  delta_upper <- delta_star * range_factor
  delta_lower <- delta_star / range_factor

  # Other settings
  psi_fn_list <- generate_sub_G_fn(sig_pre)
  v_min <- 1
  k_max <- 1e+3

  # Build CP detectors
  # When delta_lower = delta_upper = delta_star
  stcp_star <- build_stcp_exp(
    alpha,
    m_bound,
    delta_star,
    delta_star,
    is_test = FALSE,
    psi_fn_list,
    s_fn = function(x) {
      x - m_bound
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
    m_bound,
    delta_lower,
    delta_upper,
    is_test = FALSE,
    psi_fn_list,
    s_fn = function(x) {
      x - m_bound
    },
    v_fn = function(x) {
      1
    },
    v_min,
    k_max,
    tol = 1e-6
  )

  compute_single_run_stat <- function() {
    # Generate sample
    x_vec <- generator(
      max_sample,
      v,
      pre_sampler = function(n) {
        rnorm(n, m_pre, sig_pre)
      },
      post_sampler = function(n) {
        rnorm(n, m_post, sig_post)
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
         sig_pre = sig_pre,
         m_post = m_post,
         sig_post = sig_post
    )
  )
}


# ARL study

m_bound_vec <- seq(-0.4, 0.4, by = 0.2)
max_sample <- 200L

ARL_fn <- function(m_bound){
  out <- gaussian_simulator(m_pre = 0 ,
                            sig_pre = 1,
                            m_post = 0.5,
                            sig_post = 1,
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

m_post_vec <- seq(0.1, 1, by = 0.2)
max_sample <- 200L

WAD_fn <- function(m_post){
  out <- gaussian_simulator(m_pre = 0 ,
                            sig_pre = 1,
                            m_post = m_post,
                            sig_post = 1,
                            max_sample = max_sample,
                            v = max_sample / 10, # WAD study
                            m_bound = 0, # Correctly specified boundary
                            alpha = 1 / max_sample, # Set ARL >= max_sample
                            range_factor = 10,
                            num_repeat = 100L)

  return(c(
    m_post = m_post,
    WAD = out$run_avg,
    sd = out$run_sd
  ))
}


WAD_result <- lapply(m_post_vec, WAD_fn)
WAD_result <- do.call('rbind', WAD_result)
