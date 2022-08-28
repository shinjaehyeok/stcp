# If stcp is not installed then run below commands
# install.packages("devtools")
# devtools::install_github("shinjaehyeok/stcp")

library(stcp)

# Bernoulli case ----

ber_simulator <- function(p_pre = 0.3 ,
                          p_post = 0.5,
                          max_sample = 200L,
                          v = max_sample / 10,
                          p_bound = 0.3,
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
  delta_star <- (p_post - p_pre)

  # delta range = delta_star * / range_factor
  delta_upper <- delta_star * range_factor
  delta_upper <- min(c(delta_upper, 0.99 - p_pre)) # To make delta_upper + p_pre < 1
  delta_lower <- delta_star / range_factor

  # Other settings
  psi_fn_list <- generate_sub_B_fn(p = p_pre)
  v_min <- 1
  k_max <- 1e+3

  # Build CP detectors
  # When delta_lower = delta_upper = delta_star
  stcp_star <- build_stcp_exp(
    alpha,
    p_bound,
    delta_star,
    delta_star,
    is_test = FALSE,
    psi_fn_list,
    s_fn = function(x) {
      x - p_bound
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
    p_bound,
    delta_lower,
    delta_upper,
    is_test = FALSE,
    psi_fn_list,
    s_fn = function(x) {
      x - p_bound
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
        rbinom(n, 1, p_pre)
      },
      post_sampler = function(n) {
        rbinom(n, 1, p_post)
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
         p_bound = p_bound,
         max_sample = max_sample,
         v = v,
         p_pre = p_pre,
         p_post = p_post
    )
  )
}


# ARL study

p_pre <- 0.3
p_bound_vec <- p_pre + seq(-0.2, 0.2, by = 0.1)
max_sample <- 200L

ARL_fn <- function(p_bound){
  out <- ber_simulator(p_pre = p_pre ,
                       p_post = 0.5,
                       max_sample = max_sample,
                       v = max_sample, # ARL study
                       p_bound = p_bound,
                       alpha = 10 / max_sample,
                       range_factor = 10,
                       num_repeat = 100L)

  return(c(
    p_bound = p_bound,
    ARL_level = out$one_over_alpha,
    ARL = out$run_avg,
    sd = out$run_sd
  ))
}

ARL_result <- lapply(p_bound_vec, ARL_fn)
ARL_result <- do.call('rbind', ARL_result)

# Worst average delay study
p_pre <- 0.3
p_post_vec <- p_pre + seq(0.05, 0.25, by = 0.05)
max_sample <- 200L

WAD_fn <- function(p_post){
  out <- ber_simulator(p_pre = p_pre ,
                       p_post = p_post,
                       max_sample = max_sample,
                       v = max_sample / 10, # WAD study
                       p_bound = p_pre, # Correctly specified boundary
                       alpha = 1 / max_sample, # Set ARL >= max_sample
                       range_factor = 10,
                       num_repeat = 100L)

  return(c(
    p_pre = p_pre,
    p_post = p_post,
    WAD = out$run_avg,
    sd = out$run_sd
  ))
}


WAD_result <- lapply(p_post_vec, WAD_fn)
WAD_result <- do.call('rbind', WAD_result)
