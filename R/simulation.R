#' Generate a vector of simulated observations
#'
#' Generate a simulated observations for given pre- and post-observation samplers.
#'
#' @param max_sample Total number of simulated observations
#' @param v change-point
#' @param pre_sampler R function to generate pre-change observations.
#' @param post_sampler R function to generate post-change observations.
#'
#' @return A vector of simulated observations
#' @export
#'
generator <- function(max_sample = 1000L,
                      v= 500L,
                      pre_sampler,
                      post_sampler){
  x_vec <- numeric(max_sample)
  if (v > 0){
    x_vec[1:v] <- pre_sampler(v)
  }
  if (v < max_sample){
    x_vec[seq(v+1, max_sample)] <- post_sampler(max_sample - v)
  }
  return(x_vec)
}

#' Run a quick simulation for debuging
#'
#' Run single and mixtures of e-val, e-detectors for given simulation setting.
#'
#' @param x_vec Simulated observations
#' @param v change-point
#' @param alpha Inverse of ARL
#' @param delta_star delta for sinlge e-value / e-detetor
#' @param delta_upper Upper bound of \code{delta_star}
#' @param delta_lower Lower bound of \code{delta_star}
#' @param psi_fn_list psi and psi star functions
#' @param v_min Minimum of variance process
#' @param k_max Maximum of non-trivial baselines.
#'
#' @return Stopped point for single and mixtures of SR- and CUSUM-type e-detectors.
#' @export
#'
run_quick_simulation <- function(
  x_vec,
  v,
  alpha,
  delta_star,
  delta_upper,
  delta_lower,
  psi_fn_list,
  v_min,
  k_max
){
  max_sample <- length(x_vec)
  # When delta_lower = delta_upper = delta_star
  base_param <- compute_baseline(
    alpha,
    delta_lower = delta_star,
    delta_upper = delta_star,
    psi_fn_list = psi_fn_list,
    v_min = v_min,
    k_max = k_max
  )
  log_one_over_alpha <- log(1/alpha)

  # Compute e-detectors
  log_base_fn_list <- sapply(base_param$lambda,
                             generate_log_base_fn,
                             psi_fn = base_param$psi_fn_list$psi)

  # e-value for testing (it is not the scope of this package but good for debugging)
  # log_e_val_vec <- cumsum(log_base_fn_list[[1]](x_vec))
  single_e_val <- update_log_mix_e_values(x_vec,
                                          base_param$omega,
                                          log_base_fn_list)


  plot(1:max_sample, single_e_val$log_mix_e_val, type = "l",
       xlab = "n", ylab = "e-value")
  graphics::abline(h = log_one_over_alpha, col = 2)
  graphics::abline(v = v, col = 1, lty = 2)

  # e-detector 1. SR-type
  single_SR <- update_log_mix_e_detectors(x_vec,
                                          base_param$omega,
                                          log_base_fn_list)

  # e-detector 2. CUSUM-type
  single_CS <- update_log_mix_e_detectors(x_vec,
                                          base_param$omega,
                                          log_base_fn_list,
                                          is_SR_type = FALSE)


  plot(1:max_sample, single_e_val$log_mix_e_val, type = "l",
       xlim = c(0, max_sample),
       ylim = c(min(single_e_val$log_mix_e_val), max(single_SR$log_mix_e_detect_val)),
       xlab = "n", ylab = "log e-val", main = paste0("v = ", v))
  graphics::abline(h = log_one_over_alpha, col = 2)
  graphics::abline(v = v, col = 1, lty = 2)
  graphics::lines(1:max_sample, single_SR$log_mix_e_detect_val, type = "l", col = 3)
  graphics::lines(1:max_sample, single_CS$log_mix_e_detect_val, type = "l", col = 4)

  # When delta_lower < delta_star < delta_upper
  base_param <- compute_baseline(
    alpha,
    delta_lower = delta_lower,
    delta_upper = delta_upper,
    psi_fn_list = psi_fn_list,
    v_min = v_min,
    k_max = k_max
  )

  # Compute e-detectors
  log_base_fn_list <- sapply(base_param$lambda,
                             generate_log_base_fn,
                             psi_fn = base_param$psi_fn_list$psi)

  # e-value for testing (it is not the scope of this package but good for debugging)

  mix_e_val <- update_log_mix_e_values(x_vec,
                                       base_param$omega,
                                       log_base_fn_list)

  plot(1:max_sample, mix_e_val$log_mix_e_val, type = "l",
       xlab = "n", ylab = "e-value")
  graphics::abline(h = log_one_over_alpha, col = 2)
  graphics::abline(v = v, col = 1, lty = 2)

  # e-detector 1. SR-type
  mix_SR <- update_log_mix_e_detectors(x_vec,
                                       base_param$omega,
                                       log_base_fn_list)


  # e-detector 2. CUSUM-type
  mix_CS <- update_log_mix_e_detectors(x_vec,
                                       base_param$omega,
                                       log_base_fn_list,
                                       is_SR_type = FALSE)



  single_SR_stop <- min(
    which(single_SR$log_mix_e_detect_val > log_one_over_alpha)
  )
  single_CS_stop <- min(
    which(single_CS$log_mix_e_detect_val > log_one_over_alpha)
  )
  mix_SR_stop <- min(
    which(mix_SR$log_mix_e_detect_val > log_one_over_alpha)
  )
  mix_CS_stop <- min(
    which(mix_CS$log_mix_e_detect_val > log_one_over_alpha)
  )

  # Plot all
  plot(1:max_sample, mix_e_val$log_mix_e_val, type = "l",
       xlab = "n", ylab = "log e-val", main = paste0("v = ", v),
       xlim = c(0, max_sample),
       ylim = c(min(single_e_val$log_mix_e_val),
                max(mix_SR$log_mix_e_detect_val)))
  graphics::lines(1:max_sample, single_e_val$log_mix_e_val, type = "l", lty = 2)
  graphics::lines(1:max_sample, mix_SR$log_mix_e_detect_val, type = "l", col = 3)
  graphics::lines(1:max_sample, single_SR$log_mix_e_detect_val, type = "l", col = 3, lty =2)
#  graphics::lines(1:max_sample, mix_CS$log_mix_e_detect_val, type = "l", col = 4)
#  graphics::lines(1:max_sample, single_CS$log_mix_e_detect_val, type = "l", col = 4, lty = 2)

  graphics::abline(h = log_one_over_alpha, col = 2)
  graphics::abline(h = log_one_over_alpha, col = 2, lty = 2)
  graphics::abline(v = v, col = 1, lty = 2)
  graphics::abline(v = mix_SR_stop, col = 3)
  graphics::abline(v = single_SR_stop, col = 3, lty = 2)
#  graphics::abline(v = mix_CS_stop, col = 4)
#  graphics::abline(v = single_CS_stop, col = 4, lty = 2)

  return(
    list(
      single_SR_stop = single_SR_stop,
      single_CS_stop = single_CS_stop,
      mix_SR_stop = mix_SR_stop,
      mix_CS_stop = mix_CS_stop
      )
  )
}
