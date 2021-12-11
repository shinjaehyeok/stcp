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
                      v = 500L,
                      pre_sampler,
                      post_sampler) {
  x_vec <- numeric(max_sample)
  if (v > 0) {
    x_vec[1:v] <- pre_sampler(v)
  }
  if (v < max_sample) {
    x_vec[seq(v + 1, max_sample)] <- post_sampler(max_sample - v)
  }
  return(x_vec)
}

#' Run a quick simulation for debuging
#'
#' Run single and mixtures of e-val, e-detectors for given simulation setting.
#'
#' @param x_vec Simulated observations
#' @param v change-point
#' @param baseline_star baseline object for the delta_star
#' @param baseline_mix  baseline object for the mixture case
#'
#' @return Stopped point for single and mixtures of SR- and CUSUM-type e-detectors.
#' @export
#'
run_quick_simulation <- function(x_vec,
                                 v,
                                 baseline_star,
                                 baseline_mix) {
  max_sample <- length(x_vec)

  # When delta_lower = delta_upper = delta_star
  # e-value for testing
  single_e_val <- edcp(x_vec,
                       baseline_star,
                       is_test = TRUE)

  plot(
    1:max_sample,
    single_e_val$log_mix_e_val,
    type = "l",
    xlab = "n",
    ylab = "e-value"
  )
  graphics::abline(h = single_e_val$threshold, col = 2)
  graphics::abline(v = v, col = 1, lty = 2)

  # e-detector 1. SR-type
  single_SR <- edcp(x_vec, baseline_star)

  # e-detector 2. CUSUM-type
  single_CS <- edcp(x_vec, baseline_star, is_SR_type = FALSE)


  plot(
    1:max_sample,
    single_e_val$log_mix_e_val,
    type = "l",
    xlim = c(0, max_sample),
    ylim = c(
      min(single_e_val$log_mix_e_val),
      max(single_SR$log_mix_e_val)
    ),
    xlab = "n",
    ylab = "log e-val",
    main = paste0("v = ", v)
  )
  graphics::abline(h = single_e_val$threshold, col = 2)
  graphics::abline(v = v, col = 1, lty = 2)
  graphics::lines(1:max_sample,
                  single_SR$log_mix_e_val,
                  type = "l",
                  col = 3)
  graphics::lines(1:max_sample,
                  single_CS$log_mix_e_val,
                  type = "l",
                  col = 4)


  # When delta_lower < delta_star < delta_upper
  # e-value for testing
  mix_e_val <- edcp(x_vec,
                    baseline_mix,
                    is_test = TRUE)

  plot(
    1:max_sample,
    mix_e_val$log_mix_e_val,
    type = "l",
    xlab = "n",
    ylab = "e-value"
  )
  graphics::abline(h = mix_e_val$threshold, col = 2)
  graphics::abline(v = v, col = 1, lty = 2)

  # e-detector 1. SR-type
  mix_SR <- edcp(x_vec, baseline_mix)

  # e-detector 2. CUSUM-type
  mix_CS <- edcp(x_vec, baseline_mix, is_SR_type = FALSE)

  # Plot all
  plot(
    1:max_sample,
    mix_e_val$log_mix_e_val,
    type = "l",
    xlab = "n",
    ylab = "log e-val",
    main = paste0("v = ", v),
    xlim = c(0, max_sample),
    ylim = c(
      min(single_e_val$log_mix_e_val),
      max(mix_SR$log_mix_e_val)
    )
  )
  graphics::lines(1:max_sample,
                  single_e_val$log_mix_e_val,
                  type = "l",
                  lty = 2)
  graphics::lines(1:max_sample,
                  mix_SR$log_mix_e_val,
                  type = "l",
                  col = 3)
  graphics::lines(
    1:max_sample,
    single_SR$log_mix_e_val,
    type = "l",
    col = 3,
    lty = 2
  )
  #  graphics::lines(1:max_sample, mix_CS$log_mix_e_val, type = "l", col = 4)
  #  graphics::lines(1:max_sample, single_CS$log_mix_e_val, type = "l", col = 4, lty = 2)

  graphics::abline(h = mix_e_val$threshold, col = 2)
  graphics::abline(h = mix_e_val$threshold, col = 2, lty = 2)
  graphics::abline(v = v, col = 1, lty = 2)
  graphics::abline(v = mix_SR$stopped_ind, col = 3)
  graphics::abline(v = single_SR$stopped_ind, col = 3, lty = 2)
  #  graphics::abline(v = mix_CS_stop, col = 4)
  #  graphics::abline(v = single_CS_stop, col = 4, lty = 2)

  return(
    list(
      single_SR_stop = single_SR$stopped_ind,
      single_CS_stop = single_CS$stopped_ind,
      mix_SR_stop = mix_SR$stopped_ind,
      mix_CS_stop = mix_CS$stopped_ind
    )
  )
}
