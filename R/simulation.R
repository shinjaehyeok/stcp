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
#' @param stcp_mix  stcp object for the mixture case
#' @param v Change-point
#' @param stcp_star stcp object for the delta_star (optional)
#'
#' @return Stopped point for single and mixtures of SR- and CUSUM-type e-detectors.
#' @export
#'
run_quick_simulation <- function(x_vec,
                                 stcp_mix,
                                 v = Inf,
                                 stcp_star = NULL) {
  # Make test stcp obj
  stcp_test_mix <- stcp_mix
  stcp_test_mix$is_test <- TRUE
  if (!is.null(stcp_star)){
    stcp_test_star <- stcp_star
    stcp_test_star$is_test <- TRUE
  }

  if (!is.null(stcp_star)){
    # When delta_lower = delta_upper = delta_star

    # e-value for testing
    test_star_out <- run_stcp(x_vec, stcp_test_star)

    # e-detector 1. SR-type
    stcp_star_SR_out <- run_stcp(x_vec, stcp_star)

    # e-detector 2. CUSUM-type
    stcp_star_CS_out <- run_stcp(x_vec, stcp_star, is_SR_type = FALSE)

    plot(stcp_star_SR_out,
         xlab = "n",
         col = 2,
         ylim = c(
           min(test_star_out$log_mix_e_vec),
           max(stcp_star_SR_out$log_mix_e_vec)
         ))
    graphics::abline(v = v, col = 1, lty = 3)
    plot(test_star_out, add = TRUE, col = 1)
    plot(stcp_star_CS_out, add = TRUE, col = 3)
  }

  # When delta_lower < delta_star < delta_upper
  # e-value for testing
  test_mix_out <- run_stcp(x_vec,
                           stcp_test_mix)

  # e-detector 1. SR-type
  stcp_mix_SR_out <- run_stcp(x_vec, stcp_mix)

  # e-detector 2. CUSUM-type
  stcp_mix_CS_out <- run_stcp(x_vec, stcp_mix, is_SR_type = FALSE)

  plot(stcp_mix_SR_out,
       xlab = "n",
       ylim = c(
         min(test_mix_out$log_mix_e_vec),
         max(stcp_mix_SR_out$log_mix_e_vec)
       ),
       col = 2
  )
  graphics::abline(v = v, col = 1, lty = 3)
  plot(test_mix_out, add = TRUE, col = 1)
  plot(stcp_mix_CS_out, add = TRUE, col = 3)

  if (!is.null(stcp_star)){
    # Plot all
    plot(stcp_mix_SR_out,
         xlab = "n",
         ylim = c(
           min(test_star_out$log_mix_e_vec),
           max(stcp_mix_SR_out$log_mix_e_vec)
         ),
         col = 2
    )
    graphics::abline(v = v, col = 1, lty = 3)
    plot(test_mix_out, add = TRUE, col = 1)
    plot(stcp_mix_CS_out, add = TRUE, col = 3)

    plot(test_star_out, add = TRUE, col = 1, lty = 2)
    plot(stcp_star_SR_out, add = TRUE, col = 2, lty = 2)
    plot(stcp_star_CS_out, add = TRUE, col = 3, lty = 2)
  }

  out <-   list(
    mix_SR_stop = stcp_mix_SR_out$stopped_ind,
    mix_CS_stop = stcp_mix_CS_out$stopped_ind
  )

  if (!is.null(stcp_star)){
    out$star_SR_stop <- stcp_star_SR_out$stopped_ind
    out$star_CS_stop <- stcp_star_CS_out$stopped_ind
  }
  return(out)
}
