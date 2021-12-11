#' E-Detector based Change-Point detection and sequential tests.
#'
#' Run mixtures of e-detectors / e-values to detect the change-point or run a sequential test.
#'
#' @param x_vec Observations
#' @param baseline_obj Baseline object returned by compute_baseline functions.
#' @param is_test A Boolean to indicate whether to compute sequential test or not.
#' @inheritParams update_log_mix_e_detectors
#'
#' @return logarithm of mixture of e-detectors and the first stopped point.
#' @export
edcp <- function(x_vec,
                 baseline_obj,
                 prev_log_e_vec = numeric(length(baseline_obj$log_base_fn_list)),
                 is_SR_type = TRUE,
                 is_test = FALSE) {
  # Compute threshold
  log_one_over_alpha <- log(1 / baseline_star$alpha)

  if (!is_test) {
    # CP detection
    mix_e_list <- update_log_mix_e_detectors(
      x_vec,
      baseline_obj$omega,
      baseline_obj$log_base_fn_list,
      prev_log_e_vec,
      is_SR_type
    )

  } else{
    # Sequential test
    mix_e_list <- update_log_mix_e_values(x_vec,
                                          baseline_obj$omega,
                                          baseline_obj$log_base_fn_list)
  }

  # Update stopped index
  upcrossed_ind <-
    which(mix_e_list$log_mix_e_val > log_one_over_alpha)
  mix_e_list$stopped_ind <-
    ifelse(length(upcrossed_ind) == 0, Inf, min(upcrossed_ind))
  mix_e_list$threshold <- log_one_over_alpha

  return(mix_e_list)
}
