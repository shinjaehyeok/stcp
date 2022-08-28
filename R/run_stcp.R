#' E-Detector based Change-Point detection and sequential tests.
#'
#' Run mixtures of e-detectors / e-values to detect the change-point or run a sequential test.
#'
#' @param new_x A list or vector of observations
#' @param stcp_obj STCP object. It must contains model parameters (omega, log_base_fn_list) and previous log e-detector / values.
#' @param is_SR_type A Boolean to indicator whether the process use SR-type update. If not, the function uses the CUSUM-type update instead of SR-type one.
#'
#' @return logarithm of mixture of e-detectors, the first stopped point and updated stcp object.
#' @export
run_stcp <- function(new_x,
                     stcp_obj,
                     is_SR_type = TRUE) {
  if (length(new_x) == 0) {
    stop("new_x must have non-zero length.")
  }

  if (is.matrix(new_x) | is.data.frame(new_x)) {
    stop("new_x must be either vector or list of observations.")
  }

  # Faster vectorized-computation for the univariate sequential test case
  if (!is.list(new_x) & stcp_obj$is_test) {
    mix_e_list <- update_log_mix_e_values(new_x,
                                          stcp_obj$omega,
                                          stcp_obj$log_base_fn_list,
                                          stcp_obj$log_e_vec)
  } else {
    mix_e_list <- update_log_mix_e(
      new_x,
      stcp_obj$omega,
      stcp_obj$log_base_fn_list,
      stcp_obj$log_e_vec,
      stcp_obj$is_test,
      is_SR_type
    )

  }

  # Compute stopped index
  upcrossed_indice <-
    which(mix_e_list$log_mix_e_vec > stcp_obj$log_one_over_alpha)
  mix_e_list$stopped_ind <-
    ifelse(length(upcrossed_indice) == 0, Inf, min(upcrossed_indice))

  # Update log_e_vec
  stcp_obj$log_e_vec <- mix_e_list$last_log_e_vec
  stcp_obj$n <- stcp_obj$n + length(new_x)
  mix_e_list$last_log_e_vec <- NULL
  mix_e_list$stcp_obj <- stcp_obj
  class(mix_e_list) <- c("stcp_run", class(stcp_obj))
  return(mix_e_list)
}
