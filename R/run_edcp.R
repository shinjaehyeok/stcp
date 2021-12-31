#' E-Detector based Change-Point detection and sequential tests.
#'
#' Run mixtures of e-detectors / e-values to detect the change-point or run a sequential test.
#'
#' @param new_x A list or vector of observations
#' @param edcp_obj EDCP object. It must contains model parameters (omega, log_base_fn_list) and previous log e-detector / values.
#' @param is_SR_type A Boolean to indicator whether the process use SR-type update. If not, the function uses the CUSUM-type update instead of SR-type one.
#'
#' @return logarithm of mixture of e-detectors, the first stopped point and updated edcp object.
#' @export
run_edcp <- function(new_x,
                     edcp_obj,
                     is_SR_type = TRUE) {
  if (length(new_x) == 0) {
    stop("new_x must have non-zero length.")
  }

  if (!edcp_obj$is_test) {
    # CP detection
    mix_e_list <- update_log_mix_e_detectors(
      new_x,
      edcp_obj$omega,
      edcp_obj$log_base_fn_list,
      edcp_obj$log_e_vec,
      is_SR_type
    )

  } else{
    # Sequential test
    mix_e_list <- update_log_mix_e_values(new_x,
                                          edcp_obj$omega,
                                          edcp_obj$log_base_fn_list,
                                          edcp_obj$log_e_vec)
  }

  # Compute stopped index
  upcrossed_indice <-
    which(mix_e_list$log_mix_e_vec > edcp_obj$log_one_over_alpha)
  mix_e_list$stopped_ind <-
    ifelse(length(upcrossed_indice) == 0, Inf, min(upcrossed_indice))

  # Update log_e_vec
  edcp_obj$log_e_vec <- mix_e_list$last_log_e_vec
  mix_e_list$last_log_e_vec <- NULL
  mix_e_list$edcp_obj <- edcp_obj
  class(mix_e_list) <- c("edcp_out", class(edcp_obj))
  return(mix_e_list)
}
