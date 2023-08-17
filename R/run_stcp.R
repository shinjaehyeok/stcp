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
  # TODO: This copy is very expensive.
  # We need to make the run_stcp as a method of R6 stcp_obj object.
  mix_e_list$stcp_obj <- stcp_obj
  class(mix_e_list) <- c("stcp_run", class(stcp_obj))
  return(mix_e_list)
}


#' Combine two stcp runs for a given weight.
#'
#' Combine two stcp runs for a given weight. Only \code{stcp_run} objects of the same \code{is_test} and same number of new observations can be combined.
#'
#' @param run1 First \code{stcp_run} object to combine.
#' @param run2 Second \code{stcp_run} object to combine.
#' @param w Weight scalar to determine the importance of the first \code{stcp_run} object. Default is 0.5.
#'
#' @return Combined \code{stcp_run} object of 1. logarithm of mixture of e-detectors, 2. the first stopped point and 3. Underlying stcp objects.
#' @export
#'
combine_stcp_run <- function(run1, run2, w = 0.5) {
  if (w <= 0 || w >= 1 ) stop("w must be a number between (0,1)")
  if (run1$stcp_obj$is_test != run2$stcp_obj$is_test) {
    stop("e-value (is_test == TRUE) and e-detector (is_test != FALSE) cannot be combined")
  }
  if (length(run1$log_mix_e_vec) != length(run2$log_mix_e_vec)) {
    stop("Run lengths of two stcp runs must be equal to each other.")
  }

  # Combine e-values or e-detectors
  combined_e_vec <- log(w * exp(run1$log_mix_e_vec) + (1-w) * exp(run2$log_mix_e_vec))

  # Combine stcp objects
  combined_stcp <- combine_stcp(run1$stcp_obj, run2$stcp_obj, w)

  # Compute new stopped index
  upcrossed_indice <-
    which(combined_e_vec > combined_stcp$log_one_over_alpha)
  new_stopped_ind <-
    ifelse(length(upcrossed_indice) == 0, Inf, min(upcrossed_indice))

  out <- list(
    log_mix_e_vec =  combined_e_vec,
    stopped_ind = new_stopped_ind,
    stcp_obj = combined_stcp
  )
  class(out) <-c("stcp_run", class(combined_stcp))
  return(out)

}
