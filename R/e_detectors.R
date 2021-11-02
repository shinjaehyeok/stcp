#' Log of SR- or CUSUM-type e-detectors
#'
#' Compute logarithm of SR- or CUSUM-type e-detector based on the current observation and previous value.
#'
#' @param x_current The current observation
#' @param prev_log_e A numeric value of the mixture of e-detectors in the previous step. Default is \code{0}.
#' @param compute_log_baseline R function that compute the log of baseline process based on each observation.
#' @param is_SR_type A Boolean indicator whether the process use SR-type update. If not, the function uses the CUSUM-type update instead of SR-type one.
#'
#' @return Updated logarithm of SR- or CUSUM-type e-detector.
#' @export
#'
#' @examples
#' update_log_e_detector(1, 0, is_SR_type = TRUE)
#' update_log_e_detector(1, 0, is_SR_type = FALSE)
update_log_e_detector <- function(x_current,
                                  prev_log_e = 0,
                                  compute_log_baseline = function(x){x - 0.5},
                                  is_SR_type = TRUE){

  current_log_e_val <- compute_log_baseline(x_current)
  if (is_SR_type){
     update <- matrixStats::logSumExp(c(prev_log_e, 0))
  } else{
     update <- max(prev_log_e, 0)
  }
  return(current_log_e_val + update)
}

#' Generate log of baseline function
#'
#' Generate logarithm of the baseline function given lambda parameter, psi_star, sum and variance functions
#' @param lambda Lambda parameter of the target baseline function.
#' @param psi_fn R function that compute \eqn{\psi(\lambda)}.
#' @param s_fn R function that compute the sum process given an observation.
#' @param v_fn R function that compute the sum process given an observation.
#'
#' @return A function compute log of baseline process given an observation.
#' @export
#'
#' @examples
#' generate_log_base_fn(1)
#' generate_log_base_fn(1, psi_fn = function(x) x^2)
generate_log_base_fn <- function(lambda,
                                 psi_fn = function(x){x^2/2},
                                 s_fn = function(x){x},
                                 v_fn = function(x){1}){

  log_base_fn <- function(x){
    lambda * s_fn(x) - psi_fn(lambda) * v_fn(x)
  }
  return(log_base_fn)
}

#' Generate log of baseline function for bounded RVs
#'
#' Generate logarithm of the baseline function for bounded RVs in \eqn{[0,1]} given lambda and mean parameters
#' @param lambda Lambda parameter of the target baseline function.
#' @param m Mean parameter of the target baseline function.
#' @return A function compute log of baseline process given an observation.
#' @export
#'
#' @examples
#' generate_log_bounded_base_fn(1)
#' generate_log_bounded_base_fn(0.5)
generate_log_bounded_base_fn <- function(lambda, m = 0.5){

  if (!(m > 0 | m < 1)) stop("Mean parameter must be in (0,1).")
  log_base_fn <- function(x){
    log(1 + lambda * (x / m  - 1))
  }
  return(log_base_fn)
}
