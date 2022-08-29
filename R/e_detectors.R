#' Log of SR- or CUSUM-type e-detectors
#'
#' Compute logarithm of SR- or CUSUM-type e-detector based on the current observation and previous value.
#'
#' @param x_current The current observation
#' @param prev_log_e A numeric value of the mixture of e-detectors in the previous step. Default is \code{0}.
#' @param compute_log_baseline R function that compute the log of baseline process based on each observation.
#' @param is_test A Boolean to indicate whether this model is for a sequential test or not.
#' @param is_SR_type A Boolean indicator whether the process use SR-type update. If not, the function uses the CUSUM-type update instead of SR-type one.
#'
#' @return Updated logarithm of SR- or CUSUM-type e-detector.
#' @export
#'
update_log_e <- function(x_current,
                         prev_log_e = -Inf,
                         compute_log_baseline = function(x){x - 0.5},
                         is_test = FALSE,
                         is_SR_type = TRUE){

  current_log_e_val <- compute_log_baseline(x_current)
  if (is_test) {
    update <- prev_log_e
  } else {
    if (is_SR_type) {
      update <- matrixStats::logSumExp(c(prev_log_e, 0))
    } else{
      update <- max(prev_log_e, 0)
    }
  }
  return(current_log_e_val + update)
}

#' Generate log of baseline function
#'
#' Generate logarithm of the baseline function given lambda parameter, psi_star, sum and variance functions
#' @param lambda Lambda parameter of the target baseline function.
#' @param psi_fn R function that compute \eqn{\psi(\lambda)}.
#' @param s_fn R function that compute the sum process given an observation.
#' @param v_fn R function that compute the variance process given an observation.
#'
#' @return A function compute log of baseline process given an observation.
#' @export
#'
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
#' @param m Mean parameter of the target baseline function. It must be strictly larger than \code{bound_lower}.
#' @param bound_lower Lower bound of observations. Default is \code{0}.
#' @param bound_upper Upper bound of observations. Default is \code{1}.
#' @param is_flipped A Boolean to indicate whether the model should take a flipped input or not. If the input \eqn{X} is in \eqn{[0,1]} then the flipped input is defined by \eqn{1-X}.
#' @return A function compute log of baseline process given an observation.
#' @export
#'
generate_log_bounded_base_fn <- function(lambda,
                                         m = 0.5,
                                         bound_lower = 0,
                                         bound_upper = 1,
                                         is_flipped = FALSE){
  if (is_flipped) {
    if (m >= bound_upper){
      stop("Mean parameter must strictly smaller than bound_upper.")
    }
    m_gap <- bound_upper - m
    log_base_fn <- function(x){
      log(1 + lambda * ( (bound_upper - x) / m_gap  - 1))
    }
  } else {
    if (m <= bound_lower){
      stop("Mean parameter must strictly larger than bound_lower.")
    }
    m_gap <- m - bound_lower
    log_base_fn <- function(x){
      log(1 + lambda * ( (x - bound_lower) / m_gap  - 1))
    }
  }
  return(log_base_fn)
}
