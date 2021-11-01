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
