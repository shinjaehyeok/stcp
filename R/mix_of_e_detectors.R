#' Log of mixtures of SR- or CUSUM-type e-detectors
#'
#' Compute logarithm of mixtures of SR- or CUSUM-type e-detector based on the current observation and previous value.
#'
#' @param new_x_list A list of new observations. If each observation is univariate then it can be a vector instead of list of single observations.
#' @param weight_vec A vector of mixing weights.
#' @param log_base_fn_list R function that compute the log of baseline process based on each observation.
#' @param prev_log_e_vec A vector of logarithms of mixtures of e-detectors in the previous step. Default is \code{numeric(length(log_base_fn_list))}.
#' @param is_SR_type A Boolean to indicator whether the process use SR-type update. If not, the function uses the CUSUM-type update instead of SR-type one.
#'
#' @return Updated logarithm of mixture of SR- or CUSUM-type e-detectors and a vector of each component for the next iteration.
#' @export
#'
update_log_mix_e_detectors <- function(new_x_list,
                                       weight_vec,
                                       log_base_fn_list,
                                       prev_log_e_vec = numeric(length(log_base_fn_list)),
                                       is_SR_type = TRUE){
  # Check weights are positive and summed to one
  if (any(weight_vec <= 0)) stop("Weights must be positive")
  w <- sum(weight_vec)
  if (abs(w-1) > 1e-6) stop("Sum of weights must be equal to 1")

  # Compute the logs of weights
  log_weight_vec <- log(weight_vec)

  # Check the e_val, weights vectors and function list are consistent
  num_base <- length(log_base_fn_list)
  if (length(log_weight_vec) != num_base){
    stop("Lengths of log_weight_vec and log_base_fn_list are not matched.")
  }
  if (length(prev_log_e_vec) != num_base){
    stop("Lengths of prev_log_e_vec and log_base_fn_list are not matched.")
  }


  num_new_x <- length(new_x_list)
  current_ind <- 1
  log_mix_e_detect_val <- numeric(num_new_x)

  updater <- function(){
    update_log_e_for_mix_ind <- function(mix_ind){
      update_log_e_detector(
        x_current = new_x_list[[current_ind]],
        prev_log_e = prev_log_e_vec[mix_ind],
        compute_log_baseline = log_base_fn_list[[mix_ind]],
        is_SR_type = is_SR_type
      )
    }
    prev_log_e_vec <<- sapply(1:num_base, update_log_e_for_mix_ind)
    log_mix_e_detect_val[current_ind] <<- matrixStats::logSumExp(prev_log_e_vec + log_weight_vec)
    current_ind <<- current_ind + 1
  }
  while(current_ind <= num_new_x){
    updater()
  }

  updated_list <- list(log_mix_e_vec = log_mix_e_detect_val,
                       last_log_e_vec = prev_log_e_vec)

  return(updated_list)
}
