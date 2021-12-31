#' Log of mixtures of e-values for sequential test
#'
#' Compute logarithm of mixtures of e-values based on the current observation and previous value.
#'
#' @param new_x_vec A vector of new observations.
#' @param weight_vec A vector of mixing weights.
#' @param log_base_fn_list R function that compute the log of baseline process based on each observation.
#' @param prev_log_e_vec A vector of logarithms of mixtures of e-values in the previous step. Default is \code{numeric(length(log_base_fn_list))}.
#'
#' @return Updated logarithm of mixture of e-values and a vector of each component for the next iteration.
#' @export
#'
update_log_mix_e_values <- function(new_x_vec,
                                    weight_vec,
                                    log_base_fn_list,
                                    prev_log_e_vec = numeric(length(log_base_fn_list))
                                    ){
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

  # Compute each L_i
  log_e_val_mat <- matrix(sapply(log_base_fn_list, function(f) f(new_x_vec)), ncol = length(log_base_fn_list))

  # Compute increment to prod L_i
  log_e_val_mat <- matrix(apply(log_e_val_mat, 2, cumsum), ncol = length(log_base_fn_list))

  # Add prev prod L_i
  log_e_val_mat <- sweep(log_e_val_mat, 2, prev_log_e_vec, "+")

  # Store last row for future update
  prev_log_e_vec <- log_e_val_mat[nrow(log_e_val_mat),]


  log_mix_e_val_vec <- apply(log_e_val_mat, 1,
                              function(log_e_vec){
                                matrixStats::logSumExp(log_e_vec + log_weight_vec)
                                })

  updated_list <- list(log_mix_e_vec = log_mix_e_val_vec,
                       last_log_e_vec = prev_log_e_vec)
  return(updated_list)
}
