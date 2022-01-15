#' Compute lower bounds of confidence sequence for simple exponential e-values.
#'
#' Given a vector of observations, compute lower bounds of confidence sequence.
#'
#' @param x_vec A vector of observations
#' @param edcp_ci_obj EDCP_CI object, output of \code{build_ci_exp}
#' @param x_bar_pre Sample mean of previous observations. Default is \code{0} (No pre-samples)
#' @param n_pre Number of pre-samples. Default is \code{0}.
#' @inheritParams compute_single_ci
#'
#' @return List of 1. running sample mean, 2. sample size vector, 3. lower bounds of confidence sequence.
#' @export
#'
compute_ci <- function(x_vec,
                       edcp_ci_obj,
                       width_upper = 100,
                       ci_lower_trivial = -Inf,
                       x_bar_pre = 0,
                       n_pre = 0,
                       tol = 1e-6) {
  if (class(edcp_ci_obj) != "EDCP_CI") {
    stop("edcp_ci_obj must be an output of build_ci_exp function.")
  }

  ci_helper <- edcp_ci_obj$ci_helper

  if (n_pre < 0) {
    stop("Pre-sample size must be nonnegative.")
  }
  n_vec <- seq_along(x_vec) + n_pre
  x_sum_vec <- cumsum(x_vec) + x_bar_pre * n_pre
  x_bar_vec <- x_sum_vec / n_vec

  ci_lower_vec <-
    sapply(seq_along(n_vec), function(i) {
      compute_single_ci(n_vec[i],
                        x_bar_vec[i],
                        ci_helper,
                        width_upper,
                        ci_lower_trivial,
                        tol)
    })

  return(list(
    x_bar = x_bar_vec,
    n = n_vec,
    ci_lower = ci_lower_vec
  ))
}

#' Compute lower bounds of confidence sequence for simple exponential e-values.
#'
#' Given a vector of observations, compute lower bounds of confidence sequence.
#'
#' @param x_vec A vector of observations
#' @param edcp_bf_ci_obj EDCP_BF_CI object, output of \code{build_bf_ci_bounded}
#' @param max_num_ci The maximum number of CI to computes. Recomend to set it be less than 100.
#' @inheritParams compute_bf_ci_single
#'
#' @return List of 1. running sample mean, 2. sample size vector, 3. lower bounds of confidence sequence.
#' @export
#'
compute_bf_ci <- function(x_vec,
                          edcp_bf_ci_obj,
                          width_upper = 100,
                          ci_lower_trivial = -Inf,
                          max_num_ci = 100,
                          tol = 1e-6) {
  if (class(edcp_bf_ci_obj) != "EDCP_BF_CI") {
    stop("edcp_bf_ci_obj must be an output of the build_bf_ci_bounded function.")
  }
  if (max_num_ci < 1) {
    stop("max_num_ci must be at least 1")
  }

  bf_ci_helper <- edcp_bf_ci_obj$bf_ci_helper

  if (length(x_vec) <= max_num_ci) {
    n_ci_vec <- seq_along(x_vec)
  } else if (max_num_ci < 2) {
    n_vec_vec <- length(x_vec)
  } else {
    n_ci_vec <- unique(round(seq(1, length(x_vec) - 1, length.out = max_num_ci - 1)))
    n_ci_vec <- c(n_ci_vec, length(x_vec))
  }

  ci_lower_vec <-
    sapply(n_ci_vec, function(n) {
      compute_bf_ci_single(x_vec[1:n],
                           bf_ci_helper,
                           width_upper = width_upper,
                           ci_lower_trivial = ci_lower_trivial,
                           tol = tol)
    })

  x_bar_vec <- cumsum(x_vec) / seq_along(x_vec)


  return(list(
    x_bar = x_bar_vec[n_ci_vec],
    n = n_ci_vec,
    ci_lower = ci_lower_vec
  ))
}
