#' Compute the width of confidence interval for simple exponential e-values.
#'
#' Given the number of samples n, compute the width of confidence interval for simple exponential e-values.
#'
#' @param n Number of samples
#' @param ci_helper R function from \code{compute_baseline_simple_exp} used to compute the width.
#' @param ci_upper_const Constant to compute upper bound of the width of confidence interval of the form \code{ci_upper_const / sqrt(n)}.
#' @param tol Tolerance of root-finding, positive numeric. Default is 1e-6.
#'
#' @return Width of confidence interval given the sample size.
#' @export
#'
compute_single_ci_width <- function(n,
                                    ci_helper,
                                    ci_upper_const = 100,
                                    tol = 1e-6) {
  #TODO::IMPLETMENT check whether ci_helper is simple_exp class.

  f <- function(d) {
    ci_helper(d, n)
  }

  ci_upper <- ci_upper_const / sqrt(n)

  if (f(ci_upper) <= 0) {
    stop("ci_upper is too small to compute CI.")
  }

  ci_width <- stats::uniroot(f,
                             c(0, ci_upper), tol = tol)

  return(ci_width$root)
}

#' Compute confidence sequence for simple exponential e-values.
#'
#' Given a vector of observations, compute widths of confidence sequence.
#'
#' @param x_vec A vector of observations
#' @param x_bar_pre Sample mean of previous observations. Default is \code{0} (No pre-samples)
#' @param n_pre Number of pre-samples. Default is \code{0}.
#' @inheritParams compute_single_ci_width
#'
#' @return List of 1. running sample mean, 2. sample size vector, 3. widths of confidence sequence.
#' @export
#'
compute_ci_width <- function(x_vec,
                             ci_helper,
                             ci_upper_const = 100,
                             x_bar_pre = 0,
                             n_pre = 0,
                             tol = 1e-6) {
  #TODO::IMPLETMENT check whether ci_helper is simple_exp class.
  if (n_pre < 0) {
    stop("Pre-sample size must be nonnegative.")
  }
  n_vec <- seq_along(x_vec) + n_pre
  x_sum_vec <- cumsum(x_vec) + x_bar_pre * n_pre
  x_bar_vec <- x_sum_vec / n_vec

  ci_width_vec <-
    sapply(n_vec, function(n) {
      compute_single_ci_width (n,
                               ci_helper,
                               ci_upper_const,
                               tol)
    })

  return(list(
    x_bar = x_bar_vec,
    n = n_vec,
    ci_width = ci_width_vec
  ))
}

#' Convert target sample interval to delta bounds
#'
#' Given a vector of observations, compute widths of confidence sequence.
#'
#' @param n_upper Upper bound of the target sample interval
#' @param n_lower Lower bound of the target sample interval. Default is \code{1}.
#' @inheritParams compute_baseline
#'
#' @return \code{delta_lower}, \code{delta_upper} in compute_baseline for CI computation.
#' @export
#'
convert_time_to_delta_bound <- function(alpha,
                                        n_upper,
                                        n_lower,
                                        psi_fn_list = generate_sub_G_fn(),
                                        v_min = 1,
                                        k_max = 200,
                                        tol = 1e-6) {

  if (!(n_lower > 0 & n_upper >=  n_lower)) {
    stop("n_lower and n_upper must be positive with n_lower <= n_upper.")
  }

  if (n_lower == n_upper){
    g_alpha <- log(1/alpha)
  } else {
    delta_init <-
      psi_fn_list$psi_star_inv(log(1 / alpha)  / c(n_lower, n_upper))
    baseline_init <- compute_baseline(alpha,
                                      delta_init[2],
                                      delta_init[1],
                                      psi_fn_list,
                                      v_min,
                                      k_max)
    g_alpha <- baseline_init$g_alpha
  }

  # Compute delta bound
  delta_lower <-
    psi_fn_list$psi_star_inv(g_alpha / n_upper)
  delta_upper <-
    psi_fn_list$psi_star_inv(g_alpha  / n_lower)

  return(list(delta_lower = delta_lower, delta_upper = delta_upper))
}


#' Brute-force computation of confidence interval.
#'
#' Given a vector of observations, compute the lower bound of confidence interval at the last index.
#'
#' @param x_vec A vector of observations
#' @param bruth_force_ci_helper R function used to compute the lower bound of confidence interval.
#' @param ci_upper_const Constant to compute upper bound of the width of confidence interval of the form \code{ci_upper_const / sqrt(n)}.
#' @param tol Tolerance of root-finding, positive numeric. Default is 1e-6.
#'
#' @return Lower bound of confidence interval given the sample size.
#' @export
#'
brute_force_lower_ci <- function(x_vec,
                                 bruth_force_ci_helper,
                                 ci_upper_const = 100,
                                 tol = 1e-6) {
  f <- function(m) {
    bruth_force_ci_helper(m, x_vec)
  }
  n <- length(x_vec)
  x_bar <- mean(x_vec)
  ci_lower <- x_bar - ci_upper_const / sqrt(n)

  if (f(ci_lower) <= 0) {
    stop("ci_upper_const is too small to compute CI.")
  }

  ci_width <- stats::uniroot(f,
                             c(ci_lower, x_bar), tol = tol)

  return(ci_width$root)
}
