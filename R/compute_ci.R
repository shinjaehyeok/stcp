#' Compute baseline parameters given target sample size bounds.
#'
#' Given target sample size bounds interval for confidence sequences, compute  baseline parameters.
#'
#' @param n_upper Upper bound of the target sample interval
#' @param n_lower Lower bound of the target sample interval. Default is \code{1}.
#' @inheritParams compute_baseline
#'
#' @return A list of baseline parameters to build \code{ci_helper}
#' @export
#'
compute_baseline_for_sample_size <- function(alpha,
                                             n_upper,
                                             n_lower,
                                             psi_fn_list = generate_sub_G_fn(),
                                             v_min = 1,
                                             k_max = 200,
                                             tol = 1e-6) {
  if (!(n_lower > 0 & n_upper >=  n_lower)) {
    stop("n_lower and n_upper must be positive with n_lower <= n_upper.")
  }

  if (n_lower == n_upper) {
    g_alpha <- log(1 / alpha)
  } else {
    delta_init_upper <-
      psi_fn_list$psi_star_inv(log(1 / alpha)  / n_lower)
    delta_init_lower <-
      psi_fn_list$psi_star_inv(log(1 / alpha)  / n_upper)
    baseline_init <- compute_baseline(alpha,
                                      delta_init_lower,
                                      delta_init_upper,
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

  # Compute baseline parameters
  baseline_param <- compute_baseline(alpha,
                                     delta_lower,
                                     delta_upper,
                                     psi_fn_list,
                                     v_min,
                                     k_max)
  return(baseline_param)
}



#' Generate ci_helper function for CI computation.
#'
#' Given baseline parameters, generate ci_helper function for CI computation.
#'
#' @param alpha Confidence level in (0,1)
#' @param omega A vector of omega parameters from \code{compute_baseline}, which is a mixing weights.
#' @param lambda A vector of lambda parameters from \code{compute_baseline}, which is a set of model parameters.
#' @param psi_fn_list_generator Function operator generating \code{psi_fn_list}
#' @param is_psi_depend_on_m Indicator whether the psi function depends on the target parameter.
#'
#' @return \code{ci_helper} function to compute confidence intervals
#' @export
#'
generate_ci_helper <- function(alpha,
                               omega,
                               lambda,
                               psi_fn_list_generator,
                               is_psi_depend_on_m = FALSE) {
  # Make helper function for CI computation
  log_weight_vec <- log(omega)
  threshold <- log(1 / alpha)

  if (is_psi_depend_on_m) {
    ci_helper <- function(m, n, x_bar) {
      psi_fn_list <- psi_fn_list_generator(m)
      psi_lambda_vec <- sapply(lambda, psi_fn_list$psi)
      inner <- n * (lambda * (x_bar - m) - psi_lambda_vec)
      out <-
        matrixStats::logSumExp(inner + log_weight_vec) - threshold
      return(out)
    }
  } else {
    psi_fn_list <- psi_fn_list_generator()
    psi_lambda_vec <- sapply(lambda, psi_fn_list$psi)
    ci_helper <- function(m, n, x_bar) {
      inner <- n * (lambda * (x_bar - m) - psi_lambda_vec)
      out <-
        matrixStats::logSumExp(inner + log_weight_vec) - threshold
      return(out)
    }
  }
  return(ci_helper)
}

#' Compute the lower bound of confidence interval for simple exponential e-values.
#'
#' Given the number of samples n, compute the lower bound of confidence interval for simple exponential e-values.
#'
#' @param n Number of samples
#' @param x_bar Sample mean of \code{n} observations
#' @param ci_helper R function from \code{compute_baseline_simple_exp} used to compute the width.
#' @param width_upper Constant to compute upper bound of the width of confidence interval of the form \code{width_upper / sqrt(n)}.
#' @param ci_lower_trivial A trivial lower bound of confidence interval. Default is \code{-Inf}.
#' @param tol Tolerance of root-finding, positive numeric. Default is 1e-6.
#'
#' @return A lower bound of confidence interval given the sample size.
#' @export
#'
compute_single_ci <- function(n,
                              x_bar,
                              ci_helper,
                              width_upper = 100,
                              ci_lower_trivial = -Inf,
                              tol = 1e-6) {
  #TODO::IMPLETMENT check whether ci_helper is simple_exp class.

  f <- function(m) {
    ci_helper(m, n, x_bar)
  }

  if (f(ci_lower_trivial) <= 0) {
    # We have no parameter to be able to reject.
    # Hence CI must be as large as possible.
    return(ci_lower_trivial)
  }

  ci_lower_bound <-
    max(c(x_bar - width_upper / sqrt(n), ci_lower_trivial))

  if (f(ci_lower_bound) <= 0) {
    stop("width_upper is too small to compute CI")
  }

  ci_lower <- stats::uniroot(f,
                             c(ci_lower_bound, x_bar), tol = tol)

  return(ci_lower$root)
}

#' Compute lower bounds of confidence sequence for simple exponential e-values.
#'
#' Given a vector of observations, compute lower bounds of confidence sequence.
#'
#' @param x_vec A vector of observations
#' @param edcp_ci_obj EDCP_CI object, output of \code{build_edcp_exp}
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
    stop("ci_helper must be an output of generate_ci_helper function.")
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


#' Brute-force computation of confidence interval.
#'
#' Given a vector of observations, compute the lower bound of confidence interval at the last index.
#'
#' @param bruth_force_ci_helper R function used to compute the lower bound of confidence interval.
#' @inheritParams compute_ci
#'
#' @return Lower bound of confidence interval given the sample size.
#' @export
#'
compute_brute_force_ci <- function(x_vec,
                                   bruth_force_ci_helper,
                                   width_upper = 100,
                                   ci_lower_trivial = -Inf,
                                   tol = 1e-6) {
  f <- function(m) {
    bruth_force_ci_helper(m, x_vec)
  }
  if (f(ci_lower_trivial) <= 0) {
    # We have no parameter to be able to reject.
    # Hence CI must be as large as possible.
    return(ci_lower_trivial)
  }
  n <- length(x_vec)
  x_bar <- mean(x_vec)
  ci_lower_bound <-
    max(c(x_bar - width_upper / sqrt(n), ci_lower_trivial))

  if (f(ci_lower_bound) <= 0) {
    stop("width_upper is too small to compute CI.")
  }
  ci_lower <- stats::uniroot(f,
                             c(ci_lower_bound, x_bar), tol = tol)

  return(ci_lower$root)
}
