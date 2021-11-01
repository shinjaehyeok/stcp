#' Compute baseline processes
#'
#' Compuate parameters to build baseline processes.
#'
#' @param alpha ARL parameter in (0,1)
#' @param delta_lower Lower bound of target Delta. It must have be positive and smaller than or equal to \code{delta_upper}.
#' @param delta_upper Upper bound of target Delta. It must have be positive and larger than or equal to \code{delta_lower}.
#' @param psi_star R function that computes the convex conjugate of the psi function. Default is \code{function(x)x^2/2}.
#' @param psi_star_div R function that computes the derivative of \code{psi_star}. Default is \code{function(x)x}.
#' @param psi_star_inv R function that computes the inverse of \code{psi_star}. Default is \code{function(y)sqrt(2*y)}.
#' @param v_min A lower bound of v function in the baseline process. Default is \code{1}.
#' @param k_max Positive integer to determine the maximum number of baselines. Default is 100.
#' @param tol Tolerance of root-finding, positive numeric. Default is 1e-6.
#'
#' @return A list of 1. Parameters of baseline processes, 2. Mixing weights, 3. Auxiliary values for computation.
#' @export
#'
#' @examples
#' compute_baseline(0.01, 1, 2)
compute_baseline <- function(alpha,
                             delta_lower,
                             delta_upper,
                             psi_star = function(x){x^2/2},
                             psi_star_div = function(x){x},
                             psi_star_inv = function(y){sqrt(2*y)},
                             v_min = 1,
                             k_max = 100,
                             tol = 1e-6){
  # Type checks
  if (!(alpha > 0 | alpha < 1)) stop("alpha must be a number in (0,1).")
  if (!(delta_lower > 0 & delta_upper >=  delta_lower)){
    stop("delta_lower and delta_upper must be positive with delta_lower <= delta_upper.")
  }
  if (psi_star(0) != 0) stop("psi_star must be zero at x = 0.")
  if (psi_star_div(0) != 0) stop("psi_star_div must be zero at x = 0.")
  if (!(v_min >= 0)) stop("v_min must be non-negative.")
  k_max_raw <- k_max
  k_max <- as.integer(k_max_raw)
  if (k_max != k_max_raw) warning("k_max is coverted to an integer.")
  if (k_max < 1) stop("k_max must be larger than or equal to 1.")

  # Compute constants
  log_one_over_alpha <- log(1/alpha)
  d_l <- psi_star(delta_lower)
  d_u <- psi_star(delta_upper)
  ratio <- d_u / d_l
  lambda_l <- psi_star_div(delta_lower)
  lambda_u <- psi_star_div(delta_upper)

  # If delta_lower is small enough or equal to delta_upper
  # Return trivial single baseline
  if (log_one_over_alpha <= v_min * d_l | delta_lower == delta_upper){
    baseline_list <- list(
      lambda = lambda_l,
      omega = 1,
      g_alpha = log_one_over_alpha,
      k_alpha = 0,
      eta_alpha = 1,
      w = alpha
        )
    return(baseline_list)
  }

  # Compute the threshold g_alpha
  log_f <- function(g){
    k_vec <- 1:k_max
    log_f_val <- sapply(k_vec, function(k) log(k) - g * ratio^(-1/k))
    return(min(log_f_val))
  }

  log_f_with_exp <- function(g){
    exp_vec <- c(-g, log_f(g))
    return(matrixStats::logSumExp(exp_vec))
  }

  log_f_u <- log_f(v_min * d_u)
  if (log_f_u <= -log_one_over_alpha){
    root_out <- stats::uniroot(function(g){log_f(g) + log_one_over_alpha},
                        c(log_one_over_alpha, v_min * d_u), tol = tol)
  } else {
    root_out <- stats::uniroot(function(g){log_f_with_exp(g) + log_one_over_alpha},
                        c(v_min * d_u, ratio * log(2/alpha)), tol = tol)
  }
  g_alpha <- root_out$root

  # Compute the number of non-trival baselines k_alpha
  k_vec <- 1:k_max
  log_f_val <- sapply(k_vec, function(k) log(k) - g_alpha * ratio^(-1/k))
  k_alpha <- which.min(log_f_val)

  # Compute the spacing parameter eta_alpha
  eta_alpha <- ratio^(1/k_alpha)

  # Compute lambdas and mixing weights
  if (g_alpha > v_min * d_u){
    # Compute lambdas
    lambda_vec <- numeric(k_alpha + 1)
    lambda_vec[1] <- lambda_u
    lambda_vec[k_alpha +  1] <- lambda_l
    if (k_alpha >= 2){
      k_candidate <- seq(1, k_alpha - 1)
      delta_vec <- sapply(k_candidate, function(k) psi_star_inv(d_u / eta_alpha^k))
      lambda_vec[-c(1,k_alpha+1)] <- sapply(delta_vec, psi_star_div)
    }
    # Compute weights
    omega_vec <- c(exp(-g_alpha), rep(exp(-g_alpha / eta_alpha), k_alpha))
  } else {
    # Compute lambdas
    lambda_vec <- numeric(k_alpha)
    lambda_vec[k_alpha] <- lambda_l
    if (k_alpha >= 2){
      k_candidate <- seq(1, k_alpha - 1)
      delta_vec <- sapply(k_candidate, function(k) psi_star_inv(d_u / eta_alpha^k))
      lambda_vec[-k_alpha] <- sapply(delta_vec, psi_star_div)
    }
    # Compute weights
    omega_vec <- rep(exp(-g_alpha / eta_alpha), k_alpha)
  }

  # Normalize weights
  w <- sum(omega_vec)
  omega_normal_vec <- omega_vec / w

  # Collect all computed parameters
  baseline_list <- list(
    lambda = lambda_vec,
    omega = omega_normal_vec,
    g_alpha = g_alpha,
    k_alpha = k_alpha,
    eta_alpha = eta_alpha,
    w = w
  )
  return(baseline_list)
}
