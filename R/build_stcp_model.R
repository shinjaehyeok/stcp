#' Build STCP model for simple exponential e-detectors.
#'
#' Build STCP model for simple exponential e-detectors to detect up-crossing mean-shift.
#'
#' @param m_pre Upper bound of mean of pre-change observations
#' @param is_test A Boolean to indicate whether this model is for a sequential test or not.
#' @param s_fn R function that compute the sum process given an observation.
#' @param v_fn R function that compute the variance process given an observation.
#' @inheritParams compute_baseline
#'
#' @return \code{STCP} object of 1. Model parameters, 2. Memory for log e-detectors / values, 3. Auxiliaries
#' @export
#'
build_stcp_exp <- function(alpha,
                           m_pre,
                           delta_lower,
                           delta_upper,
                           is_test = FALSE,
                           psi_fn_list = generate_sub_G_fn(sig = 1),
                           s_fn = function(x) {
                             x - m_pre
                           },
                           v_fn = function(x) {
                             1
                           },
                           v_min = 1,
                           k_max = 200,
                           tol = 1e-6) {
  # Compute parameters
  base_param <- compute_baseline(alpha,
                                 delta_lower,
                                 delta_upper,
                                 psi_fn_list,
                                 v_min,
                                 k_max,
                                 tol)

  # Compute e-detectors
  log_base_fn_list <- sapply(
    base_param$lambda,
    generate_log_base_fn,
    psi_fn = base_param$psi_fn_list$psi,
    s_fn = s_fn,
    v_fn = v_fn
  )

  # Initialize log_e_vec
  if (is_test) {
    log_e_vec <- numeric(length(log_base_fn_list))
  } else {
    log_e_vec <- rep(-Inf, length(log_base_fn_list))
  }

  out <- list(
    # Model parameters
    omega = base_param$omega,
    log_base_fn_list =  log_base_fn_list,
    log_one_over_alpha = log(1/alpha),
    alpha = alpha,
    m_pre = m_pre,
    is_test = is_test,
    family_name = base_param$psi_fn_list$family_name,
    # Memory for log e-detectors / values
    log_e_vec = log_e_vec,
    n = 0,
    # Auxiliaries for debugging
    lambda = base_param$lambda,
    g_alpha = base_param$g_alpha,
    # Input Details
    delta_lower = delta_lower,
    delta_upper = delta_upper,
    s_fn = function(x) {
      x - m_pre
    },
    v_fn = function(x) {
      1
    },
    v_min = 1
  )
  class(out) <- c("SimpleExp", "STCP")
  return(out)
}

#' Build STCP model for bounded random variables.
#'
#' Build STCP model for bounded random variables to detect up-crossing mean-shift.
#'
#' @param m_pre Upper bound of mean of pre-change observations. Must be between \code{bound_lower} and \code{bound_upper}.
#' @param is_test A Boolean to indicate whether this model is for a sequential test or not.
#' @param bound_lower Lower bound of observations.
#' @param bound_upper Upper bound of observations.
#' @param is_test A Boolean to indicate whether this model is for a sequential test or not.
#' @param k_max Positive integer to determine the maximum number of baselines. Default is \code{1000}.
#' @param var_lower Lower bounds of variance of scaled post-change observations. Default is \code{0}.
#' @param var_upper Upper bounds of variance of scaled post-change observations. Default is \code{0.25}.
#' @param delta_lower_sub_E Lower bound of target Delta in sub-E scale.
#' If both upper and lower sub-E parameters are not null then these parameters will be directly used to build baseline process.
#' @param delta_upper_sub_E Upper bound of target Delta in sub-E scale.
#'  If both upper and lower sub-E parameters are not null then this parameter will be directly used to build baseline process.
#' @inheritParams compute_baseline
#'
#' @return \code{STCP} object of 1. Model parameters, 2. Memory for log e-detectors / values, 3. Auxiliaries
#' @export
#'
#' @examples
#' build_stcp_bounded(0.01, 0.5, 0.1, 0.4)
build_stcp_bounded <- function(alpha,
                               m_pre,
                               delta_lower,
                               delta_upper = bound_upper - m_pre,
                               is_test = FALSE,
                               bound_lower = 0,
                               bound_upper = 1,
                               k_max = 1000,
                               tol = 1e-6,
                               var_lower = 0,
                               var_upper = 0.25,
                               delta_lower_sub_E = NULL,
                               delta_upper_sub_E = NULL) {
  if (!(m_pre > bound_lower & m_pre < bound_upper)) {
    stop("m_pre must be between bound_lower and bound_upper.")
  }

  if (!(delta_lower > 0 & delta_upper >=  delta_lower)) {
    stop("delta_lower and delta_upper must be positive with delta_lower <= delta_upper.")
  }

  if (!(var_lower >= 0 & var_upper >=  var_lower)) {
    stop("var_lower and var_upper must be positive with var_lower <= var_upper.")
  }

  # Compute parameters
  bound_range <- bound_upper - bound_lower
  m <- (m_pre - bound_lower) / bound_range # scaled m_pre
  d_l <- delta_lower / bound_range  # scaled delta_lower
  d_u <- delta_upper / bound_range # scaled_delta_upper

  if (!is.null(delta_lower_sub_E) & !is.null(delta_upper_sub_E)) {
    delta_lower_val <- delta_lower_sub_E
    delta_upper_val <-  delta_upper_sub_E
    var_lower <- 0
    var_upper <- 0
    delta_lower <- bound_range * m / (delta_lower_val * delta_upper_val^2)^(1/3)
    delta_upper <- bound_range * m / (delta_lower_val^2 * delta_upper_val)^(1/3)

  } else {
    delta_lower_val <-  m * d_l / (var_upper + d_u ^ 2)
    delta_upper_val <-  m * d_u / (var_lower + d_l ^ 2)
  }

  base_param <- compute_baseline(
    alpha,
    delta_lower = delta_lower_val,
    delta_upper = delta_upper_val,
    psi_fn_list = generate_sub_E_fn(),
    v_min = 0,
    k_max,
    tol
  )

  # Compute e-detectors
  log_base_fn_list <- sapply(
    base_param$lambda,
    generate_log_bounded_base_fn,
    m = m_pre,
    bound_lower = bound_lower
  )

  out <- list(
    # Model parameters
    omega = base_param$omega,
    log_base_fn_list =  log_base_fn_list,
    log_one_over_alpha = log(1/alpha),
    alpha = alpha,
    m_pre = m_pre,
    is_test = is_test,
    family_name = "Bounded (sub-E based)",
    # Memory for log e-detectors / values
    log_e_vec = numeric(length(log_base_fn_list)),
    n = 0,
    # Auxiliaries for debugging
    lambda = base_param$lambda,
    g_alpha = base_param$g_alpha,
    # Input Details
    delta_lower = delta_lower,
    delta_upper = delta_upper,
    bound_lower = bound_lower,
    bound_upper = bound_upper,
    var_lower = var_lower,
    var_upper = var_upper
  )
  class(out) <- c("Bounded", "STCP")
  return(out)
}