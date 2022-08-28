#' Build stcp model.
#'
#' Build stcp model to detect up-crossing mean-shift.
#'
#' @param alpha ARL parameter in (0,1)
#' @param m_pre Upper bound of mean of pre-change observations
#' @param is_test A Boolean to indicate whether this model is for a sequential test or not.
#' @param omega Mixing weights for baseline processes.
#' @param lambda Parameters of baseline processes.
#' @param log_base_fn_generator Function factory generating log baseline functions.
#' @param ... arguments to be passed to \code{log_base_fn_generator}.
#'
#' @return \code{stcp} object of 1. Model parameters, 2. Memory for log e-detectors / values.
#' @export
#'
build_stcp <- function(alpha,
                       m_pre,
                       is_test,
                       omega,
                       lambda,
                       log_base_fn_generator,
                       ...) {

  if (length(omega) != length(lambda)) {
    stop("Number of weights and model parameters are not matched.")
  }

  if (min(omega) < 0) stop("Mixing weights must be nonnegative.")

  # Compute e-detectors
  log_base_fn_list <- sapply(
    lambda,
    log_base_fn_generator,
    ...
  )

    # Initialize log_e_vec
  if (is_test) {
    log_e_vec <- numeric(length(log_base_fn_list))
  } else {
    log_e_vec <- rep(-Inf, length(log_base_fn_list))
  }

  out <- list(
    # Model parameters
    omega = omega,
    log_base_fn_list =  log_base_fn_list,
    log_one_over_alpha = log(1/alpha),
    alpha = alpha,
    m_pre = m_pre,
    is_test = is_test,
    family_name = "Custom",
    # Memory for log e-detectors / values
    log_e_vec = log_e_vec,
    n = 0,
    # Auxiliaries for debugging
    lambda = lambda
  )
  class(out) <- c("stcp")
  return(out)
}

#' Combine two stcp objects for a given weight.
#'
#' Combine two stcp objects for a given weight. Only \code{stcp} objects of the same \code{is_test} can be combined.
#' Once combined, any additional details in input objects will be discarded.
#'
#' @param obj1 First \code{stcp} object to combine.
#' @param obj2 Second \code{stcp} object to combine.
#' @param w Weight scalar to determine the importance of the first \code{stcp} object. Default is 0.5.
#'
#' @return Combined \code{stcp} object of 1. Model parameters, 2. Memory for log e-detectors / values.
#' @export
#'
combine_stcp <- function(obj1, obj2, w = 0.5) {
  if (w <= 0 || w >= 1 ) stop("w must be a number between (0,1)")
  if (obj1$is_test != obj2$is_test) {
    stop("e-value (is_test == TRUE) and e-detector (is_test != FALSE) cannot be combined")
  }

  combined_omega <- c(w * obj1$omega, (1-w) * obj2$omega)
  combined_alpha <- w * obj1$alpha + (1-w) * obj2$alpha
  combined_log_base_fn_list <- c(obj1$log_base_fn_list, obj2$log_base_fn_list)
  combined_m_pre <- c(obj1$m_pre, obj2$m_pre)
  combined_family_name <- c(obj1$family_name, obj2$family_name)
  combined_log_e_vec <- c(obj1$log_e_vec, obj2$log_e_vec)
  combined_lambda <- c(obj1$lambda, obj2$lambda)

  out <- list(
    # Model parameters
    omega = combined_omega,
    log_base_fn_list =  combined_log_base_fn_list,
    log_one_over_alpha = log(1/combined_alpha),
    alpha = combined_alpha,
    m_pre = combined_m_pre,
    is_test = obj1$is_test,
    family_name = combined_family_name,
    # Memory for log e-detectors / values
    log_e_vec = combined_log_e_vec,
    n = 0,
    # Auxiliaries for debugging
    lambda = combined_lambda
  )
  class(out) <- c("stcp")
  return(out)

}


#' Build stcp model for simple exponential e-detectors.
#'
#' Build stcp model for simple exponential e-detectors to detect up-crossing mean-shift.
#'
#' @param m_pre Upper bound of mean of pre-change observations
#' @param is_test A Boolean to indicate whether this model is for a sequential test or not.
#' @param s_fn R function that compute the sum process given an observation.
#' @param v_fn R function that compute the variance process given an observation.
#' @inheritParams compute_baseline
#'
#' @return \code{stcp} object of 1. Model parameters, 2. Memory for log e-detectors / values, 3. Auxiliaries
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

  # Make stcp object
  stcp_obj <- build_stcp(
    alpha = alpha,
    m_pre = m_pre,
    is_test = is_test,
    omega = base_param$omega,
    lambda = base_param$lambda,
    log_base_fn_generator =generate_log_base_fn,
    psi_fn = base_param$psi_fn_list$psi,
    s_fn = s_fn,
    v_fn = v_fn
  )

  stcp_obj$family_name <- base_param$psi_fn_list$family_name

  out <- c(
    stcp_obj,
    list(
      g_alpha = base_param$g_alpha,
      # Input Details
      delta_lower = delta_lower,
      delta_upper = delta_upper,
      s_fn = s_fn,
      v_fn = v_fn,
      v_min = v_min
    )
  )
  class(out) <- c("stcp_exp", "stcp")
  return(out)
}

#' Build stcp model for bounded random variables.
#'
#' Build stcp model for bounded random variables to detect up-crossing mean-shift.
#'
#' @param m_pre Upper bound of mean of pre-change observations. Must be between \code{bound_lower} and \code{bound_upper}.
#' @param is_test A Boolean to indicate whether this model is for a sequential test or not.
#' @param bound_lower Lower bound of observations.
#' @param bound_upper Upper bound of observations.
#' @param is_test A Boolean to indicate whether this model is for a sequential test or not.
#' @param k_max Positive integer to determine the maximum number of baselines. Default is \code{1000}.
#' @param var_lower Lower bounds of variance of post-change observations. Default is \code{0}.
#' @param var_upper Upper bounds of variance of post-change observations. Default is \code{(bound_upper - bound_lower)^2 / 4}.
#' @param delta_lower_sub_E Lower bound of target Delta in sub-E scale.
#' If both upper and lower sub-E parameters are not null then these parameters will be directly used to build baseline process.
#' @param delta_upper_sub_E Upper bound of target Delta in sub-E scale.
#'  If both upper and lower sub-E parameters are not null then this parameter will be directly used to build baseline process.
#' @inheritParams compute_baseline
#'
#' @return \code{stcp} object of 1. Model parameters, 2. Memory for log e-detectors / values, 3. Auxiliaries
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
                               var_upper = (bound_upper - bound_lower)^2 / 4,
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

  if (abs(delta_upper + m_pre - bound_upper) < 1e-6) {
    # If delta_upper is equal to bound_upper - m_pre, the upper delta is given by var_upper = 0.
    var_upper <- 0
  }

  # Compute parameters
  bound_range <- bound_upper - bound_lower
  m_scaled <- (m_pre - bound_lower) / bound_range # scaled m_pre
  d_l <- delta_lower / bound_range  # scaled delta_lower
  d_u <- delta_upper / bound_range # scaled_delta_upper
  v_l <- var_lower / bound_range^2 # scaled var_lower
  v_u <- var_upper / bound_range^2 # scaled var_upper

  if (!is.null(delta_lower_sub_E) & !is.null(delta_upper_sub_E)) {
    delta_lower_val <- delta_lower_sub_E
    delta_upper_val <-  delta_upper_sub_E
    var_lower <- 0
    var_upper <- 0
    delta_lower <- bound_range * m_scaled / (delta_lower_val * delta_upper_val^2)^(1/3)
    delta_upper <- bound_range * m_scaled / (delta_lower_val^2 * delta_upper_val)^(1/3)

  } else {
    delta_lower_val <-  m_scaled * d_l / (v_u + d_u ^ 2)
    delta_upper_val <-  m_scaled * d_u / (v_l + d_l ^ 2)
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

  # Make stcp object
  stcp_obj <- build_stcp(
    alpha = alpha,
    m_pre = m_pre,
    is_test = is_test,
    omega = base_param$omega,
    lambda = base_param$lambda,
    log_base_fn_generator = generate_log_bounded_base_fn,
    m = m_pre,
    bound_lower = bound_lower
  )

  stcp_obj$family_name <- "Bounded (sub-E based)"

  out <- c(
    stcp_obj,
    list(
      g_alpha = base_param$g_alpha,
      # Input Details
      delta_lower = delta_lower,
      delta_upper = delta_upper,
      bound_lower = bound_lower,
      bound_upper = bound_upper,
      var_lower = var_lower,
      var_upper = var_upper
    )
  )
  class(out) <- c("stcp_bounded", "stcp")
  return(out)
}
