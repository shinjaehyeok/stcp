#' Build a generator of lower bound of always-valid confidence intervals for simple exponential e-values.
#'
#' Build a generator of lower bound of always-valid confidence intervals
#' for simple exponential e-values given target range of numbers of observations.
#'
#' @param psi_fn_list_generator Function operator generating \code{psi_fn_list}
#' @inheritParams compute_baseline_for_sample_size
#'
#' @return A list of 1. Model parameters, 2. Memory for log e-detectors / values, 3. Auxiliaries
#' @export
#'
build_ci_exp <- function(alpha,
                         n_upper,
                         n_lower,
                         psi_fn_list_generator = generate_sub_G_fn,
                         v_min = 1,
                         k_max = 200,
                         tol = 1e-6) {

  psi_fn_list <- psi_fn_list_generator()
  baseline_obj <- compute_baseline_for_sample_size(alpha,
                                                   n_upper,
                                                   n_lower,
                                                   psi_fn_list,
                                                   v_min,
                                                   k_max,
                                                   tol)
  ci_helper <- generate_ci_helper(
    baseline_obj$alpha,
    baseline_obj$omega,
    baseline_obj$lambda,
    psi_fn_list_generator,
    psi_fn_list$is_psi_depend_on_m
  )

  out <- list(ci_helper = ci_helper,
              baseline_obj = baseline_obj,
              family_name = psi_fn_list$family_name,
              # Input Details
              n_upper = n_upper,
              n_lower = n_lower,
              v_min = v_min
              )

  class(out) <- "EDCP_CI"
  return(out)
}

#' Build a brute-force generator of lower bound of always-valid confidence intervals for bounded random variables.
#'
#' Build a brute-force generator of lower bound of always-valid confidence intervals
#' for bounded random variables given target range of numbers of observations.
#'
#' @inheritParams compute_baseline_for_sample_size
#' @param  bound_lower  Lower bound of observations.
#'
#' @return A list of 1. Model parameters, 2. Memory for log e-detectors / values, 3. Auxiliaries
#' @export
#'
build_bf_ci_bounded <- function(alpha,
                                n_upper,
                                n_lower,
                                bound_lower,
                                k_max = 200,
                                tol = 1e-6) {

  baseline_obj <- compute_baseline_for_sample_size(alpha,
                                                 n_upper,
                                                 n_lower,
                                                 generate_sub_E_fn(),
                                                 v_min = 0,
                                                 k_max = k_max,
                                                 tol = tol)

  bf_ci_helper <- function(m, x_vec, tol = 1e-6) {
    if (m <= bound_lower) {
      m <- bound_lower + tol
    }
    # Compute e-detectors
    log_base_fn_list <- sapply(
      baseline_obj$lambda,
      generate_log_bounded_base_fn,
      m = m,
      bound_lower = bound_lower
    )

    edcp_model <- list(
      omega = baseline_obj$omega,
      log_base_fn_list =  log_base_fn_list,
      log_one_over_alpha = log(1/alpha),
      alpha = alpha,
      is_test = TRUE,
      family_name = "Bounded (sub-E based)",
      log_e_vec = numeric(length(log_base_fn_list)),
      n = 0
    )
    class(edcp_model) <- c("Bounded", "EDCP")

    e_val <- run_edcp(x_vec, edcp_model)
    return(e_val$log_mix_e_vec[length(e_val$log_mix_e_vec)] - e_val$edcp_obj$log_one_over_alpha)
  }

  out <- list(bf_ci_helper = bf_ci_helper,
              baseline_obj = baseline_obj,
              family_name = "Bounded (sub-E based)",
              # Input Details
              n_upper = n_upper,
              n_lower = n_lower,
              bound_lower = bound_lower
  )

  class(out) <- "EDCP_BF_CI"
  return(out)
}
