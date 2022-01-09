#' Build a generator of always-valid confidence intervals for simple exponential e-values.
#'
#' Build a generator of always-valid confidence intervals for simple exponential e-valuess
#' given target range of numbers of observations.
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
