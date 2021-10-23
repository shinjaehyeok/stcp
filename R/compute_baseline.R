#' Compute baseline processes
#'
#' Compuate parameters to build baselin processes.
#'
#' @param alpha ARL parameter in (0,1)
#' @param delta_lower Lower bound of target Delta. It must have the same sign of \code{delta_upper}.
#' @param delta_upper Upper bound of target Delta. It must have the same sign of \code{delta_lower}.
#' @param psi_star R function which computes the convex conjugate of the psi function. Default is \code{function(x) x^2/2}.
#' @param psi_star_div R function which computes the derivative of \code{psi_star}. Default is \code{function (x) x}.
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
                             v_min = 1,
                             k_max = 100,
                             tol = 1e-6){
  # Check the sign of delta bounds.
  # TODO Add to Rstudio Cloud
  print("Under Construction")
  print(psi_star)
  print(psi_star_div)
}
