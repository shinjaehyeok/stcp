# Helper functions

#' Pre-defined psi_star functions for sub-Gaussian family.
#'
#' @param sig The sigma parameter of the sub-Gaussian family (default = 1).
#'
#' @return A list of pre-defined psi_star functions for sub-Gaussian family.
#'
#' @export
generate_sub_G_fn <- function(sig = 1){
  force(sig)
  G_fn_list <- list(
    family_name = "sub-G",
    is_psi_depend_on_m = FALSE,
    psi = function(x){x^2 * sig^2 / 2},
    psi_star = function(x){x^2 / 2 / sig^2},
    psi_star_div = function(x){x / sig^2},
    psi_star_inv = function(x){sig * sqrt(2 * x)}
  )
  return(G_fn_list)
}


#' Pre-defined psi_star functions for sub-Bernoulli family.
#'
#' @param p The boundary of mean space of the pre-change distributions (default = 0.5).
#' @param tol A small number to use truncate \code{p} into (\code{tol}, \code{1-tol}) range.
#'
#' @return A list of pre-defined psi_star functions for sub-Bernoulli family.
#'
#' @export
generate_sub_B_fn <- function(p = 0.5, tol = 1e-6){
  force(p)
  if (p <= 0){
    p <- tol
  }
  if (p >= 1){
    p <- 1- tol
  }
  B_fn_list <- list(
    family_name = "sub-B",
    is_psi_depend_on_m = TRUE,
    psi = function(x){
      log(1-p + p * exp(x)) - x * p
    },
    psi_star = function(x){
      d <- x + p
      if (d >= 1) {
        val <- log(1/p)
      } else if (d <= 0) {
        val <- log(1/(1-p))
      } else {
        val <- d * log(d / p) + (1-d) * log((1-d) / (1-p))
      }
      return(val)
    },
    psi_star_div = function(x){
      d <- x + p
      if (d >= 1) {
        return(Inf)
      } else if (d <= 0) {
        return(-Inf)
      } else {
        div <- log(d * (1-p) / p / (1-d))
        return(div)
      }
    },
    psi_star_inv = function(y){
      right_end <- 1-p
      bound <- B_fn_list$psi_star(right_end)
      if (y >= bound) return(right_end)
      f <- function(x) B_fn_list$psi_star(x) - y
      x <- stats::uniroot(f,
                          c(0, 1-p),
                          tol = 1e-12)$root
      return(x)
    }
  )
  return(B_fn_list)
}


#' Pre-defined psi_star functions for sub-exponential family.
#'
#' @return A list of pre-defined psi_star functions for sub-exponential family.
#'
#' @export
generate_sub_E_fn <- function(){
  E_fn_list <- list(
    family_name = "sub-E",
    is_psi_depend_on_m = FALSE,
    psi = function(x){
      -log(1-x) - x
    },
    psi_star = function(x){
      val <- x - log(1+x)
      return(val)
    },
    psi_star_div = function(x){
      div <- x / (1+x)
      return(div)
    },
    psi_star_inv = function(y, max_bound = 1000){
      f <- function(x) E_fn_list$psi_star(x) - y
      if (f(max_bound) <= 0) return(max_bound)
      x <- stats::uniroot(f,
                          c(0, max_bound),
                          tol = 1e-12)$root
      return(x)
    }
  )
  return(E_fn_list)
}
