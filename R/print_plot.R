#' Plot observed e-value / e-detectors
#'
#' @param x stcp_run object
#' @param n A numeric vector of time indices
#' @param draw_detect_line If true then draw a vertical line to indicate the detection time.
#' @param add Indicator whether to make an overlay on the existing plot or not
#' @param ... Graphical parameters can be passed except \code{type} and \code{ylab}, which are pre-defined.
#'
#' @export
#'
plot.stcp_run <- function(x,
                             n = seq_along(x$log_mix_e_vec),
                             draw_detect_line = TRUE,
                             add = FALSE,
                             ...) {
  if (x$stcp_obj$is_test) {
    ylab_ = "log e-value"
  } else {
    ylab_ = "log e-detector"
  }
  if (!add) {
    plot(
      n,
      x$log_mix_e_vec,
      type = "l",
      ylab = ylab_,
      ...
    )
    graphics::abline(h = x$stcp_obj$log_one_over_alpha, col = 2)
  } else {
    graphics::lines(n,
                    x$log_mix_e_vec,
                    ...)
  }
  if (draw_detect_line & !is.infinite(x$stopped_ind)){
    graphics::abline(v = n[x$stopped_ind], ...)
  }
}

#' Print stcp model summary
#'
#' @param x stcp object
#' @param ... Placeholder for formatting options - not yet implemented.
#'
#' @export
#'
print.stcp <- function(x,...) {
  cat("stcp Model:\n",...)
  cat("- Family: ", x$family_name, "\n")
  cat("- alpha: ", x$alpha, "\n")
  cat("- m_pre: ", x$m_pre, "\n")
  cat("- is_test: ", x$is_test, "\n")
  cat("- Num. of mixing components: ", length(x$log_base_fn_list), "\n")
  cat("- Obs. have been passed: ", x$n, "\n")
  if (class(x)[1] == "stcp_exp") {
    cat("- delta_lower: ", x$delta_lower, "\n")
    cat("- delta_upper: ", x$delta_upper, "\n\n")
    cat("- s_fn: ")
    print(x$s_fn)
    cat("\n")
    cat("- v_fn: ")
    print(x$v_fn)
    cat("\n")
    cat("- v_min: ", x$v_min, "\n")
  } else if (class(x)[1] == "stcp_bounded") {
    cat("- delta_lower: ", x$delta_lower, "\n")
    cat("- delta_upper: ", x$delta_upper, "\n")
    cat("- bound_lower: ", x$bound_lower, "\n")
    cat("- bound_upper: ", x$bound_upper, "\n")
    cat("- var_lower: ", x$var_lower, "\n")
    cat("- var_upper: ", x$var_upper, "\n")
  }
  invisible(x)
}

#' Print stcp run summary
#'
#' @param x stcp_out object
#' @param ... Placeholder for formatting options - not yet implemented.
#'
#' @export
#'
print.stcp_run <- function(x,...) {
  cat("stcp Run:\n")
  cat("- Num. of new obs: ", length(x$log_mix_e_vec), "\n")
  cat("- Stopped index: ", x$stopped_ind, "\n\n")
  print(x$stcp_obj)
  invisible(x)
}

#' Print stcp_ci model summary
#'
#' @param x stcp_ci object
#' @param ... Placeholder for formatting options - not yet implemented.
#'
#' @export
#'
print.stcp_ci <- function(x,...) {
  cat("stcp_ci Model:\n",...)
  cat("- Family: ", x$family_name, "\n")
  cat("- alpha: ", x$baseline_obj$alpha, "\n")
  cat("- Num. of mixing components: ", length(x$baseline_obj$lambda), "\n")
  cat("- n_upper: ", x$n_upper, "\n")
  cat("- n_lower: ", x$n_lower, "\n\n")
  cat("- s_fn: ")
  print(function(x){x})
  cat("\n")
  cat("- v_fn: ")
  print(function(x){1})
  cat("\n")
  cat("- v_min: ", x$v_min, "\n")
  invisible(x)
}

#' Print stcp_bf_ci model summary
#'
#' @param x stcp_bf_ci object
#' @param ... Placeholder for formatting options - not yet implemented.
#'
#' @export
#'
print.stcp_bf_ci <- function(x,...) {
  cat("stcp_bf_ci Model:\n",...)
  cat("- Family: ", x$family_name, "\n")
  cat("- alpha: ", x$baseline_obj$alpha, "\n")
  cat("- Num. of mixing components: ", length(x$baseline_obj$lambda), "\n")
  cat("- n_upper: ", x$n_upper, "\n")
  cat("- n_lower: ", x$n_lower, "\n")
  if(!is.null(x$bound_lower)){
    cat("- bound_lower: ", x$bound_lower, "\n")
  }
  if(!is.null(x$v_min)){
    cat("- v_min: ", x$v_min, "\n")
  }

  invisible(x)
}


