#' Plot observed e-value / e-detectors
#'
#' @param x EDCP_RUNNED object
#' @param n A numeric vector of time indices
#' @param draw_detect_line If true then draw a vertical line to indicate the detection time.
#' @param add Indicator whether to make an overlay on the existing plot or not
#' @param ... Graphical parameters can be passed except \code{type} and \code{ylab}, which are pre-defined.
#'
#' @export
#'
plot.EDCP_RUNNED <- function(x,
                          n = seq_along(x$log_mix_e_vec),
                          draw_detect_line = TRUE,
                          add = FALSE,
                          ...) {
  if (x$edcp_obj$is_test) {
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
    graphics::abline(h = x$edcp_obj$log_one_over_alpha, col = 2)
  } else {
    graphics::lines(n,
                    x$log_mix_e_vec,
                    ...)
  }
  if (draw_detect_line & !is.infinite(x$stopped_ind)){
    graphics::abline(v = n[x$stopped_ind], ...)
  }
}

#' Print EDCP model summary
#'
#' @param x EDCP object
#' @param ... Placeholder for formatting options - not yet implemented.
#'
#' @export
#'
print.EDCP <- function(x,...) {
  # TODO::Implement detailed print for exp and bounded sub models.
  cat("EDCP Model:\n",...)
  cat("- alpha: ", x$alpha, "\n")
  cat("- is_test: ", x$is_test, "\n")
  cat("- Num. of mixing components: ", length(x$log_base_fn_list), "\n")
  cat("- Obs. have been passed: ", x$n, "\n")
  invisible(x)
}

#' Print EDCP run summary
#'
#' @param x edcp_out object
#' @param ... Placeholder for formatting options - not yet implemented.
#'
#' @export
#'
print.EDCP_RUNNED <- function(x,...) {
  cat("EDCP Run:\n")
  cat("- Num. of new obs: ", length(x$log_mix_e_vec), "\n")
  cat("- Stopped index: ", x$stopped_ind, "\n\n")
  print(x$edcp_obj)
  invisible(x)
}




