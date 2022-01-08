#' Plot observed e-value / e-detectors
#'
#' @param x edcp_out object
#' @param y A numeric vector of time indices
#' @param draw_detect_line If true then draw a vertical line to indicate the detection time.
#' @param add Indicator whether to make an overlay on the existing plot or not
#' @param ... Graphical parameters can be passed except \code{type} and \code{ylab}, which are pre-defined.
#'
#' @export
#'
plot.edcp_out <- function(x,
                          y = seq_along(x$log_mix_e_vec),
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
      y,
      x$log_mix_e_vec,
      type = "l",
      ylab = ylab_,
      ...
    )
    graphics::abline(h = x$edcp_obj$log_one_over_alpha, col = 2)
  } else {
    graphics::lines(y,
                    x$log_mix_e_vec,
                    ...)
  }
  if (draw_detect_line & !is.infinite(x$stopped_ind)){
    graphics::abline(v = y[x$stopped_ind], ...)
  }
}
