#' Plot observed e-value / e-detectors
#'
#' @param edcp_out edcp_out object
#' @param n A numeric vector of time indices
#' @param add Indicator whether to make an overlay on the existing plot or not
#'
#' @export
#'
plot.edcp_out <- function(edcp_out,
                          n = seq_along(edcp_out$log_mix_e_vec),
                          draw_detect_line = TRUE,
                          add = FALSE,
                          lty = 1,...) {
  if (edcp_out$edcp_obj$is_test) {
    ylab_ = "log e-value"
  } else {
    ylab_ = "log e-detector"
  }
  if (!add) {
    plot(
      n,
      edcp_out$log_mix_e_vec,
      type = "l",
      xlab = "n",
      ylab = ylab_,
      lty = lty,
      ...
    )
    graphics::abline(h = edcp_out$edcp_obj$log_one_over_alpha, col = 2)
  } else {
    graphics::lines(n,
                    edcp_out$log_mix_e_vec,
                    lty = lty,
                    ...)
  }
  if (draw_detect_line & !is.infinite(edcp_out$stopped_ind)){
    graphics::abline(v = n[edcp_out$stopped_ind], lty = lty + 1, ...)
  }
}
