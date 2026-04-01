# =============================================================================
# fapa_plots.R
# Visualisation functions for FAPA results.
# =============================================================================


#' Plot a core profile with BCa confidence intervals
#'
#' Displays the \eqn{i}th core profile from a \code{\link{fapa_bca}} result,
#' split at a variable boundary (e.g., pre vs post) with BCa CI bands.
#' Variables before \code{split_at} are drawn in red (dashed), variables from
#' \code{split_at + 1} onward in blue (solid).
#'
#' @param bca      A list returned by \code{\link{fapa_bca}}.
#' @param i        Integer. Which core profile to plot. Default \code{1}.
#' @param split_at Integer. Index at which to switch colour/line-type.
#'   Default \code{11} (11 pre + 11 post EDI-2 subscales).
#' @param main     Character. Plot title. Default auto-generated.
#' @param ylab     Character. Y-axis label.
#'
#' @return Invisibly returns \code{NULL}. Called for its side-effect.
#' @export
plot_fapa_core <- function(bca, i = 1, split_at = 11,
                            main = NULL,
                            ylab = "Core-Profile Coordinate") {
  if (is.null(main))
    main <- sprintf("FAPA Core Profile %d  (95%% BCa CI)", i)
  s  <- bca$ci[[i]]
  op <- graphics::par(mar = c(2.1, 4.1, 2.5, 1)); on.exit(graphics::par(op))
  yr <- range(s$BCaLower, s$BCaUpper, na.rm = TRUE)

  idx_before <- seq_len(split_at)
  idx_after  <- seq(split_at + 1L, nrow(s))

  plot(s$Ori[idx_before], type = "b", lty = 2, lwd = 2,
       col = "red", ylim = yr, xlab = "", ylab = ylab,
       main = main, xaxt = "n")
  graphics::axis(1, at = seq_along(idx_before),
                 labels = bca$varnames[idx_before], las = 2, cex.axis = 0.6)
  graphics::lines(s$BCaLower[idx_before], lty = 3, col = "red")
  graphics::lines(s$BCaUpper[idx_before], lty = 3, col = "red")
  graphics::abline(h = 0, lty = 3, col = "grey50")
  graphics::lines(s$Ori[idx_after], type = "b", lwd = 2, col = "blue")
  graphics::lines(s$BCaLower[idx_after], lty = 3, col = "blue")
  graphics::lines(s$BCaUpper[idx_after], lty = 3, col = "blue")
  graphics::legend("bottom", inset = -0.02, xpd = NA, bty = "n",
                   legend = c("Before", "After"),
                   col = c("red", "blue"), lty = c(2, 1), lwd = 2, ncol = 2)
  invisible(NULL)
}


#' Scree plot for Horn's parallel analysis
#'
#' Plots observed \eqn{\sigma^2_k} versus the random 95th-percentile reference
#' line, with a vertical cut at the retention boundary.
#'
#' @param pa   A list returned by \code{\link{fapa_pa}}.
#' @param main Character. Plot title.
#'
#' @return Invisibly returns \code{NULL}.
#' @export
plot_pa_scree <- function(pa,
                           main = "Horn's Parallel Analysis \u2014 Scree") {
  K    <- length(pa$obs_sv2)
  op   <- graphics::par(mar = c(4, 4.5, 2.5, 1)); on.exit(graphics::par(op))
  ylim <- range(c(pa$obs_sv2, pa$thresh)) * c(0.95, 1.05)
  plot(seq_len(K), pa$obs_sv2, type = "b", lwd = 2, col = "steelblue",
       pch = 16, ylim = ylim,
       main = main, xlab = "Component",
       ylab = expression(sigma^2))
  graphics::lines(seq_len(K), pa$thresh, type = "b", lwd = 2, col = "red",
                  lty = 2, pch = 1)
  graphics::abline(v = pa$n_retain + 0.5, lty = 3, col = "grey40")
  graphics::legend("topright", bty = "n",
                   legend = c("Observed", "Random 95th pct", "Retention cut"),
                   col    = c("steelblue", "red", "grey40"),
                   lty    = c(1, 2, 3), lwd = 2, pch = c(16, 1, NA))
  invisible(NULL)
}


#' Overlay a person's ipsatized profile against FAPA core profiles
#'
#' Plots the standardized ipsatized profile of person \code{p} alongside each
#' of the first \code{K} core profiles (also standardized), one panel per
#' dimension.
#'
#' @param bca    A list returned by \code{\link{fapa_bca}}.
#' @param Xtilde Numeric matrix (persons \eqn{\times} variables), ipsatized.
#' @param p      Integer. Row index of the focal person. Default \code{1}.
#' @param K      Integer. Number of core profiles to overlay. Default \code{2}.
#'
#' @return Invisibly returns \code{NULL}.
#' @export
plot_person_match <- function(bca, Xtilde, p = 1, K = 2) {
  op <- graphics::par(mfrow = c(1L, K), mar = c(4, 4.1, 2.5, 1))
  on.exit(graphics::par(op))
  y  <- scale(as.numeric(Xtilde[p, ]))
  for (i in seq_len(K)) {
    cp <- scale(bca$X0[, i])
    plot(y, type = "b", col = "steelblue",
         main  = sprintf("FAPA CP_%d vs Person #%d", i, p),
         xlab  = "Subscale", ylab = "Standardized")
    graphics::lines(cp, lty = 2, col = "tomato", lwd = 2)
    graphics::abline(h = 0, lty = 3, col = "grey70")
    graphics::legend("topright", bty = "n",
                     legend = c("Person", sprintf("CP%d", i)),
                     col = c("steelblue", "tomato"), lty = c(1, 2), lwd = 2)
  }
  invisible(NULL)
}


#' Distribution plots for Stage 2 principal angles
#'
#' Draws one histogram per dimension showing the bootstrap distribution of
#' principal angles, with the stability threshold marked as a vertical line.
#'
#' @param pr A list returned by \code{\link{fapa_procrustes}}.
#'
#' @return Invisibly returns \code{NULL}.
#' @export
plot_principal_angles <- function(pr) {
  K  <- pr$K
  op <- graphics::par(mfrow = c(1L, K), mar = c(4, 4.1, 2.5, 1))
  on.exit(graphics::par(op))
  for (k in seq_len(K)) {
    graphics::hist(pr$angles_mat[, k], breaks = 40,
                   main  = sprintf("Principal Angle \u2014 Dimension %d", k),
                   xlab  = "Degrees",
                   col   = "lightsteelblue", border = "white")
    graphics::abline(v = pr$angle_thresh, col = "red", lty = 2, lwd = 2)
    graphics::legend("topright", bty = "n",
                     legend = sprintf("< %.0f\u00b0 stability bound",
                                      pr$angle_thresh),
                     col = "red", lty = 2, lwd = 2)
  }
  invisible(NULL)
}


#' Distribution plots for Stage 3 Tucker congruence coefficients
#'
#' Draws one histogram per core profile showing the bootstrap distribution of
#' Tucker CCs, with reference lines at the fair (default 0.85) and high (0.95)
#' thresholds.
#'
#' @param tc        A list returned by \code{\link{fapa_tucker}}.
#' @param cc_thresh Numeric. Fair-similarity reference line. Default
#'   \code{0.85}.
#'
#' @return Invisibly returns \code{NULL}.
#' @export
plot_tucker_cc <- function(tc, cc_thresh = 0.85) {
  K  <- tc$K
  op <- graphics::par(mfrow = c(1L, K), mar = c(4, 4.1, 2.5, 1))
  on.exit(graphics::par(op))
  for (k in seq_len(K)) {
    graphics::hist(tc$cc_mat[, k], breaks = 40, xlim = c(-1, 1),
                   main  = sprintf("Tucker CC \u2014 Core Profile %d", k),
                   xlab  = "Congruence Coefficient",
                   col   = "lightsteelblue", border = "white")
    graphics::abline(v = cc_thresh, col = "red",    lty = 2, lwd = 2)
    graphics::abline(v = 0.95,      col = "orange", lty = 3, lwd = 2)
    graphics::legend("topleft", bty = "n",
                     legend = c(sprintf("%.2f (fair)", cc_thresh),
                                "0.95 (high)"),
                     col = c("red", "orange"), lty = c(2, 3), lwd = 2)
  }
  invisible(NULL)
}
