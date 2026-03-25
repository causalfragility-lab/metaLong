#' Combined Publication-Ready Trajectory Figure
#'
#' Produces a multi-panel figure combining the pooled trajectory, confidence
#' band, spline fit (if supplied), and ITCV sensitivity profile.  Designed for
#' direct inclusion in manuscripts.
#'
#' @param meta_obj  Output from [ml_meta()] (required).
#' @param sens_obj  Output from [ml_sens()] (optional; adds ITCV panel).
#' @param bench_obj Output from [ml_benchmark()] (optional; adds benchmark marks).
#' @param spline_obj Output from [ml_spline()] (optional; overlays spline).
#' @param frag_obj  Output from [ml_fragility()] (optional; adds fragility panel).
#' @param ncol      Number of columns in the panel layout.  Default auto.
#' @param main      Overall figure title.
#' @param col_effect Colour for the pooled effect trajectory.
#' @param col_sens   Colour for the ITCV line.
#' @param col_spline Colour for the spline curve.
#' @param delta     Fragility benchmark line on the ITCV panel.  Inherits from
#'   `sens_obj` if available.
#'
#' @return Invisibly returns a list of the objects passed in.
#'
#' @examples
#' dat  <- sim_longitudinal_meta(k = 10, times = c(0, 6, 12), seed = 1)
#' meta <- ml_meta(dat, yi = "yi", vi = "vi", study = "study", time = "time")
#' sens <- ml_sens(dat, meta, yi = "yi", vi = "vi", study = "study", time = "time")
#' ml_plot(meta, sens_obj = sens)
#'
#' \donttest{
#' spl  <- ml_spline(meta, df = 2)
#' ml_plot(meta, sens_obj = sens, spline_obj = spl,
#'         main = "Longitudinal Meta-Analysis Profile")
#' }
#'
#' @export
ml_plot <- function(meta_obj,
                     sens_obj   = NULL,
                     bench_obj  = NULL,
                     spline_obj = NULL,
                     frag_obj   = NULL,
                     ncol       = NULL,
                     main       = NULL,
                     col_effect = "#2166ac",
                     col_sens   = "#d73027",
                     col_spline = "#1a9641",
                     delta      = NULL) {

  if (!inherits(meta_obj, "ml_meta"))
    stop("`meta_obj` must be output from ml_meta().", call. = FALSE)

  # Determine panels
  panels <- "effect"
  if (!is.null(sens_obj))   panels <- c(panels, "sens")
  if (!is.null(frag_obj))   panels <- c(panels, "frag")

  n_panels <- length(panels)
  ncol     <- ncol %||% min(n_panels, 2L)
  nrow     <- ceiling(n_panels / ncol)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(nrow, ncol), mar = c(4.5, 4.5, 3, 1.5), oma = c(0, 0, 3, 0))

  # ---- Panel 1: Pooled trajectory ----------------------------------------
  ok    <- !is.na(meta_obj$theta)
  ylim  <- .safe_range(c(meta_obj$ci_lb[ok], meta_obj$ci_ub[ok]))
  if (!is.null(spline_obj))
    ylim <- .safe_range(c(ylim, spline_obj$pred$ci_lb, spline_obj$pred$ci_ub))

  plot(meta_obj$time[ok], meta_obj$theta[ok],
       type = "n", ylim = ylim,
       xlab = "Follow-up time", ylab = "Pooled effect",
       main = "Pooled Trajectory")
  abline(h = 0, lty = 2, col = "grey70")

  # CI band
  .filled_band(meta_obj$time[ok], meta_obj$ci_lb[ok], meta_obj$ci_ub[ok],
               col_effect, alpha = 0.15)

  # Spline overlay
  if (!is.null(spline_obj)) {
    pr <- spline_obj$pred
    .filled_band(pr$time, pr$ci_lb, pr$ci_ub, col_spline, alpha = 0.12)
    lines(pr$time, pr$fit, col = col_spline, lwd = 2.2, lty = 1)
    nl_p <- spline_obj$p_nonlinear
    if (!is.na(nl_p))
      mtext(paste0("Nonlinearity p=", round(nl_p, 3)),
            side = 3, adj = 1, cex = 0.70, col = col_spline)
  }

  # Observed points + error bars
  suppressWarnings(
    arrows(meta_obj$time[ok], meta_obj$ci_lb[ok],
           meta_obj$time[ok], meta_obj$ci_ub[ok],
           angle = 90, code = 3, length = 0.04, col = "grey50", lwd = 0.8)
  )

  sig <- ok & !is.na(meta_obj$p_val) & meta_obj$p_val < (attr(meta_obj, "alpha") %||% 0.05)
  points(meta_obj$time[ok & !sig], meta_obj$theta[ok & !sig],
         pch = 21, bg = "white", col = col_effect, cex = 1.1)
  points(meta_obj$time[ok & sig],  meta_obj$theta[ok & sig],
         pch = 19, col = col_effect, cex = 1.1)

  # Benchmark exceedance marks
  if (!is.null(bench_obj)) {
    fs <- attr(bench_obj, "fragile_summary")
    if (!is.null(fs) && any(fs$any_beats, na.rm = TRUE)) {
      beat_times <- fs$time[!is.na(fs$any_beats) & fs$any_beats]
      beat_y     <- meta_obj$theta[meta_obj$time %in% beat_times]
      points(beat_times, beat_y, pch = 4, col = "#d73027", cex = 1.4, lwd = 2)
    }
  }

  # ---- Panel 2: ITCV sensitivity -----------------------------------------
  if ("sens" %in% panels && !is.null(sens_obj)) {
    delta <- delta %||% attr(sens_obj, "delta") %||% 0.15
    sok   <- !is.na(sens_obj$itcv_alpha)
    ylim2 <- c(0, max(c(sens_obj$itcv_alpha[sok], delta), na.rm = TRUE) * 1.2)
    col_pts <- ifelse(sens_obj$itcv_alpha[sok] < delta, "#d73027", col_sens)

    plot(sens_obj$time[sok], sens_obj$itcv_alpha[sok],
         type = "n", ylim = ylim2,
         xlab = "Follow-up time",
         ylab = "ITCV_alpha(t)",
         main = "Sensitivity Profile")

    lines(sens_obj$time[sok], sens_obj$itcv_alpha[sok],
          col = "grey60", lwd = 1.2)
    points(sens_obj$time[sok], sens_obj$itcv_alpha[sok],
           pch = 19, col = col_pts, cex = 1.2)
    abline(h = delta, lty = 2, col = "grey30")
    text(max(sens_obj$time[sok]), delta,
         labels = paste0("delta=", delta),
         adj = c(1.1, -0.4), cex = 0.75, col = "grey30")
  }

  # ---- Panel 3: Fragility ------------------------------------------------
  if ("frag" %in% panels && !is.null(frag_obj)) {
    fok  <- !is.na(frag_obj$fragility_index)

    if (!any(fok)) {
      # Nothing estimable -- show informative empty panel
      plot.new()
      title(main = "Fragility by Time Point",
            xlab = "Follow-up time", ylab = "Fragility index")
      text(0.5, 0.5, "No fragility indices estimable\n(increase max_k or study count)",
           adj = c(0.5, 0.5), col = "grey50")
    } else {
      cols <- ifelse(frag_obj$sig_original[fok], col_effect, "#d73027")
      barplot(
        frag_obj$fragility_index[fok],
        names.arg = frag_obj$time[fok],
        col       = cols,
        xlab      = "Follow-up time",
        ylab      = "Fragility index",
        main      = "Fragility by Time Point"
      )
      abline(h = 1, lty = 2, col = "grey50")
    }
  }

  # Overall title
  if (!is.null(main))
    mtext(main, outer = TRUE, cex = 1.1, font = 2)

  invisible(list(meta = meta_obj, sens = sens_obj,
                 bench = bench_obj, spline = spline_obj, frag = frag_obj))
}


# ---- Internal plotting helpers -------------------------------------------

.safe_range <- function(x) {
  r <- range(x, na.rm = TRUE)
  if (diff(r) == 0) r + c(-0.1, 0.1) else r * c(ifelse(r[1] < 0, 1.1, 0.9),
                                                   ifelse(r[2] > 0, 1.1, 0.9))
}

.filled_band <- function(x, lb, ub, col, alpha = 0.15) {
  rgb_v <- col2rgb(col) / 255
  fill  <- rgb(rgb_v[1], rgb_v[2], rgb_v[3], alpha = alpha)
  polygon(c(x, rev(x)), c(lb, rev(ub)), col = fill, border = NA)
}
