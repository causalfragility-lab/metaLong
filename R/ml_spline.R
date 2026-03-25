#' Spline-Based Nonlinear Time Trend in Longitudinal Meta-Analysis
#'
#' Fits a natural cubic spline meta-regression over follow-up time using the
#' pooled time-point estimates from [ml_meta()].  Produces a smooth pooled
#' trajectory with simultaneous pointwise confidence bands and tests for
#' nonlinearity.
#'
#' @details
#' The spline is fit by weighted least squares on the `ml_meta()` estimates,
#' using `1 / se^2` as weights (i.e., inverse squared SE weighting to reflect
#' the precision of each time-point estimate).  This is a second-stage model.
#'
#' For a fully joint spline model at the individual-effect level, users should
#' call [metafor::rma.mv()] directly with `mods = ~ ns(time, df)`.  This
#' function is primarily intended for visualisation and trajectory testing.
#'
#' @param meta_obj   Output from [ml_meta()].
#' @param df         Degrees of freedom for the natural cubic spline.  Default
#'   `3`.  A value of `1` recovers a linear fit.
#' @param n_pred     Number of prediction points for the smooth curve.  Default
#'   `200`.
#' @param alpha      Confidence level (inherits from `meta_obj` if `NULL`).
#' @param test_linear Logical.  If `TRUE`, performs an F-test of nonlinearity
#'   (spline df > 1 vs linear fit).  Default `TRUE`.
#'
#' @return Object of class `ml_spline` with elements:
#'   \describe{
#'     \item{`pred`}{`data.frame` with `time`, `fit`, `ci_lb`, `ci_ub` for
#'       the smooth prediction grid.}
#'     \item{`coef`}{Spline coefficient estimates.}
#'     \item{`vcov`}{Coefficient covariance matrix.}
#'     \item{`r_squared`}{Weighted R-squared of the spline fit.}
#'     \item{`p_nonlinear`}{p-value for nonlinearity test (if requested).}
#'     \item{`df`}{Spline degrees of freedom used.}
#'     \item{`meta_obj`}{The original `ml_meta` object (for plotting).}
#'   }
#'
#' @seealso [ml_meta()], [ml_plot()]
#'
#' @examples
#' dat  <- sim_longitudinal_meta(k = 10, times = c(0, 6, 12, 24), seed = 3)
#' meta <- ml_meta(dat, yi = "yi", vi = "vi", study = "study", time = "time")
#' spl  <- ml_spline(meta, df = 2)
#' print(spl)
#' plot(spl)
#'
#' @export
ml_spline <- function(meta_obj,
                      df          = 3L,
                      n_pred      = 200L,
                      alpha       = NULL,
                      test_linear = TRUE) {
  
  if (!inherits(meta_obj, "ml_meta"))
    stop("`meta_obj` must be output from ml_meta().", call. = FALSE)
  if (!requireNamespace("splines", quietly = TRUE))
    stop("Package 'splines' required.", call. = FALSE)
  
  alpha <- alpha %||% attr(meta_obj, "alpha") %||% 0.05
  
  # Remove time points with missing or degenerate estimates
  dat <- meta_obj[
    !is.na(meta_obj$theta) &
      !is.na(meta_obj$se)    &
      is.finite(meta_obj$se) &
      meta_obj$se > 0,
  ]
  
  if (nrow(dat) < (df + 1))
    stop("Too few estimable time points (", nrow(dat), ") for spline df = ", df,
         ".  Reduce `df` or add more time points.", call. = FALSE)
  
  t_obs <- dat$time
  y_obs <- dat$theta
  w_obs <- 1 / dat$se^2   # precision weights
  
  # ---- Build spline basis ------------------------------------------------
  # Knots placed at quantiles of observed time points
  knots <- stats::quantile(t_obs, probs = seq(0, 1, length.out = df + 1)[c(-1, -(df + 1))])
  basis <- splines::ns(t_obs, df = df, knots = if (df > 1) knots else NULL)
  
  X <- cbind(1, basis)
  
  # ---- Weighted least squares --------------------------------------------
  W      <- diag(w_obs)
  XtWX   <- t(X) %*% W %*% X
  XtWy   <- t(X) %*% W %*% y_obs
  
  # Robust solve
  beta   <- tryCatch(
    solve(XtWX, XtWy),
    error = function(e) .shim_ginv(XtWX) %*% XtWy
  )
  Vbeta  <- solve(XtWX)   # covariance of beta under WLS
  
  # ---- Prediction grid ---------------------------------------------------
  t_pred <- seq(min(t_obs), max(t_obs), length.out = n_pred)
  basis_pred <- splines::ns(t_pred, df = df,
                            knots  = if (df > 1) knots else NULL,
                            Boundary.knots = range(t_obs))
  X_pred <- cbind(1, basis_pred)
  
  fit_pred   <- as.vector(X_pred %*% beta)
  se_pred    <- sqrt(pmax(0, rowSums((X_pred %*% Vbeta) * X_pred)))
  crit_val   <- qt(1 - alpha / 2, df = max(nrow(dat) - ncol(X), 1))
  ci_lb_pred <- fit_pred - crit_val * se_pred
  ci_ub_pred <- fit_pred + crit_val * se_pred
  
  pred_df <- data.frame(
    time  = t_pred,
    fit   = fit_pred,
    ci_lb = ci_lb_pred,
    ci_ub = ci_ub_pred
  )
  
  # ---- Weighted R-squared -------------------------------------------------------
  y_hat   <- as.vector(X %*% beta)
  ss_tot  <- sum(w_obs * (y_obs - weighted.mean(y_obs, w_obs))^2)
  ss_res  <- sum(w_obs * (y_obs - y_hat)^2)
  r2      <- if (ss_tot > 0) 1 - ss_res / ss_tot else NA_real_
  
  # ---- Nonlinearity test (spline vs linear) ------------------------------
  p_nonlinear <- NA_real_
  if (test_linear && df > 1) {
    X_lin    <- cbind(1, t_obs)
    WLS_lin  <- tryCatch({
      b_lin  <- solve(t(X_lin) %*% W %*% X_lin, t(X_lin) %*% W %*% y_obs)
      y_lin  <- as.vector(X_lin %*% b_lin)
      ss_lin <- sum(w_obs * (y_obs - y_lin)^2)
      list(ss = ss_lin, df = nrow(dat) - 2)
    }, error = function(e) NULL)
    
    if (!is.null(WLS_lin)) {
      df_diff  <- ncol(X) - 2      # extra df for nonlinear terms
      df_resid <- nrow(dat) - ncol(X)
      if (df_diff > 0 && df_resid > 0) {
        F_stat      <- ((WLS_lin$ss - ss_res) / df_diff) / (ss_res / df_resid)
        p_nonlinear <- pf(F_stat, df1 = df_diff, df2 = df_resid, lower.tail = FALSE)
      }
    }
  }
  
  structure(
    list(
      pred        = pred_df,
      coef        = beta,
      vcov        = Vbeta,
      r_squared   = r2,
      p_nonlinear = p_nonlinear,
      df          = df,
      alpha       = alpha,
      meta_obj    = meta_obj
    ),
    class = "ml_spline"
  )
}


# ---- S3 methods -----------------------------------------------------------

#' @export
print.ml_spline <- function(x, digits = 3, ...) {
  cat("\n=== metaLong: Spline Time Trend ===\n")
  cat("Spline df:", x$df, "\n")
  cat("Weighted R-squared:", round(x$r_squared, digits), "\n")
  if (!is.na(x$p_nonlinear))
    cat("Nonlinearity test p-value:", round(x$p_nonlinear, digits), "\n")
  cat("\nPrediction range: time in [",
      round(min(x$pred$time), 2), ",",
      round(max(x$pred$time), 2), "]\n")
  invisible(x)
}

#' @export
plot.ml_spline <- function(x,
                           main  = "Spline Time Trend",
                           xlab  = "Time",
                           ylab  = "Pooled Effect",
                           col   = "#2166ac",
                           col_ci = NULL,
                           add   = FALSE,
                           ...) {
  meta <- x$meta_obj
  pred <- x$pred
  
  col_ci <- col_ci %||% {
    rgb_v <- col2rgb(col) / 255
    rgb(rgb_v[1], rgb_v[2], rgb_v[3], alpha = 0.2)
  }
  
  ok_meta <- !is.na(meta$theta)
  ylim    <- range(c(pred$ci_lb, pred$ci_ub,
                     meta$ci_lb[ok_meta], meta$ci_ub[ok_meta]), na.rm = TRUE)
  
  if (!add) {
    plot(pred$time, pred$fit, type = "n",
         xlim = range(pred$time), ylim = ylim,
         xlab = xlab, ylab = ylab, main = main, ...)
    abline(h = 0, lty = 2, col = "grey60")
  }
  
  # CI band
  polygon(
    c(pred$time, rev(pred$time)),
    c(pred$ci_lb, rev(pred$ci_ub)),
    col = col_ci, border = NA
  )
  
  # Spline curve
  lines(pred$time, pred$fit, col = col, lwd = 2)
  
  # Observed time-point estimates
  if (ok_meta[1]) {
    suppressWarnings(
      arrows(meta$time[ok_meta], meta$ci_lb[ok_meta],
             meta$time[ok_meta], meta$ci_ub[ok_meta],
             angle = 90, code = 3, length = 0.04, col = "grey40")
    )
    points(meta$time[ok_meta], meta$theta[ok_meta],
           pch = 19, col = "grey30", cex = 0.9)
  }
  
  nl_label <- if (!is.na(x$p_nonlinear))
    paste0("Nonlinearity p = ", round(x$p_nonlinear, 3))
  else
    paste0("Spline df = ", x$df)
  mtext(nl_label, side = 3, line = 0, adj = 1, cex = 0.75, col = "grey40")
  
  invisible(x)
}