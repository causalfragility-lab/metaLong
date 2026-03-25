#' Time-Varying Sensitivity Analysis via Longitudinal ITCV
#'
#' Computes the Impact Threshold for a Confounding Variable (ITCV) at each
#' follow-up time point using the pooled estimates and robust inference from
#' [ml_meta()].  Two versions are returned: the raw ITCV (threshold to nullify
#' the pooled effect) and the significance-adjusted ITCV_alpha (threshold to render
#' the result non-significant under small-sample-corrected inference).
#'
#' @section Mathematical background:
#' At each time \eqn{t}, let \eqn{\hat\theta_t} be the pooled effect,
#' \eqn{s_{y,t}^2} the weighted variance of observed effect sizes, and
#' \eqn{c_t = t_{1-\alpha/2,\nu_t} \cdot \widehat{SE}(\hat\theta_t)} the
#' minimum effect still deemed significant.  The correlation-scale pooled effect
#' is
#' \deqn{r_t = \hat\theta_t / \sqrt{\hat\theta_t^2 + s_{y,t}^2}}
#' and the raw ITCV is \eqn{\sqrt{|r_t|}}.  The significance-adjusted version
#' replaces \eqn{\hat\theta_t} with \eqn{|\hat\theta_t| - c_t}.
#'
#' @param data      A `data.frame` in long format (same as passed to [ml_meta()]).
#' @param meta_obj  Output from [ml_meta()].
#' @param yi,vi,study,time  Column names (same meaning as in [ml_meta()]).
#' @param alpha     Significance level.  Defaults to the value stored in
#'   `meta_obj` (or `0.05` if absent).
#' @param delta     Numeric. User-defined practical fragility benchmark: time
#'   points with `ITCV_alpha(t) < delta` are flagged as "practically fragile".
#'   Default `0.15`.
#'
#' @return An object of class `ml_sens` (a `data.frame`) with columns:
#'   \describe{
#'     \item{`time`}{Follow-up time.}
#'     \item{`theta`, `se`, `df`}{Copied from `meta_obj`.}
#'     \item{`sy`}{Weighted SD of observed effect sizes.}
#'     \item{`r_effect`}{Pooled effect on correlation scale.}
#'     \item{`itcv`}{Raw ITCV: confounding needed to nullify the estimate.}
#'     \item{`itcv_alpha`}{Significance-adjusted ITCV: confounding needed to
#'       make the result non-significant.}
#'     \item{`fragile`}{Logical; `TRUE` when `itcv_alpha < delta`.}
#'   }
#'   Attributes include trajectory summaries (`itcv_min`, `itcv_mean`,
#'   `fragile_prop`) and a `"fragile_times"` character vector.
#'
#' @references
#' Frank, K. A. (2000). Impact of a confounding variable on a regression
#' coefficient. *Sociological Methods & Research*, 29(2), 147-194.
#'
#' @seealso [ml_meta()], [ml_benchmark()], [ml_plot()]
#'
#' @examples
#' dat  <- sim_longitudinal_meta(k = 10, times = c(0, 6, 12), seed = 1)
#' meta <- ml_meta(dat, yi = "yi", vi = "vi", study = "study", time = "time")
#' sens <- ml_sens(dat, meta, yi = "yi", vi = "vi", study = "study", time = "time")
#' print(sens)
#' plot(sens)
#'
#' @export
ml_sens <- function(data,
                    meta_obj,
                    yi,
                    vi,
                    study,
                    time,
                    alpha = NULL,
                    delta = 0.15) {

  # ---- Validation --------------------------------------------------------
  if (!inherits(meta_obj, "ml_meta"))
    stop("`meta_obj` must be output from ml_meta().", call. = FALSE)

  alpha <- alpha %||% attr(meta_obj, "alpha") %||% 0.05

  dat          <- .prep_data(data, yi, vi, study, time)
  weights_list <- attr(meta_obj, "weights_by_time")
  time_vals    <- meta_obj$time

  out_list <- vector("list", length(time_vals))

  for (j in seq_along(time_vals)) {
    tt      <- time_vals[j]
    dtt     <- dat[dat$.time == tt, ]
    res_tt  <- meta_obj[meta_obj$time == tt, ]

    theta_hat <- res_tt$theta
    se_hat    <- res_tt$se
    df_hat    <- res_tt$df
    # tau2 from the fitted model -- used for random-effects weight alignment
    tau2_hat  <- if ("tau2" %in% names(res_tt) && !is.na(res_tt$tau2))
      as.numeric(res_tt$tau2) else 0

    skeleton <- data.frame(
      time       = tt,
      theta      = theta_hat,
      se         = se_hat,
      df         = df_hat,
      sy         = NA_real_,
      r_effect   = NA_real_,
      itcv       = NA_real_,
      itcv_alpha = NA_real_,
      fragile    = NA,
      stringsAsFactors = FALSE
    )

    if (is.na(theta_hat) || is.na(se_hat) || is.na(df_hat)) {
      out_list[[j]] <- skeleton
      next
    }

    # Weights: prefer stored tau2-adjusted weights from ml_meta(),
    # falling back to tau2-adjusted weights computed here, then plain 1/vi
    w_key <- as.character(tt)
    w <- if (!is.null(weights_list) && !is.null(weights_list[[w_key]])) {
      weights_list[[w_key]]                      # already tau2-adjusted
    } else {
      1 / (dtt$.vi + tau2_hat)                   # compute here
    }

    # Weighted dispersion of observed effects
    yi_vals <- dtt$.yi
    sy2     <- sum(w * (yi_vals - theta_hat)^2, na.rm = TRUE) / sum(w, na.rm = TRUE)
    sy      <- sqrt(sy2)

    # Correlation-scale pooled effect
    r_effect <- .to_r_scale(theta_hat, sy2)

    # Raw ITCV
    itcv <- .itcv_from_r(r_effect)

    # Significance-adjusted ITCV
    crit         <- qt(1 - alpha / 2, df = df_hat)
    theta_star   <- abs(theta_hat) - crit * se_hat

    itcv_alpha <- if (theta_star <= 0) {
      0
    } else {
      r_star <- .to_r_scale(theta_star, sy2)
      .itcv_from_r(r_star)
    }

    out_list[[j]] <- data.frame(
      time       = tt,
      theta      = theta_hat,
      se         = se_hat,
      df         = df_hat,
      sy         = sy,
      r_effect   = r_effect,
      itcv       = itcv,
      itcv_alpha = itcv_alpha,
      fragile    = itcv_alpha < delta,
      stringsAsFactors = FALSE
    )
  }

  out           <- do.call(rbind, out_list)
  rownames(out) <- NULL

  # ---- Trajectory summaries ----------------------------------------------
  valid_itcv    <- out$itcv_alpha[!is.na(out$itcv_alpha)]
  fragile_times <- out$time[!is.na(out$fragile) & out$fragile]

  attr(out, "itcv_min")      <- if (length(valid_itcv) > 0) min(valid_itcv)  else NA_real_
  attr(out, "itcv_mean")     <- if (length(valid_itcv) > 0) mean(valid_itcv) else NA_real_
  attr(out, "fragile_prop")  <- if (length(valid_itcv) > 0) mean(out$itcv_alpha < delta, na.rm = TRUE) else NA_real_
  attr(out, "fragile_times") <- fragile_times
  attr(out, "delta")         <- delta
  attr(out, "alpha")         <- alpha

  class(out) <- c("ml_sens", "data.frame")
  out
}


# ---- S3 methods -----------------------------------------------------------

#' @export
print.ml_sens <- function(x, digits = 3, ...) {
  cat("\n=== metaLong: Longitudinal Sensitivity (ITCV) ===\n")
  cat("delta (fragility benchmark):", attr(x, "delta"), "\n")
  cat("ITCV_alpha range: [",
      round(attr(x, "itcv_min"), digits), ",",
      round(max(x$itcv_alpha, na.rm = TRUE), digits), "]\n")
  cat("Fragile proportion:", round(attr(x, "fragile_prop"), digits), "\n\n")

  pr       <- as.data.frame(x)
  num_cols <- c("theta", "se", "sy", "r_effect", "itcv", "itcv_alpha")
  pr[num_cols] <- lapply(pr[num_cols], round, digits = digits)
  print(pr[, c("time", "theta", "sy", "itcv", "itcv_alpha", "fragile")],
        row.names = FALSE)
  invisible(x)
}

#' @export
summary.ml_sens <- function(object, ...) {
  cat("ITCV_alpha min (weakest point):", round(attr(object, "itcv_min"), 3), "\n")
  cat("ITCV_alpha mean:               ", round(attr(object, "itcv_mean"), 3), "\n")
  cat("Fragile time points:           ", paste(attr(object, "fragile_times"), collapse = ", "), "\n")
  invisible(object)
}

#' @export
plot.ml_sens <- function(x, delta = NULL,
                          col_robust   = "#2166ac",
                          col_fragile  = "#d73027",
                          main = "Longitudinal Sensitivity Profile",
                          xlab = "Time", ...) {
  delta <- delta %||% attr(x, "delta") %||% 0.15
  ok    <- !is.na(x$itcv_alpha)
  if (!any(ok)) { message("No estimable ITCV values to plot."); return(invisible(x)) }

  ylim <- c(0, max(x$itcv_alpha[ok], delta, na.rm = TRUE) * 1.15)

  col_pts <- ifelse(x$itcv_alpha[ok] < delta, col_fragile, col_robust)

  plot(x$time[ok], x$itcv_alpha[ok], type = "n",
       ylim = ylim, xlab = xlab, ylab = "ITCV_alpha(t)",
       main = main, ...)
  lines(x$time[ok], x$itcv_alpha[ok], col = "grey50", lwd = 1.2)
  points(x$time[ok], x$itcv_alpha[ok], pch = 19, col = col_pts, cex = 1.2)
  abline(h = delta, lty = 2, col = "grey30")
  text(max(x$time[ok]), delta, labels = paste0("delta = ", delta),
       adj = c(1.05, -0.4), cex = 0.75)

  invisible(x)
}

#' @export
tidy.ml_sens <- function(x, ...) {
  out <- as.data.frame(x)
  out$interpretation <- ifelse(
    is.na(out$fragile), "insufficient data",
    ifelse(out$fragile,
           "fragile: small confounding could overturn result",
           "robust: large confounding needed to overturn result")
  )
  out
}
