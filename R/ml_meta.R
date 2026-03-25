#' Longitudinal Meta-Analysis with Robust Variance Estimation
#'
#' Fits a random-effects meta-analytic model at each unique time point in a
#' long-format dataset of multi-wave effect sizes.  Inference uses robust
#' variance estimation (RVE) with optional Tipton (2015) small-sample
#' corrections via the `clubSandwich` package.
#'
#' @section Engine choice:
#' Two fitting engines are supported:
#' \describe{
#'   \item{`"rma.uni"` (default)}{`metafor::rma.uni()` -- appropriate when each
#'     study contributes exactly one effect size per time point.  Simpler,
#'     faster, and stores `tau2` directly from the REML estimate.}
#'   \item{`"rma.mv"`}{`metafor::rma.mv()` with a prebuilt working covariance
#'     matrix -- appropriate when studies contribute *multiple* effect sizes at
#'     the same time point (dependent effects within cluster).  Requires the
#'     `rho` argument.}
#' }
#'
#' @param data A `data.frame` in **long format**: one row per study x time point.
#' @param yi   Character. Name of the effect-size column.
#' @param vi   Character. Name of the sampling-variance column.
#' @param study Character. Name of the study-ID column (cluster variable).
#' @param time  Character. Name of the follow-up time column (numeric).
#' @param alpha Significance level for confidence intervals and p-values.
#'   Default `0.05`.
#' @param rho   Assumed within-study correlation between effect sizes (used only
#'   when `engine = "rma.mv"`).  Default `0.8`.
#' @param small_sample Logical. If `TRUE` (default), applies CR2 sandwich
#'   variance estimation with Satterthwaite degrees of freedom (Tipton, 2015).
#'   If `FALSE`, uses uncorrected z-based inference.
#' @param min_k Integer. Minimum number of studies required to fit a model at
#'   a given time point.  Default `2`.
#' @param method Character. Variance estimator passed to metafor.  Default
#'   `"REML"`.
#' @param engine Character. Fitting engine: `"rma.uni"` (default) or
#'   `"rma.mv"`.  See section *Engine choice*.
#'
#' @return An object of class `ml_meta` (a `data.frame`) with one row per time
#'   point and columns: `time`, `k`, `theta`, `se`, `df`, `t_stat`, `p_val`,
#'   `ci_lb`, `ci_ub`, `tau2`, `note`.
#'
#'   Attributes:
#'   \describe{
#'     \item{`"fits"`}{Named list of fitted model objects (one per time point).}
#'     \item{`"weights_by_time"`}{Named list of weight vectors for downstream
#'       use by [ml_sens()] and [ml_benchmark()].}
#'     \item{`"engine"`, `"alpha"`, `"rho"`, `"small_sample"`}{Call metadata.}
#'   }
#'
#' @references
#' Hedges, L. V., Tipton, E., & Johnson, M. C. (2010). Robust variance
#' estimation in meta-regression with dependent effect size estimates.
#' *Research Synthesis Methods*, 1(1), 39-65.
#'
#' Tipton, E. (2015). Small sample adjustments for robust variance estimation
#' with meta-regression. *Psychological Methods*, 20(3), 375-393.
#'
#' @seealso [ml_sens()], [ml_benchmark()], [ml_spline()]
#'
#' @examples
#' dat <- sim_longitudinal_meta(k = 10, times = c(0, 6, 12), seed = 1)
#' result <- ml_meta(dat, yi = "yi", vi = "vi", study = "study", time = "time")
#' print(result)
#' plot(result)
#'
#' # rma.mv engine for dependent effects
#' \donttest{
#' result_mv <- ml_meta(dat, yi = "yi", vi = "vi", study = "study", time = "time",
#'                      engine = "rma.mv", rho = 0.8)
#' }
#'
#' @export
ml_meta <- function(data,
                    yi,
                    vi,
                    study,
                    time,
                    alpha        = 0.05,
                    rho          = 0.8,
                    small_sample = TRUE,
                    min_k        = 2L,
                    method       = "REML",
                    engine       = c("rma.uni", "rma.mv")) {

  engine <- match.arg(engine)

  # ---- Input validation --------------------------------------------------
  if (!is.data.frame(data))
    stop("`data` must be a data.frame.", call. = FALSE)
  if (!requireNamespace("metafor", quietly = TRUE))
    stop("Package 'metafor' required. Install with install.packages('metafor').",
         call. = FALSE)

  dat       <- .prep_data(data, yi, vi, study, time)
  time_vals <- sort(unique(dat$.time))

  out_list     <- vector("list", length(time_vals))
  weights_list <- vector("list", length(time_vals))
  fits_list    <- vector("list", length(time_vals))
  names(weights_list) <- as.character(time_vals)
  names(fits_list)    <- as.character(time_vals)

  # ---- Main loop ---------------------------------------------------------
  for (j in seq_along(time_vals)) {
    tt  <- time_vals[j]
    dtt <- dat[dat$.time == tt, ]
    k_t <- length(unique(dtt$.study))

    skeleton <- data.frame(
      time    = tt,
      k       = k_t,
      theta   = NA_real_,
      se      = NA_real_,
      df      = NA_real_,
      t_stat  = NA_real_,
      p_val   = NA_real_,
      ci_lb   = NA_real_,
      ci_ub   = NA_real_,
      tau2    = NA_real_,
      note    = NA_character_,
      stringsAsFactors = FALSE
    )

    if (k_t < min_k) {
      skeleton$note <- paste0("k < min_k (", k_t, " < ", min_k, ")")
      out_list[[j]] <- skeleton
      next
    }

    res <- if (engine == "rma.uni") {
      .fit_intercept_uni(
        yi           = dtt$.yi,
        vi           = dtt$.vi,
        cluster      = dtt$.study,
        small_sample = small_sample,
        alpha        = alpha,
        method       = method
      )
    } else {
      .fit_intercept_mv(
        yi           = dtt$.yi,
        vi           = dtt$.vi,
        cluster      = dtt$.study,
        rho          = rho,
        small_sample = small_sample,
        alpha        = alpha,
        method       = method
      )
    }

    if (is.null(res)) {
      skeleton$note <- "model_failed"
      out_list[[j]] <- skeleton
      next
    }

    # Store tau2-adjusted weights -- used by ml_sens() for sy2 computation
    tau2_t <- if (!is.na(res$tau2) && res$tau2 > 0) res$tau2 else 0
    weights_list[[as.character(tt)]] <- 1 / (dtt$.vi + tau2_t)
    fits_list[[as.character(tt)]]    <- res$fit

    out_list[[j]] <- data.frame(
      time    = tt,
      k       = k_t,
      theta   = res$theta,
      se      = res$se,
      df      = res$df,
      t_stat  = res$theta / res$se,
      p_val   = res$p_val,
      ci_lb   = res$ci_lb,
      ci_ub   = res$ci_ub,
      tau2    = res$tau2,
      note    = NA_character_,
      stringsAsFactors = FALSE
    )
  }

  out           <- do.call(rbind, out_list)
  rownames(out) <- NULL

  attr(out, "fits")            <- fits_list
  attr(out, "weights_by_time") <- weights_list
  attr(out, "call")            <- match.call()
  attr(out, "alpha")           <- alpha
  attr(out, "rho")             <- rho
  attr(out, "small_sample")    <- small_sample
  attr(out, "engine")          <- engine

  class(out) <- c("ml_meta", "data.frame")
  out
}


# ---- S3 methods -----------------------------------------------------------

#' @export
print.ml_meta <- function(x, digits = 3, ...) {
  cat("\n=== metaLong: Longitudinal Pooled Effects ===\n")
  cat("Engine:", attr(x, "engine") %||% "rma.uni", "\n")
  cat("Time points:", length(unique(x$time)), "\n")
  cat("Small-sample correction:", isTRUE(attr(x, "small_sample")), "\n\n")

  pr       <- as.data.frame(x)
  num_cols <- intersect(c("theta", "se", "df", "t_stat", "p_val",
                           "ci_lb", "ci_ub", "tau2"), names(pr))
  pr[num_cols] <- lapply(pr[num_cols], round, digits = digits)
  show_cols <- intersect(c("time", "k", "theta", "se", "df",
                            "p_val", "ci_lb", "ci_ub", "tau2"), names(pr))
  print(pr[, show_cols], row.names = FALSE)
  invisible(x)
}

#' @export
summary.ml_meta <- function(object, ...) {
  cat("\n=== ml_meta summary ===\n")
  cat("Engine:          ", attr(object, "engine") %||% "rma.uni", "\n")
  cat("Time points:     ", nrow(object), "\n")
  alpha <- attr(object, "alpha") %||% 0.05
  sig   <- !is.na(object$p_val) & object$p_val < alpha
  cat("Significant:     ", sum(sig, na.rm = TRUE), "/",
      sum(!is.na(object$p_val)), "\n")
  cat("Effect range:    [", round(min(object$theta, na.rm = TRUE), 3), ",",
      round(max(object$theta, na.rm = TRUE), 3), "]\n")
  if ("tau2" %in% names(object))
    cat("tau2 range:      [", round(min(object$tau2,  na.rm = TRUE), 3), ",",
        round(max(object$tau2,  na.rm = TRUE), 3), "]\n")
  invisible(object)
}

#' @export
plot.ml_meta <- function(x, xlab = "Time", ylab = "Pooled Effect",
                          main = "Pooled Effects Over Time",
                          pch = 19, col = "#2166ac", ...) {
  ok <- !is.na(x$theta) & !is.na(x$ci_lb)
  if (!any(ok)) { message("No estimable time points to plot."); return(invisible(x)) }

  ylim <- range(c(x$ci_lb[ok], x$ci_ub[ok]), na.rm = TRUE)
  ylim <- ylim + diff(ylim) * c(-0.05, 0.05)

  plot(x$time[ok], x$theta[ok], type = "b", pch = pch, col = col,
       ylim = ylim, xlab = xlab, ylab = ylab, main = main, ...)

  polygon(
    c(x$time[ok], rev(x$time[ok])),
    c(x$ci_lb[ok], rev(x$ci_ub[ok])),
    col    = rgb(col2rgb(col)[1]/255, col2rgb(col)[2]/255,
                 col2rgb(col)[3]/255, 0.15),
    border = NA
  )

  abline(h = 0, lty = 2, col = "grey60")

  sig <- ok & !is.na(x$p_val) & x$p_val < (attr(x, "alpha") %||% 0.05)
  if (any(sig))
    points(x$time[sig], x$theta[sig], pch = 8, col = "#d73027", cex = 0.9)

  invisible(x)
}

#' Extract stored fitted model objects
#'
#' @param x An `ml_meta` object.
#' @return Named list of fitted model objects, one per estimable time point.
#' @export
fits <- function(x) {
  if (!inherits(x, "ml_meta"))
    stop("`x` must be an ml_meta object.", call. = FALSE)
  attr(x, "fits")
}

# Null coalescing operator (R < 4.4 compatibility)
`%||%` <- function(a, b) if (!is.null(a)) a else b

#' Tidy an metaLong object into a clean data frame
#'
#' @param x A `ml_sens` or `ml_benchmark` object.
#' @param ... Additional arguments (unused).
#' @return A tidy `data.frame`.
#' @export
tidy <- function(x, ...) UseMethod("tidy")
