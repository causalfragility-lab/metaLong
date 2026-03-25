#' metaLong: Longitudinal Meta-Analysis with Robust Inference and Sensitivity Analysis
#'
#' @description
#' `metaLong` provides a coherent workflow for synthesising evidence from studies
#' that report outcomes at multiple follow-up time points.  The package covers:
#'
#' * **Pooling**: `ml_meta()` fits a random-effects or fixed-effects model at each
#'   time point using robust variance estimation (RVE) and optional Tipton
#'   small-sample corrections via `clubSandwich`.
#' * **Sensitivity**: `ml_sens()` computes the time-varying Impact Threshold for a
#'   Confounding Variable (ITCV), both raw and significance-adjusted.
#' * **Benchmark calibration**: `ml_benchmark()` regresses each observed study-level
#'   covariate and compares its partial correlation against the ITCV threshold.
#' * **Nonlinear time trends**: `ml_spline()` fits natural-cubic-spline meta-regression
#'   over follow-up time with pointwise confidence bands.
#' * **Fragility**: `ml_fragility()` computes leave-one-out and leave-k-out fragility
#'   indices across the trajectory.
#' * **Visualisation**: `ml_plot()` and S3 `plot()` methods produce publication-ready
#'   trajectory, sensitivity, and benchmark figures.
#'
#' @section Typical workflow:
#' ```r
#' meta  <- ml_meta(data, yi="g", vi="vg", study="id", time="wave")
#' sens  <- ml_sens(data, meta, yi="g", vi="vg", study="id", time="wave")
#' bench <- ml_benchmark(data, meta, sens, yi="g", vi="vg",
#'                        study="id", time="wave",
#'                        covariates=c("pub_year","n","quality"))
#' spl   <- ml_spline(meta, df=3)
#' ml_plot(meta, sens, bench, spl)
#' ```
#'
#' @references
#' Frank, K. A. (2000). Impact of a confounding variable on a regression coefficient.
#' *Sociological Methods & Research*, 29(2), 147-194.
#'
#' Hedges, L. V., Tipton, E., & Johnson, M. C. (2010). Robust variance estimation
#' in meta-regression with dependent effect size estimates.
#' *Research Synthesis Methods*, 1(1), 39-65.
#'
#' Tipton, E. (2015). Small sample adjustments for robust variance estimation with
#' meta-regression. *Psychological Methods*, 20(3), 375-393.
#'
#' @keywords internal
"_PACKAGE"


# ============================================================
#  Internal helpers
# ============================================================

#' Build working covariance matrix (uses clubSandwich if available, else shim)
#'
#' @param vi  numeric vector of sampling variances
#' @param cluster factor/character of study IDs
#' @param rho  within-study working correlation
#' @keywords internal
.build_V <- function(vi, cluster, rho) {
  .cs_impute_cov(vi = vi, cluster = cluster, r = rho)
}

#' rma.uni engine: fits intercept-only model, returns tidy inference
#'
#' Used when each study has at most one effect per time point.
#' Stores tau2 directly from the REML estimate.
#' @keywords internal
.fit_intercept_uni <- function(yi, vi, cluster, small_sample, alpha,
                                method = "REML") {
  fit <- tryCatch(
    suppressWarnings(
      metafor::rma.uni(yi = yi, vi = vi, method = method)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  theta_hat <- as.numeric(stats::coef(fit)[1])
  tau2_hat  <- if (!is.null(fit$tau2)) as.numeric(fit$tau2) else NA_real_

  if (small_sample) {
    cs <- tryCatch(
      .cs_coef_test(fit, vcov = "CR2", cluster = cluster,
                    test = "Satterthwaite"),
      error = function(e) NULL
    )
    if (is.null(cs)) return(NULL)
    se_hat <- as.numeric(cs$SE[1])
    df_hat <- as.numeric(cs$df[1])
    p_val  <- as.numeric(cs$p_Satt[1])
  } else {
    se_hat <- as.numeric(fit$se)
    df_hat <- Inf
    p_val  <- as.numeric(fit$pval)
  }

  crit <- qt(1 - alpha / 2, df = df_hat)
  list(
    theta = theta_hat,
    se    = se_hat,
    df    = df_hat,
    p_val = p_val,
    ci_lb = theta_hat - crit * se_hat,
    ci_ub = theta_hat + crit * se_hat,
    tau2  = tau2_hat,
    fit   = fit
  )
}

#' rma.mv engine: builds working V matrix, fits rma.mv, returns tidy inference
#'
#' Used when studies contribute multiple effects at the same time point.
#' @keywords internal
.fit_intercept_mv <- function(yi, vi, cluster, rho, small_sample, alpha,
                               method = "REML") {
  V <- tryCatch(
    .build_V(vi, cluster, rho),
    error = function(e) NULL
  )
  if (is.null(V)) return(NULL)

  fit <- tryCatch(
    suppressWarnings(
      metafor::rma.mv(yi = yi, V = V, method = method, sparse = TRUE)
    ),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  theta_hat <- as.numeric(fit$b)
  # rma.mv stores sigma2 not tau2; use first sigma2 if available
  tau2_hat  <- if (!is.null(fit$sigma2) && length(fit$sigma2) > 0)
    as.numeric(fit$sigma2[1]) else NA_real_

  if (small_sample) {
    cs <- tryCatch(
      .cs_coef_test(fit, vcov = "CR2", cluster = cluster,
                    test = "Satterthwaite"),
      error = function(e) NULL
    )
    if (is.null(cs)) return(NULL)
    se_hat <- as.numeric(cs$SE[1])
    df_hat <- as.numeric(cs$df[1])
    p_val  <- as.numeric(cs$p_Satt[1])
  } else {
    se_hat <- as.numeric(fit$se)
    df_hat <- Inf
    p_val  <- as.numeric(fit$pval)
  }

  crit <- qt(1 - alpha / 2, df = df_hat)
  list(
    theta = theta_hat,
    se    = se_hat,
    df    = df_hat,
    p_val = p_val,
    ci_lb = theta_hat - crit * se_hat,
    ci_ub = theta_hat + crit * se_hat,
    tau2  = tau2_hat,
    fit   = fit
  )
}

# Keep legacy alias so any internal calls still work
.fit_intercept <- .fit_intercept_mv

#' Weighted mean and dispersion of y
#' @keywords internal
.weighted_stats <- function(yi, vi) {
  w      <- 1 / vi
  theta  <- sum(w * yi, na.rm = TRUE) / sum(w, na.rm = TRUE)
  sy2    <- sum(w * (yi - theta)^2, na.rm = TRUE) / sum(w, na.rm = TRUE)
  list(theta = theta, sy2 = sy2, w = w)
}

#' Convert effect to correlation scale (for ITCV)
#' @keywords internal
.to_r_scale <- function(theta, sy2) {
  theta / sqrt(theta^2 + sy2)
}

#' Compute ITCV from r-scale value
#' @keywords internal
.itcv_from_r <- function(r) sqrt(abs(r))

#' Validate required columns exist in data
#' @keywords internal
.check_cols <- function(data, ...) {
  cols    <- c(...)
  missing <- cols[!cols %in% names(data)]
  if (length(missing) > 0)
    stop("Columns not found in data: ", paste(missing, collapse = ", "), call. = FALSE)
}

#' Standardise long-format data
#' @keywords internal
.prep_data <- function(data, yi, vi, study, time) {
  .check_cols(data, yi, vi, study, time)
  dat           <- as.data.frame(data)
  dat$.yi       <- as.numeric(dat[[yi]])
  dat$.vi       <- as.numeric(dat[[vi]])
  dat$.study    <- as.character(dat[[study]])
  dat$.time     <- as.numeric(dat[[time]])
  dat[order(dat$.time, dat$.study), ]
}

#' Partial r from t-stat
#' @keywords internal
.r_partial <- function(t, df) t / sqrt(t^2 + df)
