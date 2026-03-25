#' Benchmark Calibration of Longitudinal ITCV Against Observed Covariates
#'
#' For each follow-up time point, regresses each observed study-level covariate
#' on the effect sizes using RVE meta-regression, extracts the covariate's
#' partial correlation with the outcome, and compares it to the
#' significance-adjusted ITCV threshold from [ml_sens()].  A covariate that
#' *beats* the threshold demonstrates that real-world confounding of at least
#' that magnitude exists, which is direct evidence of effect fragility.
#'
#' @section Interpretation:
#' If an observed covariate (e.g., publication year, sample quality, attrition
#' rate) has `|r_partial| >= ITCV_alpha(t)`, then an *unobserved* confounder with
#' the same relationship to exposure and outcome would be sufficient to nullify
#' the pooled effect at time \eqn{t}.  This does not prove confounding--it
#' calibrates the plausibility threshold.
#'
#' @param data       Long-format `data.frame`.
#' @param meta_obj   Output from [ml_meta()].
#' @param sens_obj   Output from [ml_sens()].
#' @param yi,vi,study,time  Column names.
#' @param covariates Character vector of observed moderator column names to benchmark.
#' @param alpha      Significance level (inherits from `meta_obj` if `NULL`).
#' @param rho        Working within-study correlation for V matrix.
#' @param small_sample Logical; use CR2 + Satterthwaite?
#' @param min_k      Minimum studies required at a time point.  Default `3L`
#'   (one extra relative to [ml_meta()] because regression needs more d.f.).
#'
#' @return Object of class `ml_benchmark` (a `data.frame`) with columns:
#'   \describe{
#'     \item{`time`}{Follow-up time.}
#'     \item{`covariate`}{Covariate name.}
#'     \item{`k`}{Number of studies.}
#'     \item{`r_partial`}{Partial correlation of covariate with effect size.}
#'     \item{`t_stat`, `df`, `p_val`}{RVE inference for the covariate slope.}
#'     \item{`itcv_alpha`}{ITCV_alpha threshold at this time point.}
#'     \item{`beats_threshold`}{Logical: does `|r_partial| >= itcv_alpha`?}
#'     \item{`skip_reason`}{Character; reason a cell was skipped, else `NA`.}
#'   }
#'   The `"fragile_summary"` attribute contains one row per time with counts.
#'
#' @seealso [ml_sens()], [ml_meta()]
#'
#' @examples
#' \donttest{
#' dat   <- sim_longitudinal_meta(k = 15, times = c(0, 6, 12), seed = 2)
#' meta  <- ml_meta(dat, yi = "yi", vi = "vi", study = "study", time = "time")
#' sens  <- ml_sens(dat, meta, yi = "yi", vi = "vi", study = "study", time = "time")
#' bench <- ml_benchmark(dat, meta, sens,
#'                        yi = "yi", vi = "vi", study = "study", time = "time",
#'                        covariates = c("pub_year", "quality"))
#' print(bench)
#' plot(bench)
#' }
#'
#' @export
ml_benchmark <- function(data,
                          meta_obj,
                          sens_obj,
                          yi,
                          vi,
                          study,
                          time,
                          covariates,
                          alpha        = NULL,
                          rho          = 0.8,
                          small_sample = TRUE,
                          min_k        = 3L) {

  # ---- Validation --------------------------------------------------------
  if (!inherits(meta_obj, "ml_meta"))
    stop("`meta_obj` must be output from ml_meta().", call. = FALSE)
  if (!inherits(sens_obj, "ml_sens"))
    stop("`sens_obj` must be output from ml_sens().", call. = FALSE)
  if (!is.character(covariates) || length(covariates) == 0)
    stop("`covariates` must be a non-empty character vector.", call. = FALSE)

  alpha <- alpha %||% attr(meta_obj, "alpha") %||% 0.05

  missing_covs <- covariates[!covariates %in% names(data)]
  if (length(missing_covs) > 0)
    stop("Covariates not found in data: ", paste(missing_covs, collapse = ", "),
         call. = FALSE)

  dat       <- .prep_data(data, yi, vi, study, time)
  time_vals <- sort(unique(dat$.time))
  results   <- vector("list", length(time_vals))

  for (j in seq_along(time_vals)) {
    tt          <- time_vals[j]
    dtt         <- dat[dat$.time == tt, ]
    k_t         <- length(unique(dtt$.study))
    sens_row    <- sens_obj[sens_obj$time == tt, ]
    itcv_alpha_t <- if (nrow(sens_row) == 1) sens_row$itcv_alpha else NA_real_

    cov_results <- vector("list", length(covariates))

    for (cv_idx in seq_along(covariates)) {
      cov_name <- covariates[cv_idx]
      cov_vals <- dtt[[cov_name]]

      skeleton <- data.frame(
        time             = tt,
        covariate        = cov_name,
        k                = k_t,
        r_partial        = NA_real_,
        t_stat           = NA_real_,
        df               = NA_real_,
        p_val            = NA_real_,
        itcv_alpha       = itcv_alpha_t,
        beats_threshold  = NA,
        skip_reason      = NA_character_,
        stringsAsFactors = FALSE
      )

      # Skip conditions
      skip_reason <- .check_covariate(cov_vals, k_t, min_k)
      if (!is.na(skip_reason)) {
        skeleton$skip_reason <- skip_reason
        cov_results[[cv_idx]] <- skeleton
        next
      }

      # Build V
      V <- tryCatch(
        .build_V(dtt$.vi, dtt$.study, rho),
        error = function(e) NULL
      )
      if (is.null(V)) {
        skeleton$skip_reason <- "V_matrix_failed"
        cov_results[[cv_idx]] <- skeleton
        next
      }

      # Centre covariate for interpretability
      cov_c <- cov_vals - mean(cov_vals, na.rm = TRUE)

      # Fit meta-regression
      fit <- tryCatch(
        suppressWarnings(
          metafor::rma.mv(
            yi     = dtt$.yi,
            V      = V,
            mods   = ~ cov_c,
            method = "REML",
            sparse = TRUE
          )
        ),
        error = function(e) NULL
      )

      if (is.null(fit)) {
        skeleton$skip_reason <- "rma.mv_failed"
        cov_results[[cv_idx]] <- skeleton
        next
      }

      # RVE inference
      if (small_sample) {
        cs <- tryCatch(
          .cs_coef_test(
            fit,
            vcov    = "CR2",
            cluster = dtt$.study,
            test    = "Satterthwaite"
          ),
          error = function(e) NULL
        )
        if (is.null(cs) || nrow(cs) < 2) {
          skeleton$skip_reason <- "CR2_failed"
          cov_results[[cv_idx]] <- skeleton
          next
        }
        t_stat <- cs$tstat[2]
        df_cov <- cs$df[2]
        p_val  <- cs$p_Satt[2]
      } else {
        t_stat <- fit$zval[2]
        df_cov <- Inf
        p_val  <- fit$pval[2]
      }

      r_partial_val   <- .r_partial(t_stat, df_cov)
      beats_threshold <- if (!is.na(itcv_alpha_t) && !is.na(r_partial_val)) {
        abs(r_partial_val) >= itcv_alpha_t
      } else {
        NA
      }

      cov_results[[cv_idx]] <- data.frame(
        time             = tt,
        covariate        = cov_name,
        k                = k_t,
        r_partial        = r_partial_val,
        t_stat           = t_stat,
        df               = df_cov,
        p_val            = p_val,
        itcv_alpha       = itcv_alpha_t,
        beats_threshold  = beats_threshold,
        skip_reason      = NA_character_,
        stringsAsFactors = FALSE
      )
    }

    results[[j]] <- do.call(rbind, cov_results)
  }

  out           <- do.call(rbind, results)
  rownames(out) <- NULL

  # ---- Fragility summary -------------------------------------------------
  fragile_summary <- do.call(rbind, lapply(time_vals, function(tt) {
    sub      <- out[out$time == tt & !is.na(out$beats_threshold), ]
    n_beats  <- sum(sub$beats_threshold, na.rm = TRUE)
    n_tested <- nrow(sub)
    data.frame(
      time         = tt,
      itcv_alpha   = unique(out$itcv_alpha[out$time == tt])[1],
      n_covariates = n_tested,
      n_beats      = n_beats,
      prop_beats   = if (n_tested > 0) n_beats / n_tested else NA_real_,
      any_beats    = n_beats > 0,
      stringsAsFactors = FALSE
    )
  }))

  attr(out, "fragile_summary") <- fragile_summary
  attr(out, "alpha")           <- alpha
  class(out) <- c("ml_benchmark", "data.frame")
  out
}


# ---- Internal helper -----------------------------------------------------

.check_covariate <- function(vals, k_t, min_k) {
  if (k_t < min_k)          return(paste0("k < min_k (", k_t, ")"))
  if (all(is.na(vals)))     return("all_na_covariate")
  if (sum(!is.na(vals)) < 3) return("too_few_non_na")
  if (isTRUE(var(vals, na.rm = TRUE) < .Machine$double.eps))
                              return("zero_variance_covariate")
  NA_character_
}


# ---- S3 methods -----------------------------------------------------------

#' @export
print.ml_benchmark <- function(x, digits = 3, ...) {
  cat("\n=== metaLong: Benchmark Calibration ===\n\n")
  fs <- attr(x, "fragile_summary")
  if (!is.null(fs)) {
    cat("Fragility summary by time:\n")
    fs[, -1] <- lapply(fs[, -1], function(col) {
      if (is.numeric(col)) round(col, digits) else col
    })
    print(fs, row.names = FALSE)
    cat("\n")
  }

  pr       <- as.data.frame(x)
  num_cols <- c("r_partial", "t_stat", "df", "p_val", "itcv_alpha")
  pr[num_cols] <- lapply(pr[num_cols], round, digits = digits)
  print(pr[, c("time", "covariate", "k", "r_partial", "itcv_alpha",
               "beats_threshold", "p_val")], row.names = FALSE)
  invisible(x)
}

#' @export
plot.ml_benchmark <- function(x,
                               main  = "Benchmark: |r_partial| vs ITCV_alpha",
                               col_beat   = "#d73027",
                               col_nobeat = "#4575b4",
                               ...) {
  time_vals <- sort(unique(x$time))
  n_t       <- length(time_vals)
  ncols     <- min(3L, n_t)
  nrows     <- ceiling(n_t / ncols)

  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow = c(nrows, ncols), mar = c(5, 4, 3, 1))

  for (tt in time_vals) {
    sub       <- x[x$time == tt & !is.na(x$r_partial), ]
    threshold <- unique(sub$itcv_alpha)[1]
    if (nrow(sub) == 0) next

    abs_r  <- abs(sub$r_partial)
    colors <- ifelse(isTRUE(sub$beats_threshold), col_beat, col_nobeat)
    ylim   <- c(0, max(c(abs_r, threshold), na.rm = TRUE) * 1.2)

    bp <- barplot(
      abs_r,
      names.arg = sub$covariate,
      col       = colors,
      ylim      = ylim,
      las       = 2,
      main      = paste0("t = ", tt),
      ylab      = "|r partial|",
      cex.names = 0.80
    )
    if (!is.na(threshold)) {
      abline(h = threshold, lty = 2, lwd = 1.5)
      text(max(bp) * 0.98, threshold,
           labels = paste0("ITCV_a=", round(threshold, 2)),
           adj    = c(1, -0.3), cex = 0.70)
    }
  }
  mtext(main, outer = TRUE, line = -1.2, cex = 1.0)
  invisible(x)
}

#' @export
tidy.ml_benchmark <- function(x, ...) {
  out <- as.data.frame(x)[, c("time", "covariate", "k", "r_partial",
                                "t_stat", "df", "p_val",
                                "itcv_alpha", "beats_threshold")]
  out$interpretation <- dplyr_free_case(
    is.na(out$beats_threshold), "insufficient data",
    out$beats_threshold,        "exceeds threshold -- effect fragile",
    TRUE,                       "below threshold -- effect robust"
  )
  out
}

# Minimal case_when without dplyr
dplyr_free_case <- function(...) {
  args <- list(...)
  out  <- rep(NA_character_, length(args[[1]]))
  for (i in seq(1, length(args) - 1, by = 2)) {
    cond <- args[[i]]; val <- args[[i + 1]]
    out[cond & is.na(out)] <- val
  }
  out
}
