#' Leave-One-Out and Leave-k-Out Fragility Analysis
#'
#' Computes fragility indices for each time point by systematically removing
#' studies and re-estimating the pooled effect.  The fragility index at time
#' \eqn{t} is the minimum number of studies whose removal changes the
#' statistical conclusion (significant -> non-significant or vice versa).
#'
#' @details
#' At each time point, studies are removed one at a time (or in combinations
#' for the leave-k-out version) and the model is re-fit.  The fragility index
#' is the smallest \eqn{k} such that removing any set of \eqn{k} studies
#' flips the significance of the pooled estimate.  A fragility index of 1
#' means a single study's removal changes the conclusion.
#'
#' For the leave-k-out version, a random sample of combinations is used when
#' the number of combinations is large (controlled by `max_combinations`).
#'
#' @param data      Long-format `data.frame`.
#' @param meta_obj  Output from [ml_meta()].
#' @param yi,vi,study,time  Column names.
#' @param max_k     Maximum number of studies to remove.  Default `5`.
#' @param max_combinations  Maximum number of combinations to test per \eqn{k}.
#'   Default `500`.  Larger values are more exhaustive but slower.
#' @param alpha     Significance level.
#' @param rho       Working correlation.
#' @param small_sample  Use CR2 + Satterthwaite?
#' @param seed      Random seed for sampling combinations.  Default `NULL`.
#'
#' @return Object of class `ml_fragility` (a `data.frame`) with columns:
#'   \describe{
#'     \item{`time`}{Follow-up time.}
#'     \item{`k_studies`}{Number of studies at this time point.}
#'     \item{`p_original`}{Original p-value.}
#'     \item{`sig_original`}{Was the original result significant?}
#'     \item{`fragility_index`}{Min number of removals to flip significance.
#'       `NA` if not found within `max_k`.}
#'     \item{`fragility_quotient`}{`fragility_index / k_studies` (proportion).}
#'     \item{`study_removed`}{Study ID whose removal achieved the flip
#'       (leave-one-out only).}
#'   }
#'
#' @examples
#' \donttest{
#' dat  <- sim_longitudinal_meta(k = 10, times = c(0, 6, 12), seed = 5)
#' meta <- ml_meta(dat, yi = "yi", vi = "vi", study = "study", time = "time")
#' frag <- ml_fragility(dat, meta, yi = "yi", vi = "vi",
#'                       study = "study", time = "time",
#'                       max_k = 1L, seed = 1)
#' print(frag)
#' }
#'
#' @export
ml_fragility <- function(data,
                          meta_obj,
                          yi,
                          vi,
                          study,
                          time,
                          max_k            = 5L,
                          max_combinations = 500L,
                          alpha            = NULL,
                          rho              = 0.8,
                          small_sample     = TRUE,
                          seed             = NULL) {

  if (!inherits(meta_obj, "ml_meta"))
    stop("`meta_obj` must be output from ml_meta().", call. = FALSE)
  alpha <- alpha %||% attr(meta_obj, "alpha") %||% 0.05

  dat       <- .prep_data(data, yi, vi, study, time)
  time_vals <- meta_obj$time
  out_list  <- vector("list", length(time_vals))

  for (j in seq_along(time_vals)) {
    tt      <- time_vals[j]
    dtt     <- dat[dat$.time == tt, ]
    k_t     <- length(unique(dtt$.study))
    res_tt  <- meta_obj[meta_obj$time == tt, ]

    skeleton <- data.frame(
      time               = tt,
      k_studies          = k_t,
      p_original         = res_tt$p_val,
      sig_original       = !is.na(res_tt$p_val) && res_tt$p_val < alpha,
      fragility_index    = NA_integer_,
      fragility_quotient = NA_real_,
      study_removed      = NA_character_,
      stringsAsFactors   = FALSE
    )

    if (is.na(res_tt$p_val) || k_t < 3) {
      out_list[[j]] <- skeleton
      next
    }

    sig_orig  <- res_tt$p_val < alpha
    study_ids <- unique(dtt$.study)
    flip_k    <- NA_integer_
    flip_study <- NA_character_

    # --- Leave-one-out pass first (fast + identifies specific culprit) -----
    for (s_id in study_ids) {
      dtt_sub  <- dtt[dtt$.study != s_id, ]
      k_sub    <- length(unique(dtt_sub$.study))
      if (k_sub < 2) next

      res_sub <- .fit_intercept(dtt_sub$.yi, dtt_sub$.vi, dtt_sub$.study,
                                 rho, small_sample, alpha)
      if (is.null(res_sub)) next

      p_sub <- res_sub$p_val
      if (!is.na(p_sub) && (p_sub < alpha) != sig_orig) {
        flip_k     <- 1L
        flip_study <- s_id
        break
      }
    }

    # --- Leave-k-out for k = 2 .. max_k if not found at k=1 ---------------
    if (is.na(flip_k) && max_k > 1) {
      if (!is.null(seed)) set.seed(seed)

      for (kk in 2:min(max_k, k_t - 2)) {
        combos <- .sample_combinations(study_ids, kk, max_combinations)
        found  <- FALSE

        for (combo in combos) {
          dtt_sub <- dtt[!dtt$.study %in% combo, ]
          k_sub   <- length(unique(dtt_sub$.study))
          if (k_sub < 2) next

          res_sub <- .fit_intercept(dtt_sub$.yi, dtt_sub$.vi, dtt_sub$.study,
                                     rho, small_sample, alpha)
          if (is.null(res_sub)) next

          p_sub <- res_sub$p_val
          if (!is.na(p_sub) && (p_sub < alpha) != sig_orig) {
            flip_k <- kk
            found  <- TRUE
            break
          }
        }
        if (found) break
      }
    }

    skeleton$fragility_index    <- flip_k
    skeleton$fragility_quotient <- if (!is.na(flip_k)) flip_k / k_t else NA_real_
    skeleton$study_removed      <- flip_study
    out_list[[j]]               <- skeleton
  }

  out           <- do.call(rbind, out_list)
  rownames(out) <- NULL
  class(out)    <- c("ml_fragility", "data.frame")
  out
}


# ---- Internal helper: sample combinations --------------------------------

.sample_combinations <- function(ids, k, max_n) {
  n    <- length(ids)
  total <- choose(n, k)
  if (total <= max_n) {
    # Generate all combinations
    utils::combn(ids, k, simplify = FALSE)
  } else {
    # Sample random combinations
    replicate(max_n, sample(ids, k, replace = FALSE), simplify = FALSE)
  }
}


# ---- S3 methods -----------------------------------------------------------

#' @export
print.ml_fragility <- function(x, digits = 3, ...) {
  cat("\n=== metaLong: Fragility Analysis ===\n\n")
  pr         <- as.data.frame(x)
  pr$p_original <- round(pr$p_original, digits)
  pr$fragility_quotient <- round(pr$fragility_quotient, digits)
  print(pr[, c("time", "k_studies", "p_original", "sig_original",
               "fragility_index", "fragility_quotient", "study_removed")],
        row.names = FALSE)
  invisible(x)
}

#' @export
plot.ml_fragility <- function(x,
                               main = "Fragility Index Over Time",
                               col_sig    = "#2166ac",
                               col_nonsig = "#d73027",
                               ...) {
  ok   <- !is.na(x$fragility_index)
  if (!any(ok)) { message("No fragility indices to plot."); return(invisible(x)) }

  cols <- ifelse(x$sig_original[ok], col_sig, col_nonsig)
  ylim <- c(0, max(x$fragility_index[ok], na.rm = TRUE) * 1.2)

  barplot(
    x$fragility_index[ok],
    names.arg = x$time[ok],
    col       = cols,
    ylim      = ylim,
    xlab      = "Time",
    ylab      = "Fragility Index",
    main      = main,
    ...
  )
  abline(h = 1, lty = 2, col = "grey50")
  legend("topright",
         legend = c("Originally significant", "Not significant"),
         fill   = c(col_sig, col_nonsig), bty = "n", cex = 0.8)
  invisible(x)
}
