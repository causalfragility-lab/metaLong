#' Internal clubSandwich shim
#'
#' Drop-in replacements for clubSandwich functions when the package is absent.
#' The real package is always preferred when available.
#'
#' @keywords internal
#' @noRd

# ---- Dispatchers ---------------------------------------------------------

.cs_impute_cov <- function(vi, cluster, r) {
  if (requireNamespace("clubSandwich", quietly = TRUE)) {
    return(clubSandwich::impute_covariance_matrix(
      vi = vi, cluster = cluster, r = r))
  }
  .shim_impute_covariance_matrix(vi = vi, cluster = cluster, r = r)
}

.cs_coef_test <- function(obj, vcov = "CR2", cluster,
                           test = "Satterthwaite") {
  if (requireNamespace("clubSandwich", quietly = TRUE)) {
    return(clubSandwich::coef_test(obj, vcov = vcov, cluster = cluster,
                                   test = test))
  }
  .shim_coef_test(obj, cluster = cluster)
}

# ---- Shim: impute_covariance_matrix --------------------------------------

.shim_impute_covariance_matrix <- function(vi, cluster, r) {
  cluster <- as.character(cluster)
  n       <- length(vi)
  V       <- matrix(0.0, n, n)
  for (id in unique(cluster)) {
    idx         <- which(cluster == id)
    cross       <- base::outer(sqrt(vi[idx]), sqrt(vi[idx])) * r
    diag(cross) <- vi[idx]
    V[idx, idx] <- cross
  }
  V
}

# ---- Shim: coef_test (CR1S) ----------------------------------------------
# Uses only reliable rma.mv slots: b, X, yi, vi.
# Avoids obj$V -- may be sparse S4 Matrix.

.shim_coef_test <- function(obj, cluster) {
  b      <- as.numeric(obj$b)
  p      <- length(b)
  X      <- obj$X
  yi_raw <- as.numeric(obj$yi)
  vi_raw <- as.numeric(obj$vi)

  if (is.null(X) || length(yi_raw) == 0 || length(vi_raw) == 0)
    return(NULL)

  cluster <- as.character(cluster)
  n       <- length(yi_raw)

  if (length(cluster) != n || nrow(X) != n)
    return(NULL)

  ids <- unique(cluster)
  g   <- length(ids)

  # Inverse-variance weights
  w   <- 1.0 / vi_raw

  # Residuals from REML fit
  res <- yi_raw - as.numeric(X %*% b)

  # Bread: (X'WX)^{-1}
  XtWX  <- t(X * w) %*% X
  bread <- tryCatch(solve(XtWX), error = function(e) .shim_ginv(XtWX))

  # Meat: sum_g score_g score_g'  (CR1S scaled)
  meat <- matrix(0.0, p, p)
  for (id in ids) {
    idx     <- which(cluster == id)
    score_g <- as.numeric(t(X[idx, , drop = FALSE]) %*% (w[idx] * res[idx]))
    meat    <- meat + base::outer(score_g, score_g)
  }
  if (g > 1L && n > p)
    meat <- (g / (g - 1L)) * ((n - 1L) / (n - p)) * meat

  # Sandwich variance, SE, df, p
  V_sand <- bread %*% meat %*% bread
  SE     <- sqrt(pmax(0.0, diag(V_sand)))
  df_val <- as.numeric(pmax(1L, g - p))
  t_stat <- b / SE
  p_val  <- 2.0 * stats::pt(abs(t_stat), df = df_val, lower.tail = FALSE)

  rn <- if (!is.null(rownames(obj$b))) rownames(obj$b) else paste0("b", seq_len(p))
  data.frame(
    Coef   = b,
    SE     = SE,
    tstat  = t_stat,
    df     = rep(df_val, p),
    p_Satt = p_val,
    row.names = rn,
    stringsAsFactors = FALSE
  )
}

# Moore-Penrose pseudoinverse (internal, no external dependency)
.shim_ginv <- function(X, tol = sqrt(.Machine$double.eps)) {
  s   <- base::svd(X)
  pos <- s$d > (tol * max(s$d))
  if (all(pos)) return(s$v %*% (1.0 / s$d * t(s$u)))
  s$v[, pos, drop = FALSE] %*%
    ((1.0 / s$d[pos]) * t(s$u[, pos, drop = FALSE]))
}
