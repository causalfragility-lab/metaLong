#' Simulate a Longitudinal Meta-Analytic Dataset
#'
#' Generates a synthetic long-format dataset suitable for testing and
#' illustrating all `metaLong` functions.  Studies contribute effect sizes at
#' multiple follow-up time points with within-study correlation.
#'
#' @details
#' The true effect at time \eqn{t} for study \eqn{i} is
#' \deqn{\theta_{it} = \mu_t + u_i + \epsilon_{it}}
#' where \eqn{\mu_t} is a time-varying mean effect (optionally nonlinear),
#' \eqn{u_i \sim N(0, \tau^2)} is a study-level random effect, and
#' \eqn{\epsilon_{it} \sim N(0, v_{it})} is sampling error.  Within-study
#' correlation between time points is introduced through \eqn{u_i}.
#'
#' @param k          Number of studies.  Default `20`.
#' @param times      Numeric vector of follow-up time points.
#'   Default `c(0, 6, 12, 24)`.
#' @param mu         Named numeric vector of true effects at each time point,
#'   or a single value (recycled).  Default `0.4`.
#' @param tau        Between-study SD.  Default `0.2`.
#' @param v_range    Two-element vector for the uniform sampling variance range.
#'   Default `c(0.02, 0.12)`.
#' @param missing_prop  Proportion of study x time combinations to set missing
#'   (simulates unbalanced follow-up).  Default `0.0`.
#' @param add_covariates Logical.  If `TRUE`, adds study-level covariates
#'   `pub_year`, `quality`, and `n` for use with [ml_benchmark()].
#'   Default `TRUE`.
#' @param seed       Random seed.  Default `NULL`.
#'
#' @return A `data.frame` in long format with columns:
#'   \describe{
#'     \item{`study`}{Study identifier (character).}
#'     \item{`time`}{Follow-up time.}
#'     \item{`yi`}{Observed effect size.}
#'     \item{`vi`}{Sampling variance.}
#'     \item{`pub_year`, `quality`, `n`}{Study-level covariates
#'       (if `add_covariates = TRUE`).}
#'   }
#'
#' @examples
#' dat <- sim_longitudinal_meta(k = 10, times = c(0, 6, 12), seed = 42)
#' head(dat)
#'
#' # Nonlinear true trajectory
#' \donttest{
#' mu_t <- c("0" = 0.2, "6" = 0.5, "12" = 0.4, "24" = 0.1)
#' dat2 <- sim_longitudinal_meta(k = 10, times = c(0, 6, 12, 24), mu = mu_t,
#'                                missing_prop = 0.1, seed = 99)
#' }
#'
#' @export
sim_longitudinal_meta <- function(k              = 20L,
                                   times          = c(0, 6, 12, 24),
                                   mu             = 0.4,
                                   tau            = 0.2,
                                   v_range        = c(0.02, 0.12),
                                   missing_prop   = 0.0,
                                   add_covariates = TRUE,
                                   seed           = NULL) {

  if (!is.null(seed)) set.seed(seed)

  times <- sort(unique(as.numeric(times)))
  T_    <- length(times)

  # Recycle or name mu across time points
  if (length(mu) == 1) {
    mu <- setNames(rep(mu, T_), as.character(times))
  } else {
    if (is.null(names(mu))) names(mu) <- as.character(times)
    mu <- mu[as.character(times)]
    if (anyNA(mu))
      stop("`mu` names must match `times`.", call. = FALSE)
  }

  # Study-level random effects
  u_i <- stats::rnorm(k, 0, tau)

  # Build long-format grid
  grid <- base::expand.grid(study = seq_len(k), time = times,
                      KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE)
  grid$study <- paste0("s", sprintf("%02d", grid$study))

  # Sampling variances
  grid$vi <- stats::runif(nrow(grid), v_range[1], v_range[2])

  # True effects + sampling error
  u_vec      <- u_i[as.integer(sub("s", "", grid$study))]
  mu_vec     <- mu[as.character(grid$time)]
  grid$yi    <- mu_vec + u_vec + stats::rnorm(nrow(grid), 0, sqrt(grid$vi))

  # Add study-level covariates
  if (add_covariates) {
    study_ids <- unique(grid$study)
    cov_df    <- data.frame(
      study    = study_ids,
      pub_year = as.integer(stats::runif(length(study_ids), 2000, 2022)),
      quality  = round(stats::rnorm(length(study_ids), 0, 1), 2),
      n        = as.integer(stats::runif(length(study_ids), 30, 500)),
      stringsAsFactors = FALSE
    )
    grid <- merge(grid, cov_df, by = "study", sort = FALSE)
  }

  # Introduce missing observations
  if (missing_prop > 0) {
    n_miss <- round(nrow(grid) * missing_prop)
    drop   <- sample(nrow(grid), n_miss)
    grid   <- grid[-drop, ]
  }

  grid <- grid[order(grid$study, grid$time), ]
  rownames(grid) <- NULL
  grid
}
