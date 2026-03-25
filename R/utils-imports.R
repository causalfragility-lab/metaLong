# utils-imports.R
# All importFrom declarations live here so roxygen2 never drops them
# when regenerating the NAMESPACE.

#' @importFrom stats qt pt pf var sd weighted.mean na.omit setNames
#'   coef vcov predict model.matrix quantile runif rnorm rbinom
#' @importFrom graphics plot barplot abline lines polygon points legend
#'   mtext par text title axis arrows
#' @importFrom grDevices rgb col2rgb colorRampPalette
#' @importFrom utils head tail combn
#' @importFrom splines ns bs
#' @importFrom metafor rma.mv rma.uni
NULL

#' @importFrom graphics plot.new
NULL
