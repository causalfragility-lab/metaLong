# metaLong

**Longitudinal Meta-Analysis with Robust Variance Estimation and Sensitivity Analysis**

## Installation

```r
# From CRAN (once released)
install.packages("metaLong")

# Development version from GitHub
remotes::install_github("causalfragility-lab/metaLong")
```

## Overview

metaLong provides tools for synthesising evidence from studies that report
outcomes at multiple follow-up time points. Core functions:

| Function | Purpose |
|---|---|
| `ml_meta()` | Pool effects at each time point with RVE + Tipton correction |
| `ml_sens()` | Time-varying ITCV sensitivity analysis |
| `ml_benchmark()` | Benchmark ITCV against observed covariates |
| `ml_spline()` | Nonlinear spline time trend |
| `ml_fragility()` | Leave-k-out fragility analysis |
| `ml_plot()` | Combined publication figure |
| `sim_longitudinal_meta()` | Simulate longitudinal meta-analytic data |

## Basic usage

```r
library(metaLong)

dat  <- sim_longitudinal_meta(k = 20, times = c(0, 6, 12, 24), seed = 42)
meta <- ml_meta(dat, yi = "yi", vi = "vi", study = "study", time = "time")
sens <- ml_sens(dat, meta, yi = "yi", vi = "vi", study = "study", time = "time")
ml_plot(meta, sens_obj = sens)
```

## References

- Frank (2000) <doi:10.1177/0049124100029002003>
- Hedges, Tipton & Johnson (2010) <doi:10.1002/jrsm.5>
- Tipton (2015) <doi:10.1037/met0000011>

## Author

Subir Hait, Michigan State University (<haitsubi@msu.edu>)
ORCID: 0009-0004-9871-9677
