## R CMD check results

Status: 0 errors | 0 warnings | 0 notes

Tested on:
- Windows 11 x64, R 4.5.1

## Test environments

- Local: Windows 11 x64 (build 26200), R 4.5.1 (2025-06-13 ucrt)
- win-builder: (to be run before submission)

## Reverse dependencies

This is a new package. There are no existing reverse dependencies.

## Comments to CRAN

This is a first submission of a new package.

The package implements longitudinal meta-analysis methods including:
- Robust variance estimation with Tipton small-sample corrections
- Time-varying sensitivity analysis via ITCV (Frank, 2000)
- Benchmark calibration of confounding thresholds
- Spline-based nonlinear time-trend modeling

All examples run in under 1 second.
Slow functions (ml_benchmark, ml_fragility) are wrapped in \donttest{}.

The NOTE 'unable to verify current time' appears only in local checks
due to a network restriction and is not present on win-builder.

The quarto warning in devtools::check() output is a local Windows
environment issue (path with spaces) and is not part of the package.
