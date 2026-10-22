#' SDR: Formulas for Successive Difference Replication Under SRSWOR
#'
#' @description
#' This package offers formulas for evaluating bias and variance of Successive Difference Replication under stratified simple random sampling without replacement. The D-cycling definitions, groups of sample units sharing the same weight-replicate factors, matrix M defining SDR-variance-estimator as a quadratic form, and various formulas, can all be found in the main manuscript.
#'
#' Additional examples can be found in the `SDRpkg.RLog` file, either in the inst/ folder of the SDR package or within R via `system.file("SDRpkg.RLog", package="SDR")` after package installation.
#' 
#' @section Notation:
#'
#' | symbol/notation | description |
#' | :---   | :--- |
#' | H | Total number of stratum in population |
#' | h | index given for a particular stratum |
#' | Nh | Stratum total population for stratum h |
#' | nh | Sample size for stratum h |
#' | muh | Superpopulation parameter mean for stratum h |
#' | sigmah | Superpopulation parameter standard deviation for stratum h |
#' | alphah | Domain superpopulation mean for stratum h |
#'
#' @references
#' Slud, E., & Trudell, T. (2026). SDR Variance Estimates in Small Domains. arXiv preprint arXiv:2608.17353. doi:10.48550/arXiv.2608.17353
#'
#' @docType package
#' @name SDR
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats aggregate quantile var rnorm
## usethis namespace: end
NULL
