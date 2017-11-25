#' locus: a package for combined predictor and outcome selection in
#' high-dimensional set-ups using variational inference
#'
#' The locus package provides an efficient variational algorithm for
#' simultaneous variable selection of predictors and associated outcomes based
#' on a sparse multivariate regression model (H. Ruffieux, A. C. Davison,
#' J. Hager, I. Irincheeva, Efficient inference for genetic association studies
#' with multiple outcomes, Biostatistics, 2017). The methods from this package
#' have been used on large genetic datasets from molecular quantitative trait
#' locus (QTL) problems with over 200K single nucleotide polymorphisms (SNPs),
#' hundreds of molecular expression levels and hundreds of samples.
#'
#' @section locus functions: set_hyper, set_init, generate_null, locus,
#'   set_blocks, set_cv, set_groups, set_struct.
#'
#' @docType package
#' @name locus-package
#' @useDynLib locus, .registration = TRUE
#' @import RcppEigen
#' @importFrom Rcpp evalCpp
#' @importFrom stats cor dnorm median pnorm qnorm rbeta rbinom rgamma rnorm setNames uniroot var
NULL
