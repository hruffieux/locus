#' locus: a package for combined predictor and outcome selection in
#' high-dimensional set-ups using variational inference
#'
#' The locus package provides an efficient variational algorithm for
#' simultaneous variable selection of predictors and associated outcomes based
#' on a sparse multivariate regression model (Helene Ruffieux, Anthony C.
#' Davison, Jorg Hager, Irina Irincheeva, 2016, arXiv:1609.03400). This software
#' on large genetic datasets from molecular quantitative trait locus (QTL)
#' problems with over 200K single nucleotide polymorphisms (SNPs), hundreds of
#' molecular expression levels and hundreds of samples.
#'
#' @section locus functions: feed_hyperparam, feed_init_param,
#' generate_dependence, generate_null, generate_phenos, generate_snps, locus,
#' replicate_real_phenos, replicate_real_snps, set_blocks, set_cv.
#'
#' @docType package
#' @name locus-package
#' @importFrom stats cor median qnorm rbeta rbinom rgamma rnorm runif setNames var
NULL
