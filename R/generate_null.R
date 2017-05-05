# This file is part of the `locus` R package:
#     https://github.com/hruffieux/locus
#

#' Generate empirical null distribution using permuted data.
#'
#' This function runs \code{\link{locus}} on data with permuted response sample
#' labels. A common use is for estimatimation of data-driven false discovery
#' rates of given thresholds on the posterior probabilities of inclusion.
#'
#' @param n_perm Number of permuted datasets on which \code{\link{locus}} is
#'   run.
#' @param Y Response data matrix (without permuted sample indices) of size
#'   n x d, where n is the number of observations and d is the number of
#'   response variables.
#' @param X Input matrix of size n x p, where p is the number of candidate
#'   predictors. \code{X} cannot contain NAs. No intercept must be supplied.
#' @param p0_av Prior average number of predictors expected to be included in
#'   the model (if \code{list_blocks} is \code{TRUE}, average number per
#'   predictor blocks). Must be \code{NULL} if \code{list_init} and
#'   \code{list_hyper} are both non-\code{NULL} or if \code{list_cv} is
#'   non-\code{NULL}. If \code{list_blocks} is \code{NULL}, can also be a vector
#'   of length p with entry s corresponding to the prior probability that
#'   candidate predictor s is associated with at least one response.
#' @param Z Covariate matrix of size n x q, where q is the number of
#'   covariates. \code{NULL} if no covariate. Factor covariates must be supplied
#'   after transformation to dummy coding. No intercept must be supplied.
#' @param link Response link. Must be "\code{identity}" for linear regression,
#'   "\code{logit}" for logistic regression, "\code{probit}" for probit
#'   regression, or "\code{mix}" for a mix of identity and probit link functions
#'   (in this case, the indices of the binary responses must be gathered in
#'   argument \code{ind_bin}, see below).
#' @param ind_bin If \code{link = "mix"}, vector of indices corresponding to
#'   the binary variables in \code{Y}. Must be \code{NULL} if
#'   \code{link != "mix"}.
#' @param list_hyper An object of class "\code{hyper}" containing the model
#'   hyperparameters. Must be filled using the \code{\link{set_hyper}}
#'   function or must be \code{NULL} for default hyperparameters.
#' @param list_init An object of class "\code{init}" containing the initial
#'   variational parameters. Must be filled using the \code{\link{set_init}}
#'   function or be \code{NULL} for a default initialization.
#' @param list_blocks An object of class "\code{blocks}" containing setting for
#'   parallel inference on a partitioned predictor space. Must be filled using
#'   the \code{\link{set_blocks}} function or must be \code{NULL} for no
#'   partitioning.
#' @param user_seed Seed set for reproducible default choices of hyperparameters
#'   (if \code{list_hyper} is \code{NULL}) and inital variational parameters (if
#'   \code{list_init} is \code{NULL}). Default is \code{NULL}, no seed set.
#' @param tol Tolerance for the stopping criterion.
#' @param maxit Maximum number of iterations allowed.
#' @param verbose If \code{TRUE}, messages are displayed during execution.
#' @param results_dir Path where the output of each of the \code{n_perm} runs
#'   will be saved. Default is \code{NULL}, the output is not saved to files and
#'   a list of objects of class "\code{perm}" is returned instead.
#' @param n_cpus Number of CPUs to be used. If large, one should ensure that
#'   enough RAM will be available for parallel execution. Set to 1 for serial
#'   execution.
#'
#' @return If \code{results_dir} is \code{NULL}, list of length \code{n_perm}
#'   with each element corresponding to the output of \code{\link{locus}} on one
#'   permuted dataset. Each contains:
#'  \item{ind_perm}{Vector of length n containing the permuted response sample
#'                  labels.}
#'  \item{gam_vb}{Posterior inclusion probability matrix of size p x d.
#'                Entry (s, t) corresponds to the posterior probability of
#'                association between predictor s and response t.}
#'  \item{om_vb}{Vector of length p containing the posterior mean of omega.
#'               Entry s controls the proportion of responses associated with
#'               candidate predictor s.}
#
#' @seealso \code{\link{locus}}
#'
#' @examples
#' seed <- 123; set.seed(seed)
#'
#' ###################
#' ## Simulate data ##
#' ###################
#'
#' ## Example using small problem sizes:
#' ##
#' n <- 200; p <- 250; p0 <- 20; d <- 25; d0 <- 20
#'
#' ## Candidate predictors (subject to selection)
#' ##
#' # Here we simulate common genetic variants (but any type of candidate
#' # predictors can be supplied).
#' # 0 = homozygous, major allele, 1 = heterozygous, 2 = homozygous, minor allele
#' #
#' X_act <- matrix(rbinom(n * p0, size = 2, p = 0.25), nrow = n)
#' X_inact <- matrix(rbinom(n * (p - p0), size = 2, p = 0.25), nrow = n)
#'
#' shuff_x_ind <- sample(p)
#' X <- cbind(X_act, X_inact)[, shuff_x_ind]
#'
#' bool_x_act <- shuff_x_ind <= p0
#'
#' pat_act <- beta <- matrix(0, nrow = p0, ncol = d0)
#' pat_act[sample(p0*d0, floor(p0*d0/5))] <- 1
#' beta[as.logical(pat_act)] <-  rnorm(sum(pat_act))
#'
#' ## Gaussian responses
#' ##
#' Y_act <- matrix(rnorm(n * d0, mean = X_act %*% beta, sd = 0.5), nrow = n)
#' Y_inact <- matrix(rnorm(n * (d - d0), sd = 0.5), nrow = n)
#' shuff_y_ind <- sample(d)
#' Y <- cbind(Y_act, Y_inact)[, shuff_y_ind]
#'
#' ########################
#' ## Infer associations ##
#' ########################
#'
#' list_perm <- generate_null(n_perm = 10, Y, X, p0_av = p0, link = "identity",
#'                            user_seed = seed, verbose = FALSE)
#'
#' @export
#'
generate_null <- function(n_perm, Y, X, p0_av, Z = NULL, link = "identity",
                          ind_bin = NULL, list_hyper = NULL, list_init = NULL,
                          list_blocks = NULL, list_groups = NULL,
                          user_seed = NULL, tol = 1e-3, maxit = 1000,
                          verbose = TRUE, results_dir = NULL, n_cpus = 1) {

  if (!is.null(user_seed)){
    RNGkind("L'Ecuyer-CMRG") # ensure reproducibility when using mclapply
    set.seed(user_seed)
    if (verbose) cat(paste("Seed set to ", user_seed, ". \n", sep = ""))
  }

  check_structure_(n_perm, "vector", "numeric", 1)
  check_natural_(n_perm)

  check_structure_(n_cpus, "vector", "numeric", 1)
  check_natural_(n_cpus)

  if (n_cpus > 1) {

    n_cpus_avail <- parallel::detectCores()
    if (n_cpus > n_cpus_avail) {
      n_cpus <- n_cpus_avail
      warning(paste("The number of CPUs specified exceeds the number of CPUs ",
                    "available on the machine. The latter has been used instead.", sep=""))
    }
    if (verbose) cat(paste("Permutations with ", n_cpus, " CPUs.\n",
                           "Please make sure that enough RAM is available. \n", sep=""))
  }

  check_structure_(results_dir, "vector", "string", 1, null_ok = TRUE)

  n <- nrow(Y)

  permute <- function(i) {
    ind_perm <- sample(1:n) # random permutation
    rownames(Y) <- NULL

    # user_seed must be NULL here otherwise always the same permutation
    res_perm <- locus(Y = Y[ind_perm, ], X = X, p0_av = p0_av, Z = Z,
                      link = link, ind_bin = ind_bin,
                      list_hyper = list_hyper, list_init = list_init,
                      list_cv = NULL, list_blocks = list_blocks,
                      list_groups = list_groups, user_seed = NULL, tol = tol,
                      maxit = maxit, verbose = verbose)

    om_vb <- res_perm$om_vb
    gam_vb <- res_perm$gam_vb
    rm(res_perm)

    if (is.null(results_dir)) {

      create_named_list_(ind_perm, gam_vb, om_vb)

    } else {

      save(ind_perm, gam_vb, om_vb,
           file = paste(results_dir, "vb_real_data_", i, ".RData", sep=""))
      rm(om_vb)
      rm(gam_vb)
      NULL

    }
  }

  parallel::mclapply(1:n_perm, function(i) permute(i), mc.cores = n_cpus)

}
