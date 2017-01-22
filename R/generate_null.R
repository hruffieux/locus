#' Generate empirical null distribution using permuted data.
#'
#' This function runs \code{\link{locus}} on data with permuted response sample
#' labels. A common use is for estimatimation of data-driven false discovery
#' rates of given thresholds on the posterior probabilities of inclusion.
#'
#' @param n_perm Number of permuted datasets on which \code{\link{locus}} is run.
#' @param Y Continuous response data matrix (without permuted sample indices)
#'   of dimension n x d, where n is the number of observations and d is the
#'   number of response variables.
#' @param X Input matrix of dimension n x p, where p is the number of candidate
#'   predictors. \code{X} cannot contain NAs. No intercept must be supplied.
#' @param p0_av Prior average number of predictors expected to be included in
#'   the model (if \code{list_blocks} is \code{TRUE}, average number per
#'   predictor blocks). Must be \code{NULL} if \code{list_init} and
#'   \code{list_hyper} are both non-\code{NULL} or if \code{list_cv} is
#'   non-\code{NULL}. If \code{list_blocks} is \code{NULL}, can also be a vector
#'   of length p with entry s corresponding to the prior probability that
#'   candidate predictor s is associated with at least one response. If
#'   \code{list_blocks} is non-\code{NULL}, can be a vector of size given by the
#'   number of blocks, with each entry corresponding to the prior average number
#'   of predictors from each block expected to be included in the model.
#' @param Z Covariate matrix of dimension n x q, where q is the number of
#'   covariates. \code{NULL} if no covariate. Factor covariates must be supplied
#'   after transformation to dummy coding. No intercept must be supplied.
#' @param list_hyper An object of class "\code{hyper}" containing the model
#'   hyperparameters. Must be filled using the \code{\link{feed_hyperparam}}
#'   function or must be \code{NULL} for default hyperparameters.
#' @param list_init An object of class "\code{init}" containing the initial
#'   variational parameters. Must be filled using the
#'   \code{\link{feed_init_param}} function or be \code{NULL} for a default
#'   initialization.
#' @param list_blocks An object of class "\code{blocks}" containing setting for
#'   parallel inference on a partitioned predictor space. Must be filled using
#'   the \code{\link{set_blocks}} function or must be \code{NULL} for no
#'   partitioning.
#' @param user_seed Seed set for reproducible default choices of hyperparameters
#'   (if \code{list_hyper} is \code{NULL}) and inital variational parameters (if
#'   \code{list_init} is \code{NULL}). Also used at the cross-validation stage
#'   (if \code{list_cv} is non-\code{NULL}). Default is \code{NULL}, no seed set.
#' @param tol Tolerance for the stopping criterion.
#' @param maxit Maximum number of iterations allowed.
#' @param batch If \code{TRUE} a fast batch updating scheme is used (recommended).
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
#'  \item{ind_perm}{Vector of size n containing the permuted response sample
#'                  labels.}
#'  \item{gam_vb}{Posterior inclusion probability matrix of dimension p x d.
#'                Entry (s, t) corresponds to the posterior probability of
#'                association between predictor s and response t.}
#'  \item{om_vb}{Vector of size p containing the posterior mean of omega. Entry
#'              s controls the proportion of responses associated with predictor
#'              s.}
#
#' @seealso \code{\link{locus}}
#'
#' @examples
#'
#' user_seed <- 123; set.seed(user_seed)
#' n <- 200; p <- 250; p0 <- 20; d <- 25; d0 <- 20
#' list_X <- generate_snps(n = n, p = p)
#' list_Y <- generate_phenos(n = n, d = d, var_err = 0.25)
#'
#' dat <- generate_dependence(list_X, list_Y, ind_d0 = sample(1:d, d0),
#'                            ind_p0 = sample(1:p, p0), vec_prob_sh = 0.05,
#'                            max_tot_pve = 0.5)
#'
#' list_perm <- generate_null(n_perm = 10, dat$phenos, dat$snps, p0_av = p0,
#'                            user_seed = user_seed, verbose = FALSE)
#'
#' @export
#'
generate_null <- function(n_perm, Y, X, p0_av, Z = NULL, list_hyper = NULL,
                          list_init = NULL, list_blocks = NULL, user_seed = NULL,
                          tol = 1e-3, maxit = 1000, batch = TRUE, verbose = TRUE,
                          results_dir = NULL, n_cpus = 1) {

  if (!is.null(user_seed)){
    RNGkind("L'Ecuyer-CMRG") # ensure reproducibility when using mclapply
    set.seed(user_seed)
    if (verbose) cat(paste("Seed set to ", user_seed, ". \n", sep = ""))
  }

  n <- nrow(Y)

  permute <- function(i) {
    ind_perm <- sample(1:n) # random permutation
    rownames(Y) <- NULL

    # user_seed must be NULL here otherwise always the same permutation
    res_perm <- locus(Y = Y[ind_perm, ], X = X, p0_av = p0_av, Z = Z,
                      list_hyper = list_hyper, list_init = list_init,
                      list_cv = NULL, list_blocks = list_blocks, user_seed = NULL,
                      tol = tol, maxit = maxit, batch = batch, verbose = verbose)

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
