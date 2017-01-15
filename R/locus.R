#' Fit a sparse multivariate regression model using variational inference.
#'
#' Variational approximation procedure fitting a sparse multivariate regression
#' model for combined selection of predictors and associated responses in
#' high-dimensional set-ups. Dependence across responses linked to the same
#' predictors is captured through the model hierarchical structure.
#'
#' @param Y Continuous response data matrix of dimension n x d, where n is the
#'   number of observations and d is the number of response variables.
#' @param X Input matrix of dimension n x p, where p is the number of candidate
#'   predictors. \code{X} cannot contain NAs. No intercept must be supplied.
#' @param p0_av Prior average number of predictors expected to be included in
#'   the model. Must be \code{NULL} if \code{list_init} and \code{list_hyper}
#'   are both non-\code{NULL} or if \code{list_cv} is non-\code{NULL}.
#' @param Z Covariate matrix of dimension n x q, where q is the number of
#'   covariates. \code{NULL} if no covariate. Factor covariates must be supplied
#'   after transformation to dummy coding. No intercept must be supplied.
#' @param list_hyper List containing the model hyperparameters. Must be filled
#'   using the \code{\link{feed_hyperparam}} function or must be \code{NULL} for
#'   default hyperparameters.
#' @param list_init List containing the variational initial parameters. Must be
#'   be filled using the \code{\link{feed_init_param}} function or be
#'   \code{NULL} for a default initialization.
#' @param list_cv List containing settings for choosing the prior average number
#'   of predictors expected to be included in the model, \code{p0_av}, by
#'   cross-validation. Must be filled using the \code{\link{set_cv}} function or
#'   must be \code{NULL} for no cross-validation. If non-\code{NULL},
#'   \code{p0_av}, \code{list_init} and \code{list_hyper} must all be
#'   \code{NULL}.
#' @param list_blocks List containing setting for parallel inference on a
#'   partitioned predictor space. Must be filled using the
#'   \code{\link{set_blocks}} function or must be \code{NULL} for no
#'   partitioning.
#' @param user_seed Seed set for reproducible default choices of hyperparameters
#'   (if \code{list_hyper} is \code{NULL}) and inital variational parameters (if
#'   \code{list_init} is \code{NULL}). Also used at the cross-validation stage
#'   (if \code{list_cv} is non-\code{NULL}).
#' @param tol Tolerance for the stopping criterion.
#' @param maxit Maximum number of iterations allowed.
#' @param batch If TRUE a fast batch updating scheme is used (recommended).
#' @param save_hyper If TRUE, the hyperparameters used for the model are saved
#'   as output.
#' @param save_init If TRUE, the initial variational parameters used for the
#'   inference are saved as output.
#' @param verbose If TRUE, messages are displayed during execution.
#'
#' @return A list containing the following variational estimates and settings:
#'  \item{lb_opt}{Optimized variational lower bound for the marginal
#'                log-likelihood.}
#'  \item{gam_vb}{Posterior inclusion probability matrix of dimension p x d.
#'                Entry (s, t) corresponds to the posterior probability of
#'                association between predictor s and response t.}
#'  \item{mu_alpha_vb}{Matrix of dimension q x d whose entries are the posterior
#'                     mean regression coefficients for the covariates provided
#'                     in \code{Z}. \code{NULL} if \code{Z} is \code{NULL}.}
#'  \item{om_vb}{Vector of size p containing the posterior mean of omega. Entry
#'              s controls the proportion of responses associated with predictor
#'              s.}
#'  \item{x_prpnst}{Vector of size p containing the sums over the colums of
#'                  \code{gam_vb}, used to evaluate the propensity for the
#'                  candidate predictors to be involved in associations.}
#'  \item{y_prpnst}{Vector of size d containing the sums over the rows of
#'                  \code{gam_vb}, used to evaluate the propensity for the
#'                  responses to be involved in associations.}
#'  \item{rmvd_cst_x, rmvd_cst_z}{Vectors containing the indices of constant
#'                                variables in X (resp. Z) removed prior to the
#'                                analysis.}
#'  \item{rmvd_coll_x, rmvd_coll_z}{Vectors containing the indices of variables
#'                                  in X (resp. Z) removed prior to the
#'                                  analysis because collinear to other
#'                                  variables.}
#'  \item{list_hyper, list_init}{If \code{save_hyper}, resp. \code{save_init},
#'                               TRUE, hyperparameters, resp. initial
#'                               variational parameters, used for inference are
#'                               saved as output.}
#' @examples
#'
#' user_seed <- 123; set.seed(user_seed)
#' n <- 200; p <- 500; p0 <- 100; d <- 50; d0 <- 40
#' list_X <- generate_snps(n = n, p = p)
#' list_Y <- generate_phenos(n = n, d = d, var_err = 0.25)
#'
#' dat <- generate_dependence(list_snps = list_X, list_phenos = list_Y,
#'                            ind_d0 = sample(1:d, d0), ind_p0 = sample(1:p, p0),
#'                            vec_prob_sh = 0.1, pve_per_snp = 0.05)
#'
#' vb <- locus(Y = dat$phenos, X = dat$snps, p0_av = p0, user_seed = user_seed)
#'
#'
#' @seealso \code{\link{feed_hyperparam}}, \code{\link{feed_init_param}},
#'   \code{\link{set_cv}}, \code{\link{set_blocks}}
#'
#' @export
#'
locus <- function(Y, X, p0_av, Z = NULL,
                  list_hyper = NULL, list_init = NULL, list_cv = NULL,
                  list_blocks = NULL, user_seed = NULL, tol = 1e-4, maxit = 1000,
                  batch = T, save_hyper = F, save_init = F, verbose = T) { ##


  if (verbose) cat("== Preparing the data ... \n")
  dat <- prepare_data_(Y, X, Z, user_seed, tol, maxit, batch, verbose)

  bool_rmvd_x <- dat$bool_rmvd_x
  bool_rmvd_z <- dat$bool_rmvd_z

  X <- dat$X
  Y <- dat$Y
  Z <- dat$Z

  n <- nrow(X)
  p <- ncol(X)
  d <- ncol(Y)

  names_x <- colnames(X)
  names_y <- colnames(Y)
  if (!is.null(Z)) {
    q <- ncol(Z)
    names_z <- colnames(Z)
  } else {
    q <- NULL
    names_z <- NULL
  }

  if (verbose) cat("... done. == \n\n")



  if (!is.null(list_cv) & is.null(list_blocks)) { ## TODO: allow cross-validation when list_blocks is used.

    if (verbose){
      cat("=============================== \n")
      cat("===== Cross-validation... ===== \n")
      cat("=============================== \n")
    }
    list_cv <- prepare_cv_(list_cv, n, p, bool_rmvd_x, p0_av, list_hyper,
                           list_init, verbose)

    p_star <- cross_validate_(Y, X, Z, d, n, p, q, list_cv, user_seed, verbose)

  } else if (!is.null(list_blocks)) {

    list_blocks <- prepare_blocks_(list_blocks, bool_rmvd_x, p0_av, list_hyper,
                                   list_init, list_cv)

    p_star <- list_blocks$p_star
    n_bl <- list_blocks$n_bl
    n_cpus <- list_blocks$n_cpus
    vec_fac_bl <- list_blocks$vec_fac_bl
    vec_p_bl <- list_blocks$vec_p_bl


  } else {

    if (is.null(list_hyper) | is.null(list_init)) {

      p_star <- convert_p0_av_(p0_av, p, verbose)

      # remove the entries corresponding to the removed constant covariates in X (if any)
      if (length(p_star) > 1) p_star <- p_star[!bool_rmvd_x]

    } else {

      if (!is.null(p0_av))
        warning(paste("Provided argument p0_av not used, as both list_hyper ",
                      "and list_init were provided.", sep = ""))

      p_star <- NULL

    }

  }


  if (verbose) cat("== Preparing the hyperparameters ... \n")
  list_hyper <- prepare_list_hyper_(list_hyper, Y, d, p, p_star, q, bool_rmvd_x,
                                    bool_rmvd_z, names_x, names_y, names_z, verbose)
  if (verbose) cat("... done. == \n\n")

  if (verbose) cat("== Preparing the parameter initialization ... \n")
  list_init <- prepare_list_init_(list_init, Y, d, p, p_star, q, bool_rmvd_x,
                                  bool_rmvd_z, names_x, names_y, names_z,
                                  user_seed, verbose)

  if (verbose) cat("... done. == \n\n")

  if (verbose){
    cat(paste("========================================================== \n",
              "== Variational algorithm sparse multivariate regression == \n",
              "========================================================== \n\n",
              sep = ""))
  }

  if (is.null(list_blocks)) {

    if (is.null(q))
      vb <- locus_core_(Y, X, d, n, p, list_hyper, list_init$gam_vb,
                        list_init$mu_beta_vb,
                        list_init$sig2_beta_vb,
                        list_init$tau_vb, tol, maxit,
                        batch, verbose)
    else
      vb <- locus_z_core_(Y, X, Z, d, n, p, q, list_hyper,
                          list_init$gam_vb, list_init$mu_beta_vb,
                          list_init$sig2_beta_vb,
                          list_init$tau_vb,
                          list_init$mu_alpha_vb, list_init$sig2_alpha_vb,
                          list_init$zeta2_inv_vb, tol, maxit, batch, verbose)

  } else {


    list_pos_bl <- split(1:p, vec_fac_bl)

    split_bl_hyper <- lapply(list_pos_bl, function(pos_bl) {
      list_hyper$p_hyper <- length(n_bl)
      list_hyper$a <- list_hyper$a[pos_bl]
      list_hyper$b <- list_hyper$b[pos_bl]
      list_hyper
    })

    split_bl_init <- lapply(list_pos_bl, function(pos_bl) {
      list_init$p_init <- length(pos_bl)
      list_init$gam_vb <- list_init$gam_vb[pos_bl,, drop = F]
      list_init$mu_beta_vb <- list_init$mu_beta_vb[pos_bl,, drop = F]
      list_init
    })

    locus_bl_ <- function(k) {

      X_bl <- X[, list_pos_bl[[k]]]

      list_hyper_bl <- split_bl_hyper[[k]]
      list_init_bl <- split_bl_init[[k]]

      if (is.null(q))

        vb_bl <- locus_core_(Y, X_bl, d, n, vec_p_bl[k], list_hyper_bl,
                             list_init_bl$gam_vb, list_init_bl$mu_beta_vb,
                             list_init_bl$sig2_beta_vb,
                             list_init_bl$tau_vb,
                             tol, maxit, batch, verbose)
      else
        vb_bl <- locus_z_core_(Y, X_bl, Z, d, n, vec_p_bl[k], q, list_hyper_bl,
                               list_init_bl$gam_vb, list_init_bl$mu_beta_vb,
                               list_init_bl$sig2_beta_vb,
                               list_init_bl$tau_vb,
                               list_init_bl$mu_alpha_vb, list_init_bl$sig2_alpha_vb,
                               list_init_bl$zeta2_inv_vb, tol, maxit, batch, verbose)
      vb_bl
    }

    list_vb <- parallel::mclapply(1:n_bl, function(k) locus_bl_(k), mc.cores = n_cpus)

    names_vec <- c("lb_opt", "om_vb", "x_prpnst", "p_star")
    names_mat <- c("gam_vb", "y_prpnst")

    vb <- c(lapply(names_vec, function(key) do.call(c, lapply(list_vb, `[[`, key))),
            lapply(names_mat, function(key) do.call(rbind, lapply(list_vb, `[[`, key))))

    names(vb) <- c(names_vec, names_mat)
    vb$y_prpnst <- colSums(vb$y_prpnst)

  }

  vb$p_star <- p_star

  vb$rmvd_cst_x <- dat$rmvd_cst_x
  vb$rmvd_coll_x <- dat$rmvd_coll_x
  if (!is.null(Z)) {
    vb$rmvd_cst_z <- dat$rmvd_cst_z
    vb$rmvd_coll_z <- dat$rmvd_coll_z
  }

  if (save_hyper) vb$list_hyper <- list_hyper
  if (save_init) vb$list_init <- list_init

  class(vb) <- "locus"


  if (verbose) cat("... done. == \n\n")

  vb

}

