#' Fit a sparse multivariate regression model using variational inference.
#'
#' \code{locus} applies a variational approximation procedure to fit a sparse
#' multivariate regression model allowing combined selection of predictors and
#' associated responses in high-dimensional set-ups. Dependence across responses
#' linked to the same predictors is captured through the model hierarchical
#' structure.
#'
#' @param Y continuous response data matrix of dimension n x d, where n is the
#'   number of observations and d is the number of response variables.
#' @param X input matrix of dimension n x p, where p is the number of candidate
#'   predictors. X cannot contain NAs.
#' @param p0_av is the 'a priori' average number of predictors included in the
#'   model.
#' @param Z covariate matrix of dimension n x q, where q is the number of
#'   covariates. NULL if no covariate. Factor covariates must be supplied after
#'   transformation to dummy coding. \code{locus} standardizes the variables so
#'   no intercept must be supplied.
#'
#' @export
locus <- function(Y, X, p0_av, Z = NULL,
                   list_hyper = NULL, list_init = NULL, list_cv = NULL,
                   user_seed = NULL, tol = 1e-3, maxit = 1000, batch = T,
                   save_hyper = F, save_init = F, verbose = T) { ##


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



  if (!is.null(list_cv)) {

    if (verbose){
      cat("=============================== \n")
      cat("===== Cross-validation... ===== \n")
      cat("=============================== \n")
    }
    list_cv <- prepare_cv_(list_cv, n, p, bool_rmvd_x, p0_av, list_hyper,
                           list_init, verbose)

    p_star <- cross_validate_(Y, X, Z, d, n, p, q, list_cv, user_seed, verbose)

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
  list_hyper <- prepare_list_hyper_(list_hyper, Y, d, p, p_star, q,
                                    bool_rmvd_x, bool_rmvd_z,
                                    names_x, names_y, names_z, verbose)
  if (verbose) cat("... done. == \n\n")

  if (verbose) cat("== Preparing the parameter initialization ... \n")
  list_init <- prepare_list_init_(list_init, Y, d, p, p_star, q,
                                  bool_rmvd_x, bool_rmvd_z,
                                  names_x, names_y, names_z,
                                  user_seed, verbose)

  if (verbose) cat("... done. == \n\n")

  if (verbose){
    cat(paste("========================================================== \n",
              "== Variational algorithm sparse multivariate regression == \n",
              "========================================================== \n\n",
              sep = ""))
  }


  if (is.null(q))
    vb <- locus_core_(Y, X, d, n, p, list_hyper, list_init, tol, maxit, batch,
                       verbose, full_output = F)
  else
    vb <- locus_z_core_(Y, X, Z, d, n, p, q, list_hyper, list_init, tol, maxit,
                         batch, verbose, full_output = F)

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

  vb

}
