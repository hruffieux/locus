prepare_data_ <- function(Y, X, Z, user_seed, tol, maxit, batch, verbose) {

  check_structure_(user_seed, "vector", "numeric", 1, null_ok = T)

  check_structure_(tol, "vector", "numeric", 1)
  check_positive_(tol, eps=.Machine$double.eps)

  check_structure_(maxit, "vector", "numeric", 1)
  check_natural_(maxit)

  check_structure_(batch, "vector", "logical", 1)

  check_structure_(verbose, "vector", "logical", 1)

  check_structure_(X, "matrix", "numeric")
  check_structure_(Y, "matrix", "double")

  n <- nrow(X)
  p <- ncol(X)
  d <- ncol(Y)

  if (nrow(Y) != n) stop("X and Y must have the same number of observations.")

  if (is.null(rownames(X)) & is.null(rownames(Y)))
    rownames(X) <- rownames(Y) <- paste("Ind_", 1:n, sep="")
  else if (is.null(rownames(X))) rownames(X) <- rownames(Y)
  else if (is.null(rownames(Y))) rownames(Y) <- rownames(X)
  else if (any(rownames(X) != rownames(Y)))
    stop("The provided rownames of X and Y must be the same.")

  if (is.null(colnames(X))) colnames(X) <- paste("Cov_x_", 1:p, sep="")
  if (is.null(colnames(Y))) colnames(Y) <- paste("Resp_", 1:d, sep="")

  if (!is.null(Z)) {

    if (verbose) cat(paste("Each factor variables must be provided by adding ",
                           "to Z (nb levels - 1) variables representing their ",
                           "levels.", sep=""))

    check_structure_(Z, "matrix", "numeric")

    q <- ncol(Z)

    if (nrow(Z) != n) stop("Z must have the same number of observations as Y and X.")

    if (is.null(rownames(Z))) rownames(Z) <- rownames(X)
    else if(any(rownames(Z) != rownames(X)))
      stop("The provided rownames of Z must be the same than those of X and Y or NULL.")

    if (is.null(colnames(Z))) colnames(Z) <- paste("Cov_z_", 1:q, sep="")

    Z <- scale(Z)

    list_Z_cst <- rm_constant_(Z, verbose)
    Z <- list_Z_cst$mat
    q <- ncol(Z)
    bool_cst_z <- list_Z_cst$bool_cst
    rmvd_cst_z <- list_Z_cst$rmvd_cst

    list_Z_coll <- rm_collinear_(Z, verbose)
    Z <- list_Z_coll$mat
    q <- ncol(Z)
    bool_coll_z <- list_Z_coll$bool_coll
    rmvd_coll_z <- list_Z_coll$rmvd_coll

    bool_rmvd_z <- bool_cst_z
    bool_rmvd_z[!bool_cst_z] <- bool_coll_z

  } else {
    q <- NULL
    bool_rmvd_z <- NULL
    rmvd_cst_z <- NULL
    rmvd_coll_z <- NULL
  }


  X <- scale(X)

  list_X_cst <- rm_constant_(X, verbose)
  X <- list_X_cst$mat
  p <- ncol(X)
  bool_cst_x <- list_X_cst$bool_cst
  rmvd_cst_x <- list_X_cst$rmvd_cst

  list_X_coll <- rm_collinear_(X, verbose)
  X <- list_X_coll$mat
  p <- ncol(X)
  bool_coll_x <- list_X_coll$bool_coll
  rmvd_coll_x <- list_X_coll$rmvd_coll

  bool_rmvd_x <- bool_cst_x
  bool_rmvd_x[!bool_cst_x] <- bool_coll_x

  Y <- scale(Y, center = T, scale = F)

  if (p < 1) stop(paste("There must be at least 1 non-constant covariate ",
                        " stored in X.", sep=""))
  if (is.null(q) || q < 1) Z <- NULL

  create_named_list_(Y, X, Z, bool_rmvd_x, bool_rmvd_z,
                     rmvd_cst_x, rmvd_cst_z, rmvd_coll_x, rmvd_coll_z)

}


convert_p0_av_ <- function(p0_av, p, verbose, eps = .Machine$double.eps^0.5) {

  check_structure_(p0_av, "vector", "numeric", c(1, p))

  if (length(p0_av) == 1) {

    if (verbose) cat(paste("Provided p0_av = ", p0_av, " interpreted as ",
                           "the prior number of covariates associated with at ",
                           "least one response. \n\n", sep = ""))

    if (p0_av / p < eps)
      stop(paste("p0_av = ", p0_av, ": \n",
                 "invalid provided value of p0_av.\n",
                 "The prior sparsity level, p0_av / p, must be larger than ",
                 "zero. \n",
                 "Please increase p0_av. \n",
                 sep = ""))

    if (p0_av / p > 1 - eps)
      stop(paste("p0_av = ", p0_av, ": \n",
                 "invalid provided value of p0_av.\n",
                 "The prior sparsity level, p0_av / p, must be smaller than ",
                 "one. \n",
                 "Please decrease p0_av. \n",
                 sep = ""))

    if (p0_av > ceiling(p / 4))
      warning(paste("p0_av = ", p0_av, ": \n",
                    "lower prior sparsity levels, p0_av / p, are commonly ",
                    "assumed for genetic association studies. \n",
                    "You may want to consider a smaller value of p0_av. \n",
                    sep=""))

    p_star <- p0_av

  } else {

    if (verbose) cat(paste("- The sth entry of the provided p0_av ",
                           "interpreted as the prior probability that ",
                           "covariate s is associated with at least one ",
                           "response. \n\n",
                           sep = ""))

    if (any(p0_av) < eps | any(p0_av) > 1 - eps)
      stop(paste("Invalid provided vector of p0_av.\n",
                 "All entries must lie between 0 and 1 (strictly). \n",
                 sep = ""))

    if (median(p0_av) > 1 / 4)
      warning(paste("A lower number of covariates is commonly assumed to have ",
                    "a significant prior probability of association with ",
                    "with responses in genetic association studies. \n",
                    "You may want to decrease the value of several ",
                    "entries of p0_av. \n",
                    sep=""))

    p_star <- p0_av * p

  }

  p_star
}


prepare_list_hyper_ <- function(list_hyper, Y, d, p, p_star, q,
                                bool_rmvd_x, bool_rmvd_z,
                                names_x, names_y, names_z, verbose) {

  if (is.null(list_hyper)) {

    if (verbose) cat("- list_hyper set automatically. \n")

    list_hyper <- auto_set_hyperparam_(Y, p, p_star, q)

  } else {

    if (!(class(list_hyper) %in% c("hyper", "out_hyper")))
      stop(paste("The provided list_hyper must be an object of class ``hyper'' ",
                 "or ``out_hyper''. \n",
                 "*** you must either use the function feed_hyperparam to ",
                 "set your own hyperparameters or use list_hyper from a ``vb'' ",
                 "object or set the argument list_hyper to NULL for automatic choice. ***",
                 sep=""))

    if (class(list_hyper) == "hyper") {
      p_hyper_match <- length(bool_rmvd_x)
    } else {
      p_hyper_match <- p
    }


    if (list_hyper$d_hyper != d)
      stop(paste("The dimensions of the provided hyperparameters ",
                 "(list_hyper) are not consistent with those of Y.\n", sep=""))

    if (list_hyper$p_hyper != p_hyper_match)
      stop(paste("The dimensions of the provided hyperparameters ",
                 "(list_hyper) are not consistent with those of X.\n", sep=""))

    if (class(list_hyper) == "hyper") {
      # remove the entries corresponding to the removed constant covariates in X
      # (if any)
      list_hyper$a <- list_hyper$a[!bool_rmvd_x]
      list_hyper$b <- list_hyper$b[!bool_rmvd_x]
    }

    if (!is.null(names(list_hyper$a)) && names(list_hyper$a) != names_x)
      stop("Provided names for the entries of a do not match the colnames of X")

    if (!is.null(names(list_hyper$b)) && names(list_hyper$b) != names_x)
      stop("Provided names for the entries of b do not match the colnames of X")

    if (!is.null(names(list_hyper$eta)) && names(list_hyper$eta) != names_y)
      stop("Provided names for the entries of eta do not match the colnames of Y")

    if (!is.null(names(list_hyper$kappa)) && names(list_hyper$kappa) != names_y)
      stop("Provided names for the entries of kappa do not match the colnames of Y")



    if (!is.null(q)) {

      if (class(list_hyper) == "hyper") {
        q_hyper_match <- length(bool_rmvd_z)
        # remove the entries corresponding to the removed constant covariates in X
        # (if any)
        list_hyper$phi <- list_hyper$phi[!bool_rmvd_z]
        list_hyper$xi <- list_hyper$xi[!bool_rmvd_z]
      } else {
        q_hyper_match <- q
      }

      if (list_hyper$q_hyper != q_hyper_match)
        stop(paste("The dimensions of the provided hyperparameters ",
                   "(list_hyper) are not consistent with those of Z.\n", sep=""))

      if (!is.null(names(list_hyper$phi)) && names(list_hyper$phi) != names_z)
        stop("Provided names for the entries of phi do not match the colnames of Z")

      if (!is.null(names(list_hyper$xi)) && names(list_hyper$xi) != names_z)
        stop("Provided names for the entries of xi do not match the colnames of Z")


    }

  }

  list_hyper

}



prepare_list_init_ <- function(list_init, Y, d, p, p_star, q, bool_rmvd_x, bool_rmvd_z,
                               names_x, names_y, names_z, user_seed, verbose) {

  if (is.null(list_init)) {

    if (!is.null(user_seed) & verbose) cat(paste("- Seed set to user_seed ",
                                                 user_seed,". \n", sep=""))

    if (verbose) cat(paste("list_init set automatically. \n", sep=""))

    list_init <- auto_init_param_(Y, p, p_star, user_seed, q)

  } else {

    if (!is.null(user_seed))
      warning("user_seed not used since a non-NULL list_init was provided. \n")

    if (!(class(list_init) %in% c("init", "out_init")))
      stop(paste("The provided list_init must be an object of class ``init'' or ",
                 " `` out_init''. \n",
                 "*** you must either use the function feed_init_param to ",
                 "set your own initialization or use list_init from a ``vb'' ",
                 "object or  set the argument list_init to NULL for automatic ",
                 "initialization. ***",
                 sep=""))

    if (class(list_init) == "init") {
      p_init_match <- length(bool_rmvd_x)
    } else {
      p_init_match <- p
    }


    if (list_init$d_init != d)
      stop(paste("The dimensions of the provided initial parameters ",
                 "(list_init) are not consistent with those of Y.\n", sep=""))

    if (list_init$p_init != p_init_match)
      stop(paste("The dimensions of the provided initial parameters ",
                 "(list_init) are not consistent with those of X.\n", sep=""))

    if (class(list_init) == "init") {
      # remove the entries corresponding to the removed constant covariates in X
      # (if any)
      list_init$gam_vb <- list_init$gam_vb[!bool_rmvd_x,, drop=F]
      list_init$mu_beta_vb <- list_init$mu_beta_vb[!bool_rmvd_x,, drop=F]
      list_init$om_vb <- list_init$om_vb[!bool_rmvd_x]
    }

    if (!is.null(q)) {

      if (class(list_init) == "init") {
        q_init_match <- length(bool_rmvd_z)
        # remove the entries corresponding to the removed constant covariates in X
        # (if any)
        list_init$mu_alpha_vb <- list_init$mu_alpha_vb[!bool_rmvd_z,, drop=F]
        list_init$sig2_alpha_vb <- list_init$sig2_alpha_vb[!bool_rmvd_z,, drop=F]
        list_init$zeta2_inv_vb <- list_init$zeta2_inv_vb[!bool_rmvd_z]
      } else {
        q_init_match <- q
      }

      if (list_init$q_init != q_init_match)
        stop(paste("The dimensions of the provided initial parameters ",
                   "(list_init) are not consistent with those of Z.\n", sep=""))
    }

  }

  list_init
}


prepare_cv_ <- function(list_cv, n, p, bool_rmvd_x, p0_av, list_hyper, list_init, verbose) {

  if (class(list_cv) != "cv")
    stop(paste("The provided list_cv must be an object of class ``cv''. \n",
               "*** you must either use the function set_cv to give the settings ",
               "for the cross-validation or set list_cv to NULL to skip the ",
               "cross-validation step. ***",
               sep=""))

  if (!is.null(p0_av) | !is.null(list_hyper) | !is.null(list_init))
    stop(paste("p0_av, list_hyper and list_init must all be NULL if if non NULL ",
               "list_cv is provided (cross-validation).", sep = ""))

  if (list_cv$n_cv != n)
    stop(paste("The number of observations n provided to the function set_cv",
               "is not consistent with that of the data.", sep=""))

  if (list_cv$p_cv != length(bool_rmvd_x))
    stop(paste("The number of covariate p provided to the function set_cv ",
               "is not consistent with X.\n", sep=""))

  if (any(list_cv$p0_av_grid > p)) { # p has potentially been reduced because
                                       # of constant covariates

    list_cv$p0_av_grid <- create_grid_(p, list_cv$size_p0_av_grid)

    new_size <- length(list_cv$p0_av_grid)
    if (list_cv$size_p0_av_grid > new_size) {
      if (verbose) cat(paste("Cross-validation p0_av_grid reduced to ", new_size,
                             " elements as p is small.\n", sep = ""))
      list_cv$size_p0_av_grid <- new_size
    }

    message <- paste("The cross-validation grid has been readjusted because to ",
                     "account for the removal of constant covariates. Grid used: ",
                     list_cv$p0_av_grid, ". \n", sep = "")

    if (verbose) cat(message)
    else warning(message)
  }

  list_cv

}


