prepare_data_ <- function(Y, X, Z, family, user_seed, tol, maxit, batch, verbose) {

  stopifnot(family %in% c("gaussian", "binomial"))

  check_structure_(user_seed, "vector", "numeric", 1, null_ok = TRUE)

  check_structure_(tol, "vector", "numeric", 1)
  check_positive_(tol, eps=.Machine$double.eps)

  check_structure_(maxit, "vector", "numeric", 1)
  check_natural_(maxit)

  check_structure_(batch, "vector", "logical", 1)

  check_structure_(verbose, "vector", "logical", 1)
#
  check_structure_(X, "matrix", "numeric")
  if (family == "gaussian") check_structure_(Y, "matrix", "double")
  else { check_structure_(Y, "matrix", "numeric")
    if(!identical(as.vector(Y), as.numeric(as.logical(Y))))
      stop("Y must be a binary matrix for logistic regression.")
  }
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

    if (family == "binomial")
      stop("Logistic regression with covariates Z not implemented yet.")

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

  if (family == "gaussian") Y <- scale(Y, center = TRUE, scale = FALSE)
  else Y <- Y - 1 / 2

  if (p < 1) stop(paste("There must be at least 1 non-constant candidate predictor ",
                        " stored in X.", sep=""))
  if (is.null(q) || q < 1) Z <- NULL

  create_named_list_(Y, X, Z, bool_rmvd_x, bool_rmvd_z,
                     rmvd_cst_x, rmvd_cst_z, rmvd_coll_x, rmvd_coll_z)

}


convert_p0_av_ <- function(p0_av, p, verbose, eps = .Machine$double.eps^0.5) {

  check_structure_(p0_av, "vector", "numeric", c(1, p))

  if (length(p0_av) == 1) {

    if (verbose) cat(paste("Provided p0_av = ", p0_av, " interpreted as ",
                           "the prior number of predictors associated with at ",
                           "least one response. \n\n", sep = ""))

    if (p0_av / p < eps)
      stop(paste("p0_av = ", p0_av, ": \n",
                 "invalid provided value of p0_av.\n",
                 "The prior sparsity level, p0_av / p, must be larger than ",
                 "zero. \n",
                 "Please increase p0_av. \n",
                 sep = ""))

    if (p0_av / p > 0.95)
      stop(paste("p0_av = ", p0_av, ": \n",
                 "invalid provided value of p0_av.\n",
                 "Induces a non-sparse formulation. Please decrease p0_av. \n",
                 sep = ""))

    if (p0_av > ceiling(4 * p / 5))
      warning(paste("Prior model size p0_av = ", p0_av, ": \n",
                    "p0_av / p is large, so multiplicity control may be weak. ",
                    "You may want to consider a smaller p0_av. \n", sep=""))

    p_star <- p0_av

  } else {

    if (verbose) cat(paste("- The sth entry of the provided p0_av ",
                           "interpreted as the prior probability that ",
                           "predictor s is associated with at least one ",
                           "response. \n\n",
                           sep = ""))

    if (any(p0_av) < eps | any(p0_av) > 1 - eps)
      stop(paste("Invalid provided vector of p0_av.\n",
                 "All entries must lie between 0 and 1 (strictly). \n",
                 sep = ""))

    if (median(p0_av) > 1 / 2)
      warning(paste("The number of predictors with large prior inclusion ",
                    "probability is large, so multiplicity control may be weak. \n",
                    "You may want to decrease the values of several ",
                    "entries of p0_av. \n",
                    sep=""))

    p_star <- p0_av * p

  }

  p_star
}


prepare_list_hyper_ <- function(list_hyper, Y, p, p_star, q, family,
                                bool_rmvd_x, bool_rmvd_z,
                                names_x, names_y, names_z, verbose) {

  d <- ncol(Y)

  if (is.null(list_hyper)) {

    if (verbose) cat("- list_hyper set automatically. \n")

    list_hyper <- auto_set_hyper_(Y, p, p_star, q, family)

  } else {

    if (!inherits(list_hyper, c("hyper", "out_hyper")))
      stop(paste("The provided list_hyper must be an object of class ``hyper'' ",
                 "or ``out_hyper''. \n",
                 "*** you must either use the function set_hyper to ",
                 "set your own hyperparameters or use list_hyper from a ``vb'' ",
                 "object or set the argument list_hyper to NULL for automatic choice. ***",
                 sep=""))

    if (inherits(list_hyper, "hyper")) {
      p_hyper_match <- length(bool_rmvd_x)
    } else {
      p_hyper_match <- p
    }


    if (list_hyper$d_hyper != d)
      stop(paste("The dimensions (d) of the provided hyperparameters ",
                 "(list_hyper) are not consistent with that of Y.\n", sep=""))

    if (list_hyper$p_hyper != p_hyper_match)
      stop(paste("The dimensions (p) of the provided hyperparameters ",
                 "(list_hyper) are not consistent with that of X.\n", sep=""))

    if (list_hyper$family_hyper != family)
      stop(paste("The argument family is not consistent with the variable
                 family_hyper in list_hyper", sep=""))

    if (inherits(list_hyper, "hyper")) {
      # remove the entries corresponding to the removed constant predictors in X
      # (if any)
      list_hyper$a <- list_hyper$a[!bool_rmvd_x]
      list_hyper$b <- list_hyper$b[!bool_rmvd_x]
    }

    if (!is.null(names(list_hyper$a)) && names(list_hyper$a) != names_x)
      stop("Provided names for the entries of a do not match the colnames of X")

    if (!is.null(names(list_hyper$b)) && names(list_hyper$b) != names_x)
      stop("Provided names for the entries of b do not match the colnames of X")

    if (family == "gaussian") {
      if (!is.null(names(list_hyper$eta)) && names(list_hyper$eta) != names_y)
        stop("Provided names for the entries of eta do not match the colnames of Y")

      if (!is.null(names(list_hyper$kappa)) && names(list_hyper$kappa) != names_y)
        stop("Provided names for the entries of kappa do not match the colnames of Y")
    }


    if (!is.null(q)) {

      if (inherits(list_hyper, "hyper")) {
        q_hyper_match <- length(bool_rmvd_z)
        # remove the entries corresponding to the removed constant predictors in X
        # (if any)
        list_hyper$phi <- list_hyper$phi[!bool_rmvd_z]
        list_hyper$xi <- list_hyper$xi[!bool_rmvd_z]
      } else {
        q_hyper_match <- q
      }

      if (list_hyper$q_hyper != q_hyper_match)
        stop(paste("The dimensions of the provided hyperparameters ",
                   "(list_hyper) are not consistent with that of Z.\n", sep=""))

      if (!is.null(names(list_hyper$phi)) && names(list_hyper$phi) != names_z)
        stop("Provided names for the entries of phi do not match the colnames of Z")

      if (!is.null(names(list_hyper$xi)) && names(list_hyper$xi) != names_z)
        stop("Provided names for the entries of xi do not match the colnames of Z")


    }

  }

  class(list_hyper) <- "out_hyper"

  list_hyper
}



prepare_list_init_ <- function(list_init, Y, p, p_star, q, family,
                               bool_rmvd_x, bool_rmvd_z, user_seed, verbose) {

  d <- ncol(Y)
  n <- nrow(Y)

  if (is.null(list_init)) {

    if (!is.null(user_seed) & verbose) cat(paste("- Seed set to user_seed ",
                                                 user_seed,". \n", sep=""))

    if (verbose) cat(paste("list_init set automatically. \n", sep=""))

    list_init <- auto_set_init_(Y, p, p_star, q, user_seed, family)

  } else {

    if (!is.null(user_seed))
      warning("user_seed not used since a non-NULL list_init was provided. \n")

    if (!inherits(list_init, c("init", "out_init")))
      stop(paste("The provided list_init must be an object of class ``init'' or ",
                 " `` out_init''. \n",
                 "*** you must either use the function set_init to ",
                 "set your own initialization or use list_init from a ``vb'' ",
                 "object or  set the argument list_init to NULL for automatic ",
                 "initialization. ***",
                 sep=""))

    if (inherits(list_init, "init")) {
      p_init_match <- length(bool_rmvd_x)
    } else {
      p_init_match <- p
    }


    if (list_init$d_init != d)
      stop(paste("The dimensions (d) of the provided initial parameters ",
                 "(list_init) are not consistent with that of Y.\n", sep=""))

    if (list_init$p_init != p_init_match)
      stop(paste("The dimensions (p) of the provided initial parameters ",
                 "(list_init) are not consistent with that of X.\n", sep=""))

    if (family == "binomial" && list_init$n_init != n)
      stop(paste("The number of observations provided when setting the initial ",
                 "parameters (list_init) is not consistent with that of X.\n", sep=""))

    if (list_init$family_init != family)
      stop(paste("The argument family is not consistent with the variable
                 family_init in list_init", sep=""))

    if (inherits(list_init, "init")) {
      # remove the entries corresponding to the removed constant predictors in X
      # (if any)
      list_init$gam_vb <- list_init$gam_vb[!bool_rmvd_x,, drop = FALSE]
      list_init$mu_beta_vb <- list_init$mu_beta_vb[!bool_rmvd_x,, drop = FALSE]

      if (family == "binomial")
        list_init$sig2_beta_vb <- list_init$sig2_beta_vb[!bool_rmvd_x,, drop = FALSE]

    }

    if (!is.null(q)) {

      if (inherits(list_init, "init")) {
        q_init_match <- length(bool_rmvd_z)
        # remove the entries corresponding to the removed constant predictors in X
        # (if any)
        list_init$mu_alpha_vb <- list_init$mu_alpha_vb[!bool_rmvd_z,, drop = FALSE]
        list_init$sig2_alpha_vb <- list_init$sig2_alpha_vb[!bool_rmvd_z,, drop = FALSE]
      } else {
        q_init_match <- q
      }

      if (list_init$q_init != q_init_match)
        stop(paste("The dimensions of the provided initial parameters ",
                   "(list_init) are not consistent with that of Z.\n", sep=""))
    }

  }

  class(list_init) <- "out_init"

  list_init
}


prepare_cv_ <- function(list_cv, n, p, bool_rmvd_x, p0_av, list_hyper, list_init,
                        verbose) {

  if (!inherits(list_cv, "cv"))
    stop(paste("The provided list_cv must be an object of class ``cv''. \n",
               "*** you must either use the function set_cv to give the settings ",
               "for the cross-validation or set list_cv to NULL to skip the ",
               "cross-validation step. ***",
               sep=""))

  if (!is.null(p0_av) | !is.null(list_hyper) | !is.null(list_init))
    stop(paste("p0_av, list_hyper and list_init must all be NULL if non NULL ",
               "list_cv is provided (cross-validation).", sep = ""))

  if (list_cv$n_cv != n)
    stop(paste("The number of observations n provided to the function set_cv",
               "is not consistent with those of the data.", sep=""))

  if (list_cv$p_cv != length(bool_rmvd_x))
    stop(paste("The number of candidate predictor p provided to the function set_cv ",
               "is not consistent with X.\n", sep=""))

  if (any(list_cv$p0_av_grid > p)) { # p has potentially been reduced because
    # of constant candidate predictors

    list_cv$p0_av_grid <- create_grid_(p, list_cv$size_p0_av_grid)

    new_size <- length(list_cv$p0_av_grid)
    if (list_cv$size_p0_av_grid > new_size) {
      if (verbose) cat(paste("Cross-validation p0_av_grid reduced to ", new_size,
                             " elements as p is small.\n", sep = ""))
      list_cv$size_p0_av_grid <- new_size
    }

    message <- paste("The cross-validation grid has been readjusted because to ",
                     "account for the removal of constant candidate predictors. Grid used: ",
                     list_cv$p0_av_grid, ". \n", sep = "")

    if (verbose) cat(message)
    else warning(message)
  }

  list_cv

}



prepare_blocks_ <- function(list_blocks, bool_rmvd_x, list_cv) {

  if (!inherits(list_blocks, "blocks"))
    stop(paste("The provided list_blocks must be an object of class ``blocks''. \n",
               "*** you must either use the function set_blocks to give the settings ",
               "for parallels applications of locus on blocks of candidate ",
               "predictors or set list_blocks to NULL to apply locus jointly on ",
               "all the candidate predictors (sufficient RAM required). ***",
               sep=""))

  if (!is.null(list_cv))
    stop(paste("list_cv must be NULL if non NULL ",
               "list_blocks is provided (parallel applications of locus on blocks of candidate predictors).\n",
               "Cross-validation for block-wise applications will be enabled soon.",sep = ""))

  if (list_blocks$p_blocks != length(bool_rmvd_x))
    stop(paste("The number of candidate predictors p provided to the function set_blocks ",
               "is not consistent with X.\n", sep=""))

  vec_fac_bl <- list_blocks$vec_fac_bl[!bool_rmvd_x]

  tab_bl <- table(list_blocks$vec_fac_bl)
  pres_bl <- tab_bl > 0

  # in case a block was removed due to the above because of bool_rmvd_x
  n_bl  <- sum(pres_bl)
  if(list_blocks$n_cpus > n_bl) n_cpus <- n_bl
  else n_cpus <- list_blocks$n_cpus

  create_named_list_(n_bl, n_cpus, vec_fac_bl)

}

#' Gather settings for parallel inference on partitioned predictor space.
#'
#' Parallel applications of the method on blocks of candidate predictors for
#' large datasets allows faster and less RAM-greedy executions.
#'
#' @param p Number of candidate predictors.
#' @param pos_bl Vector gathering the predictor block positions (first index of
#'   each block).
#' @param n_cpus Number of CPUs to be used. If large, one should ensure that
#'   enough RAM will be available for parallel execution. Set to 1 for serial
#'   execution.
#' @param verbose If \code{TRUE}, messages are displayed when calling
#'   \code{set_blocks}.
#'
#' @return An object of class "\code{blocks}" preparing the settings for parallel
#'   inference in a form that can be passed to the \code{\link{locus}}
#'   function.
#'
#' @examples
#' user_seed <- 123; set.seed(user_seed)
#' n <- 200; p <- 1200; p0 <- 200; d <- 50; d0 <- 40
#' list_X <- generate_snps(n = n, p = p)
#' list_Y <- generate_phenos(n = n, d = d, var_err = 0.25)
#'
#' dat <- generate_dependence(list_snps = list_X, list_phenos = list_Y,
#'                            ind_d0 = sample(1:d, d0), ind_p0 = sample(1:p, p0),
#'                            vec_prob_sh = 0.05, max_tot_pve = 0.9)
#' n_bl <- 6
#' pos_bl <- seq(1, p, by = ceiling(p/n_bl))
#' list_blocks <- set_blocks(p, pos_bl, n_cpus = 2)
#'
#' vb <- locus(Y = dat$phenos, X = dat$snps, p0_av = p0,
#'             family = "gaussian", list_blocks = list_blocks,
#'             user_seed = user_seed)
#'
#' @seealso \code{\link{locus}}
#'
#' @export
set_blocks <- function(p, pos_bl, n_cpus, verbose = TRUE) {

  check_structure_(verbose, "vector", "logical", 1)

  check_structure_(pos_bl, "vector", "numeric")
  check_natural_(pos_bl)


  if (any(pos_bl < 1) | any(pos_bl > p))
    stop("The positions provided in pos_bl must range between 1 and total number of variables in X, p.")

  if (any(duplicated(pos_bl)))
    stop("The positions provided in pos_bl must be unique.")

  if (any(pos_bl != cummax(pos_bl)))
    stop("The positions provided in pos_bl must be monotonically increasing.")

  vec_fac_bl <- as.factor(cumsum(seq_along(1:p) %in% pos_bl))

  n_bl <- length(unique(vec_fac_bl))

  check_structure_(n_cpus, "vector", "numeric", 1)
  check_natural_(n_cpus)


  if (n_cpus > 1) {

    n_cpus_avail <- parallel::detectCores()
    if (n_cpus > n_cpus_avail) {
      n_cpus <- n_cpus_avail
      warning(paste("The number of CPUs specified exceeds the number of CPUs ",
                    "available on the machine. The latter has been used instead.",
                    sep=""))
    }

    if (n_cpus > n_bl){
      message <- paste("The number of cpus in use is at most equal to the number of blocks.",
                       "n_cpus is therefore set to ", n_bl, ". \n", sep ="")
      if(verbose) cat(message)
      else warning(message)
      n_cpus <- n_bl
    }

    if (verbose) cat(paste("locus applied in parallel on ", n_bl,
                           " blocks of candidate predictors, using ", n_cpus, " CPUs.\n",
                           "Please make sure that enough RAM is available.", sep=""))
  }

  p_blocks <- p

  list_blocks <- create_named_list_(p_blocks, n_bl, n_cpus, vec_fac_bl)

  class(list_blocks) <- "blocks"

  list_blocks
}
