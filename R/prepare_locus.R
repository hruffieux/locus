# This file is part of the `locus` R package:
#     https://github.com/hruffieux/locus
#

# Internal function implementing sanity checks and needed preprocessing before
# the application of the different `locus_*_core` algorithms.
#
prepare_data_ <- function(Y, X, Z, V, link, ind_bin, user_seed, tol, maxit, verbose) {

  stopifnot(link %in% c("identity", "logit", "probit", "mix"))

  check_structure_(user_seed, "vector", "numeric", 1, null_ok = TRUE)

  check_structure_(tol, "vector", "numeric", 1)
  check_positive_(tol, eps=.Machine$double.eps)

  check_structure_(maxit, "vector", "numeric", 1)
  check_natural_(maxit)

  check_structure_(verbose, "vector", "logical", 1)

  check_structure_(X, "matrix", "numeric")

  n <- nrow(X)
  p <- ncol(X)

  check_structure_(Y, "matrix", "numeric")
  d <- ncol(Y)

  if (link == "mix") {

    ind_bin <- prepare_ind_bin_(d, ind_bin, link)

    if(!all(as.vector(Y[, ind_bin]) == as.numeric(as.logical(Y[, ind_bin]))))
      stop("The responses in Y corresponding to indices ind_bin must be a binary.")

  } else if (link != "identity"){

    if(!all(as.vector(Y) == as.numeric(as.logical(Y))))
      stop("Y must be a binary matrix for logistic/probit regression.")

  }

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
                           "levels. \n", sep=""))

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
  bool_cst_x <- list_X_cst$bool_cst
  rmvd_cst_x <- list_X_cst$rmvd_cst

  list_X_coll <- rm_collinear_(X, verbose)
  X <- list_X_coll$mat

  bool_coll_x <- list_X_coll$bool_coll
  rmvd_coll_x <- list_X_coll$rmvd_coll

  bool_rmvd_x <- bool_cst_x
  bool_rmvd_x[!bool_cst_x] <- bool_coll_x

  if (!is.null(V)) {

    check_structure_(V, "matrix", "numeric")

    r <- ncol(V)
    if (nrow(V) != p) stop("The number of rows of V must match the number of candidate predictors in X.")

    V <- V[!bool_rmvd_x, , drop = FALSE] # remove the rows corresponding to the removed candidate predictors

    if (is.null(rownames(V))) rownames(V) <- colnames(X)
    else if(any(rownames(V) != colnames(X)))
      stop("The provided rownames of Z must be the same than those of X and Y or NULL.")

    if (is.null(colnames(V))) colnames(V) <- paste("Annot_z_", 1:r, sep="")

    V <- scale(V)

    list_V_cst <- rm_constant_(V, verbose)
    V <- list_V_cst$mat
    bool_cst_v <- list_V_cst$bool_cst
    rmvd_cst_v <- list_V_cst$rmvd_cst

    list_V_coll <- rm_collinear_(V, verbose)
    V <- list_V_coll$mat
    r <- ncol(V)
    bool_coll_v <- list_V_coll$bool_coll
    rmvd_coll_v <- list_V_coll$rmvd_coll

    bool_rmvd_v <- bool_cst_v
    bool_rmvd_v[!bool_cst_v] <- bool_coll_v

    if (sum(!bool_rmvd_v) == 0)
      stop("All variables provided in V are constants and hence useless. Please set V to NULL.")

  } else {

    r <- NULL
    bool_rmvd_v <- NULL
    rmvd_cst_v <- NULL
    rmvd_coll_v <- NULL

  }

  p <- ncol(X)

  if (link == "identity") {

    Y <- scale(Y, center = TRUE, scale = FALSE)

  } else if (link == "mix") {

    Y[, -ind_bin] <- scale(Y[, -ind_bin], center = TRUE, scale = FALSE)

  } else if (link == "logit") {

    Y <- Y - 1 / 2

  }

  if (p < 1) stop(paste("There must be at least 1 non-constant candidate predictor ",
                        " stored in X.", sep=""))
  if (is.null(q) || q < 1) Z <- NULL
  if (is.null(r) || r < 1) V <- NULL # in principle useless given the above assert.

  create_named_list_(Y, X, Z, V,
                     bool_rmvd_x, bool_rmvd_z, bool_rmvd_v,
                     rmvd_cst_x, rmvd_cst_z, rmvd_cst_v,
                     rmvd_coll_x, rmvd_coll_z, rmvd_coll_v)

}


# Internal function implementing sanity checks and needed preprocessing for
# argument p0_av before the application of the different `locus_*_core` algorithms.
#
convert_p0_av_ <- function(p0_av, p, list_blocks, verbose, eps = .Machine$double.eps^0.5) {

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
                 "Please increase p0_av.",
                 sep = ""))

    if (p0_av / p > 0.95)
      stop(paste("p0_av = ", p0_av, ": \n",
                 "invalid provided value of p0_av.\n",
                 "Induces a non-sparse formulation. Please decrease p0_av.",
                 sep = ""))

    if (p0_av > ceiling(4 * p / 5))
      warning(paste("Prior model size p0_av = ", p0_av, ": \n",
                    "p0_av / p is large, so multiplicity control may be weak. ",
                    "You may want to consider a smaller p0_av.", sep=""))

    p_star <- p0_av

  } else {

    if (verbose) cat(paste("The sth entry of the provided p0_av ",
                           "interpreted as the prior probability that ",
                           "predictor s is associated with at least one ",
                           "response. \n\n",
                           sep = ""))

    if (any(p0_av < eps) | any(p0_av > 1 - eps))
      stop(paste("Invalid provided vector of p0_av.\n",
                 "All entries must lie between 0 and 1 (strictly).",
                 sep = ""))

    if (median(p0_av) > 1 / 2)
      warning(paste("The number of predictors with large prior inclusion ",
                    "probability is large, so multiplicity control may be weak. \n",
                    "You may want to decrease the values of several ",
                    "entries of p0_av.",
                    sep=""))

    p_star <- p0_av * p

  }

  # the sparsity level needs to be adapted when block-wise inference is used
  # otherwise the selected models may be too small (empirical considerations here)
  if (!is.null(list_blocks)) {
    p_star <- sapply(p_star, function(p_star_j) min(p_star_j * list_blocks$n_bl, 0.975 * p))

    if (verbose) cat(paste("The sparsity level is adapted for block-wise inference ",
                           "to ensure only sufficiently large models are selected.\n\n", sep = ""))
  }

  p_star
}


# Internal function implementing sanity checks and needed preprocessing for the
# model hyperparameters before the application of the different `locus_*_core`
# algorithms.
#
prepare_list_hyper_ <- function(list_hyper, Y, p, p_star, q, r, link, ind_bin,
                                bool_rmvd_x, bool_rmvd_z, bool_rmvd_v, names_x,
                                names_y, names_z, verbose) {

  d <- ncol(Y)

  if (is.null(list_hyper)) {

    if (verbose) cat("list_hyper set automatically. \n")

    list_hyper <- auto_set_hyper_(Y, p, p_star, q, r, link, ind_bin)

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

    if (list_hyper$link_hyper != link)
      stop(paste("The argument link is not consistent with the variable
                 link_hyper in list_hyper", sep=""))

    if(link == "mix") {
      if (!all(list_hyper$ind_bin_hyper == ind_bin))
        stop(paste("The argument ind_bin is not consistent with the variable
                   ind_bin_hyper in list_hyper", sep=""))
    }

    if (is.null(r)) {

      if (!is.null(list_hyper$r_hyper))
        stop(paste("The dimension (r) of the provided hyperparameters ",
                   "(list_hyper) is not consistent is V being NULL.\n", sep=""))

      if (inherits(list_hyper, "hyper")) {
        # remove the entries corresponding to the removed constant predictors in X
        # (if any)
        list_hyper$a <- list_hyper$a[!bool_rmvd_x]
        list_hyper$b <- list_hyper$b[!bool_rmvd_x]
      }

      if (!is.null(names(list_hyper$a)) && names(list_hyper$a) != names_x)
        stop("Provided names for the entries of a do not match the colnames of X.")

      if (!is.null(names(list_hyper$b)) && names(list_hyper$b) != names_x)
        stop("Provided names for the entries of b do not match the colnames of X.")

    } else {

      if (inherits(list_hyper, "hyper")) {
        # remove the entries corresponding to the removed constant predictors in X
        # (if any)
        r_hyper_match <- length(bool_rmvd_v)
        list_hyper$m0 <- list_hyper$m0[!bool_rmvd_x]
      } else {
        r_hyper_match <- r
      }

      if (list_hyper$r_hyper != r_hyper_match)
        stop(paste("The dimensions of the provided hyperparameters ",
                   "(list_hyper) are not consistent with that of V.", sep=""))

      if (!is.null(names(list_hyper$m0)) && names(list_hyper$m0) != names_x)
        stop("Provided names for the entries of m0 do not match the colnames of X.")

    }

    if (link %in% c("identity", "mix")) {

      if (link == "mix") names_y <- names_y[-ind_bin]

      if (!is.null(names(list_hyper$eta)) && names(list_hyper$eta) != names_y)
        stop("Provided names for the entries of eta do not match the colnames of the continuous variables in Y")

      if (!is.null(names(list_hyper$kappa)) && names(list_hyper$kappa) != names_y)
        stop("Provided names for the entries of kappa do not match the colnames of the continuous variables in Y")
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
                   "(list_hyper) are not consistent with that of Z.", sep=""))

      if (!is.null(names(list_hyper$phi)) && names(list_hyper$phi) != names_z)
        stop("Provided names for the entries of phi do not match the colnames of Z.")

      if (!is.null(names(list_hyper$xi)) && names(list_hyper$xi) != names_z)
        stop("Provided names for the entries of xi do not match the colnames of Z.")

    }

  }

  class(list_hyper) <- "out_hyper"

  list_hyper
}


# Internal function implementing sanity checks and needed preprocessing for the
# starting values before the application of the different `locus_*_core`
# algorithms.
#
prepare_list_init_ <- function(list_init, Y, p, p_star, q, r, link, ind_bin,
                               bool_rmvd_x, bool_rmvd_z, bool_rmvd_v, user_seed,
                               verbose) {

  d <- ncol(Y)
  n <- nrow(Y)

  if (is.null(list_init)) {

    if (!is.null(user_seed) & verbose) cat(paste("Seed set to user_seed ",
                                                 user_seed,". \n", sep=""))

    if (verbose) cat(paste("list_init set automatically. \n", sep=""))

    list_init <- auto_set_init_(Y, p, p_star, q, r, user_seed, link, ind_bin)

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

    if (list_init$link_init != link)
      stop(paste("The argument link is not consistent with the variable
                 link_init in list_init", sep=""))

    if(link == "mix") {
      if (!all(list_init$ind_bin_init == ind_bin))
        stop(paste("The argument ind_bin is not consistent with the variable
                   ind_bin_init in list_init", sep=""))
    }

    if (inherits(list_init, "init")) {
      # remove the entries corresponding to the removed constant predictors in X
      # (if any)
      list_init$gam_vb <- list_init$gam_vb[!bool_rmvd_x,, drop = FALSE]
      list_init$mu_beta_vb <- list_init$mu_beta_vb[!bool_rmvd_x,, drop = FALSE]

      if (link == "logit")
        list_init$sig2_beta_vb <- list_init$sig2_beta_vb[!bool_rmvd_x,, drop = FALSE]

    }

    if (!is.null(q)) {

      if (inherits(list_init, "init")) {
        q_init_match <- length(bool_rmvd_z)
        # remove the entries corresponding to the removed constant predictors in X
        # (if any)
        list_init$mu_alpha_vb <- list_init$mu_alpha_vb[!bool_rmvd_z,, drop = FALSE]

        if (link == "probit"){

          list_init$sig2_alpha_vb <- list_init$sig2_alpha_vb[!bool_rmvd_z]

        } else {

          list_init$sig2_alpha_vb <- list_init$sig2_alpha_vb[!bool_rmvd_z,, drop = FALSE]

        }

      } else {
        q_init_match <- q
      }

      if (list_init$q_init != q_init_match)
        stop(paste("The dimensions of the provided initial parameters ",
                   "(list_init) are not consistent with that of Z.", sep=""))
    }

    if (!is.null(r)) {

      if (inherits(list_init, "init")) {
        r_init_match <- length(bool_rmvd_v)
      } else {
        r_init_match <- r
      }

      if (list_init$r_init != r_init_match)
        stop(paste("The dimensions of the provided initial parameters ",
                   "(list_init) are not consistent with that of V.", sep=""))
    } else {

      if (!is.null(list_init$r_init))
        stop(paste("The dimension (r) of the provided initial parameters ",
                   "(list_init) is not consistent is V being NULL.\n", sep=""))
    }

  }

  class(list_init) <- "out_init"

  list_init
}


# Internal function implementing sanity checks and needed preprocessing for the
# model hyperparameters before the application of the cross-validation procedure
# for parameter p0_av.
#
prepare_cv_ <- function(list_cv, n, p, r, bool_rmvd_x, p0_av, link, list_hyper,
                        list_init, verbose) {

  if (!inherits(list_cv, "cv"))
    stop(paste("The provided list_cv must be an object of class ``cv''. \n",
               "*** you must either use the function set_cv to give the settings ",
               "for the cross-validation or set list_cv to NULL to skip the ",
               "cross-validation step. ***",
               sep=""))

  if (link == "logit")
    stop("Cross-validation not implemented only for logistic regression. Please, set list_cv to NULL or use probit regression.")

  if (!is.null(r))
    stop("Cross-validation implemented only models with no external information (V set to NULL). Please, set list_cv to NULL.")

  if (!is.null(p0_av) | !is.null(list_hyper) | !is.null(list_init))
    stop(paste("p0_av, list_hyper and list_init must all be NULL if non NULL ",
               "list_cv is provided (cross-validation).", sep = ""))

  if (list_cv$n_cv != n)
    stop(paste("The number of observations n provided to the function set_cv",
               "is not consistent with those of the data.", sep=""))

  if (list_cv$p_cv != length(bool_rmvd_x))
    stop(paste("The number of candidate predictor p provided to the function set_cv ",
               "is not consistent with X.", sep=""))

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


# Internal function implementing sanity checks and needed preprocessing to the
# settings provided by the user for block-wise parallel inference.
#
prepare_blocks_ <- function(list_blocks, r, bool_rmvd_x, list_cv) {

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

  tab_bl <- table(vec_fac_bl)
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
#' seed <- 123; set.seed(seed)
#'
#' ###################
#' ## Simulate data ##
#' ###################
#'
#' ## Example using small problem sizes:
#' ##
#' n <- 200; p <- 1200; p0 <- 200; d <- 50; d0 <- 40
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
#' n_bl <- 6
#' pos_bl <- seq(1, p, by = ceiling(p/n_bl))
#' list_blocks <- set_blocks(p, pos_bl, n_cpus = 2)
#'
#' vb <- locus(Y = Y, X = X, p0_av = p0, link = "identity",
#'             list_blocks = list_blocks, user_seed = seed)
#'
#' @seealso \code{\link{locus}}
#'
#' @export
set_blocks <- function(p, pos_bl, n_cpus, verbose = TRUE) {

  check_structure_(verbose, "vector", "logical", 1)

  check_structure_(pos_bl, "vector", "numeric")
  check_natural_(pos_bl)

  if (length(pos_bl) > 25)
    warning(paste("The provided number of blocks may be too large for accurate ",
                  "inference. If possible, use less blocks.", sep = ""))

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
                           "Please make sure that enough RAM is available. \n", sep=""))
  }

  p_blocks <- p

  list_blocks <- create_named_list_(p_blocks, n_bl, n_cpus, vec_fac_bl)

  class(list_blocks) <- "blocks"

  list_blocks
}


# Internal function implementing sanity checks the index of binary responses in
# case `locus_mix_core` or `locus_mix_info_core` is used.
#
prepare_ind_bin_ <- function(d, ind_bin, link) {

  if (link == "mix") {

    check_structure_(ind_bin, "vector", "numeric")
    ind_bin <- sort(unique(ind_bin))
    if (!all(ind_bin %in% 1:d))
      stop(paste("All indices provided in ind_bin must be integers between 1 ",
                 "and the total number of responses, d = ", d, ".", sep = ""))

    if (length(ind_bin) == d)
      stop(paste("Argument ind_bin indicates that all responses are binary. \n",
                 "Please set link to logit or probit, or change ind_bin to ",
                 "the indices of the binary responses only.", sep = ""))

  } else if (!is.null(ind_bin)) {

    stop("Argument ind_bin must be NULL if link is not set to mix.")

  }

  ind_bin
}
