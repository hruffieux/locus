#' Gather settings for the cross-validation procedure used in \code{locus}.
#'
#' The cross-validation procedure uses the variational lower bound as objective
#' function and is used to select the prior average number of predictors
#' \code{p0_av} expected to be included in the model used to set the model
#' hyperparameters and ensure sparse predictor selections.
#'
#' @param n Number of observations.
#' @param p Number of candidate predictors.
#' @param n_folds Number of number of folds. Large folds are not recommended for
#'   large datasets as the procedure may become computationally expensive. Must
#'   be greater than 2 and smaller than the number of observations.
#' @param size_p0_av_grid Number of possible values of p0_av to be compared.
#'   Large numbers are not recommended for large datasets as the procedure may
#'   become computationally expensive.
#' @param n_cpus Number of CPUs to be used for the cross-validation procedure.
#'   If large, one should ensure that enough RAM will be available for parallel
#'   execution. Set to 1 for serial execution.
#' @param tol_cv Tolerance for the variational algorithm stopping criterion used
#'   within the cross-validation procedure.
#' @param maxit_cv Maximum number of iterations allowed for the variational
#'   algorithm used within the cross-validation procedure.
#' @param batch_cv If \code{TRUE}, a fast batch updating scheme is used within
#'   the cross-validation procedure (recommended).
#' @param verbose If \code{TRUE}, messages are displayed when calling
#'   \code{set_cv}.
#'
#' @return An object of class "\code{cv}" preparing the settings for the
#'   cross-validation settings in a form that can be passed to the
#'   \code{\link{locus}} function.
#'
#' @examples
#' user_seed <- 123; set.seed(user_seed)
#' n <- 150; p <- 200; p0 <- 50; d <- 25; d0 <- 20
#' list_X <- generate_snps(n = n, p = p)
#' list_Y <- generate_phenos(n = n, d = d, var_err = 0.25)
#'
#' dat <- generate_dependence(list_snps = list_X, list_phenos = list_Y,
#'                            ind_d0 = sample(1:d, d0),
#'                            ind_p0 = sample(1:p, p0), vec_prob_sh = 0.1,
#'                            max_tot_pve = 0.9)
#'
#' list_cv <- set_cv(n, p, n_folds = 3, size_p0_av_grid = 3, n_cpus = 2)
#'
#' vb <- locus(Y = dat$phenos, X = dat$snps, p0_av = NULL, family = "gaussian",
#'             list_cv = list_cv, user_seed = user_seed)
#'
#' @seealso \code{\link{locus}}
#'
#' @export
#'
set_cv <- function(n, p, n_folds, size_p0_av_grid, n_cpus, tol_cv = 1e-3,
                   maxit_cv = 1e3, batch_cv = TRUE, verbose = TRUE) {

  check_structure_(n_folds, "vector", "numeric", 1)
  check_natural_(n_folds)

  if (!(n_folds %in% 2:n))
    stop("n_folds must be a natural number greater than 2 and smaller than the number of observations.")

  # 16 may correspond to (a multiple of) the number of cores available
  if (n_folds > 16) warning("n_folds is large and may induce expensive computations.")


  check_structure_(size_p0_av_grid, "vector", "numeric", 1)
  check_natural_(size_p0_av_grid)
  if (size_p0_av_grid < 2) stop(paste("size_p0_av_grid must be at greater 1 ",
                                      "to allow for comparisons.",
                                      sep=""))
  if (size_p0_av_grid > 10) stop(paste("size_p0_av_grid is large and may ",
                                       "induce expensive computations. Choose ",
                                       "size_p0_av_grid in {2, 3, ..., 10}.",
                                       sep=""))

  p0_av_grid <- create_grid_(p, size_p0_av_grid)

  new_size <- length(p0_av_grid)
  if (size_p0_av_grid > new_size) {
    if (verbose) cat(paste("Cross-validation p0_av_grid reduced to ", new_size,
                           " elements as p is small.\n", sep = ""))
    size_p0_av_grid <- new_size
  }

  check_structure_(tol_cv, "vector", "numeric", 1)
  check_positive_(tol_cv, eps=.Machine$double.eps)

  check_structure_(maxit_cv, "vector", "numeric", 1)
  check_natural_(maxit_cv)

  check_structure_(batch_cv, "vector", "logical", 1)

  check_structure_(n_cpus, "vector", "numeric", 1)
  check_natural_(n_cpus)
  if (n_cpus > n_folds){
    message <- paste("The number of cpus in use will be at most equal to n_folds.",
                     "n_cpus is therefore set to n_folds = ", n_folds, ". \n", sep ="")
    if(verbose) cat(message)
    else warning(message)
    n_cpus <- n_folds
  }

  if (n_cpus > 1) {

    n_cpus_avail <- parallel::detectCores()
    if (n_cpus > n_cpus_avail) {
      n_cpus <- n_cpus_avail
      warning(paste("The number of CPUs specified exceeds the number of CPUs ",
                    "available on the machine. The latter has been used instead.", sep=""))
    }
    if (verbose) cat(paste("Cross-validation with ", n_cpus, " CPUs.\n",
                           "Please make sure that enough RAM is available. \n", sep=""))
  }

  n_cv <- n
  p_cv <- p

  list_cv <- create_named_list_(n_cv, p_cv, n_folds, p0_av_grid, size_p0_av_grid,
                                tol_cv, n_cpus, maxit_cv, batch_cv)
  class(list_cv) <- "cv"

  list_cv

}


create_grid_ <- function(p, size_p0_av_grid) {

  if (p < 75) { # a different treatment to avoid having a single element in the grid
    p0_av_grid <- unique(round(seq(max(floor(p/4), 1), max(ceiling(p/3), 2),
                                   length.out = size_p0_av_grid), 0))
  } else {

    p0_av_grid <- seq(max(min(1000, p/4), 1),
                      max(min(1500, p/2), 1),
                      length.out = size_p0_av_grid)

    base_round_ <- function(x, base){
      sapply( round(x / base) * base, function(el) max(el, 1) )
    }

    p0_av_grid <- unique(base_round_(p0_av_grid, 25))
  }

  p0_av_grid

}


cross_validate_ <- function(Y, X, Z, list_cv, user_seed, verbose) {


  d <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)

  if (!is.null(Z)) q <- ncol(Z)
  else q <- NULL

  with(list_cv, {

    folds <- rep_len(1:n_folds, n)

    evaluate_fold_ <- function(k) {
      if (verbose) { cat(paste("Evaluating fold k = ", k, "... \n", sep=""))
        cat("-------------------------\n")
      }

      current <- which(folds == k)

      Y_tr <- Y[-current,, drop = FALSE]
      Y_test <- Y[current,, drop = FALSE] # drop = FALSE for the case where n_folds = n

      X_tr <- X[-current,, drop = FALSE]
      X_test <- X[current,, drop = FALSE]

      # rescale X_tr as the algorithm then uses sample variance 1
      X_tr <- scale(X_tr)
      bool_cst_x <- is.nan(colSums(X_tr))
      if (any(bool_cst_x)) {
        X_tr <- X_tr[, !bool_cst_x, drop = FALSE]
        # remove the corresponding columns in X_test too so that the lower bound
        # obtained with X_tr can be evaluated on X_test
        X_test <- X_test[, !bool_cst_x, drop = FALSE]
        p <- ncol(X_tr)
      }

      Y_tr <- scale(Y_tr, center = TRUE, scale = FALSE)

      #if (family == "gaussian") Y_test <- scale(Y_test, center = TRUE, scale = FALSE)
      Y_test <- scale(Y_test, center = TRUE, scale = FALSE)

      if (!is.null(Z)) {
        Z_tr <- Z[-current,, drop = FALSE]
        Z_test <- Z[current,, drop = FALSE]
        Z_tr <- scale(Z_tr)
        bool_cst_z <- is.nan(colSums(Z_tr))
        if (any(bool_cst_z)) {
          Z_tr <- Z_tr[, !bool_cst_z, drop = FALSE]
          Z_test <- Z_test[, !bool_cst_z, drop = FALSE]
          q <- ncol(Z_tr)
        }
      } else {
        Z_tr <- Z_test <- NULL
      }

      lb_vec <- vector(length = size_p0_av_grid)

      for(ind_pg in 1:size_p0_av_grid) {

        pg <-  p0_av_grid[ind_pg]

        if (verbose) cat(paste("Evaluating p0_av = ", pg, "... \n", sep=""))

        list_hyper_pg <- auto_set_hyper_(Y_tr, p, pg, q, family = "gaussian")
        list_init_pg <- auto_set_init_(Y_tr, p, pg, q, user_seed, family = "gaussian")

        if (is.null(q)) {
          vb_tr <- locus_core_(Y_tr, X_tr, list_hyper_pg,
                               list_init_pg$gam_vb, list_init_pg$mu_beta_vb,
                               list_init_pg$sig2_beta_vb, list_init_pg$tau_vb,
                               tol_cv, maxit_cv, batch_cv, verbose = FALSE,
                               full_output = TRUE)

          lb_vec[ind_pg] <- with(vb_tr, {
            lower_bound_(Y_test, X_test, a, a_vb, b, b_vb, eta, gam_vb, kappa,
                         lambda, nu, sig2_beta_vb, sig2_inv_vb, tau_vb,
                         m1_beta, m2_beta, sum_gam)

          })
        } else {
          vb_tr <- locus_z_core_(Y_tr, X_tr, Z_tr, list_hyper_pg,
                                 list_init_pg$gam_vb, list_init_pg$mu_alpha_vb,
                                 list_init_pg$mu_beta_vb, list_init_pg$sig2_alpha_vb,
                                 list_init_pg$sig2_beta_vb, list_init_pg$tau_vb,
                                 tol_cv, maxit_cv, batch_cv, verbose = FALSE,
                                 full_output = TRUE)

          lb_vec[ind_pg] <- with(vb_tr, {
            lower_bound_z_(Y_test, X_test, Z_test, a, a_vb, b, b_vb, eta,
                           gam_vb, kappa, lambda, mu_alpha_vb, nu, phi, phi_vb,
                           sig2_alpha_vb, sig2_beta_vb, sig2_inv_vb, tau_vb,
                           xi, zeta2_inv_vb, m2_alpha, m1_beta, m2_beta, sum_gam)
          })
        }


        if (verbose) { cat(paste("Lower bound on test set, fold ", k, ", p0_av ",
                                 pg, ": ", lb_vec[ind_pg], ". \n", sep = ""))
          cat("-------------------------\n") }
      }
      lb_vec
    }

    RNGkind("L'Ecuyer-CMRG") # ensure reproducibility when using mclapply

    lb_mat <- parallel::mclapply(1:n_folds, function(k) evaluate_fold_(k),
                                 mc.cores = n_cpus)
    lb_mat <- do.call(rbind, lb_mat)

    rownames(lb_mat) <- paste("fold_", 1:n_folds, sep = "")
    colnames(lb_mat) <- paste("p0_av_", p0_av_grid, sep = "")

    p0_av_opt <- p0_av_grid[which.max(colMeans(lb_mat))]

    if (verbose) {
      cat("Lower bounds on test sets for each fold and each grid element: \n")
      print(lb_mat)
      cat(paste("===== ...end of cross-validation with selected p0_av = ",
                p0_av_opt, " ===== \n", sep=""))
    }

    p0_av_opt
  })
}

