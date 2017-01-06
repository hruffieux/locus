#' @export
set_cv <- function(n, p, n_folds, size_p_guess_grid, n_cpus, tol_cv, maxit_cv = 1000,
                   batch_cv = T, verbose = T) {

  check_structure_(n_folds, "vector", "numeric", 1)
  check_natural_(n_folds)
  if (n_folds > n) stop("n_folds must not exceed the number of observations.")

  # 16 may correspond to (a multiple of) the number of cores available
  if (n_folds > 16) warning("n_folds is large and may induce expensive computations.")


  check_structure_(size_p_guess_grid, "vector", "numeric", 1)
  check_natural_(size_p_guess_grid)
  if (size_p_guess_grid < 2) stop(paste("size_p_guess_grid must be at greater 1 ",
                                        "to allow for comparisons.",
                                        sep=""))
  if (size_p_guess_grid > 10) stop(paste("size_p_guess_grid is large and may ",
                                         "induce expensive computations. Choose ",
                                         "size_p_guess_grid in {2, 3, ..., 10}.",
                                         sep=""))

  p_guess_grid <- create_grid_(p, size_p_guess_grid)

  new_size <- length(p_guess_grid)
  if (size_p_guess_grid > new_size) {
    if (verbose) cat(paste("Cross-validation p_guess_grid reduced to ", new_size,
                           " elements as p is small.\n", sep = ""))
    size_p_guess_grid <- new_size
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
                           "Please make sure that enough RAM is available.", sep=""))
  }

  n_cv <- n
  p_cv <- p

  list_cv <- create_named_list_(n_cv, p_cv, n_folds, p_guess_grid, size_p_guess_grid,
                                tol_cv, n_cpus, maxit_cv, batch_cv)
  class(list_cv) <- "cv"

  list_cv

}


create_grid_ <- function(p, size_p_guess_grid) {

  if (p < 75) { # a different treatment to avoid having a single element in the grid
    p_guess_grid <- unique(round(seq(max(floor(p/4), 1), max(ceiling(p/3), 2),
                                     length.out = size_p_guess_grid), 0))
  } else {

    p_guess_grid <- seq(max(min(1000, p/4), 1),
                        max(min(1500, p/2), 1),
                        length.out = size_p_guess_grid)

    base_round_ <- function(x, base){
      sapply( round(x / base) * base, function(el) max(el, 1) )
    }

    p_guess_grid <- unique(base_round_(p_guess_grid, 25))
  }

  p_guess_grid

}


cross_validate_ <- function(Y, X, Z, d, n, p, q, list_cv, user_seed, verbose) {

  list2env(list_cv, envir = environment())
  rm(list_cv)

  folds <- rep_len(1:n_folds, n)

  evaluate_fold_ <- function(k) {
    if (verbose) { cat(paste("Evaluating fold k = ", k, "... \n", sep=""))
      cat("-------------------------\n")
    }

    current <- which(folds == k)

    n_test <- length(current)
    n_tr <- n - n_test

    Y_tr <- Y[-current,, drop = F]
    Y_test <- Y[current,, drop = F] # drop = F for the case where n_folds = n

    X_tr <- X[-current,, drop = F]
    X_test <- X[current,, drop = F]

    # rescale X_tr as the algorithm then uses sample variance 1
    X_tr <- scale(X_tr)
    bool_cst_x <- is.nan(colSums(X_tr))
    if (any(bool_cst_x)) {
      X_tr <- X_tr[, !bool_cst_x, drop = F]
      # remove the corresponding columns in X_test too so that the lower bound
      # obtained with X_tr can be evaluated on X_test
      X_test <- X_test[, !bool_cst_x, drop = F]
      p <- ncol(X_tr)
    }

    Y_tr <- scale(Y_tr, center = T, scale = F)
    Y_test <- scale(Y_test, center = T, scale = F)


    if (!is.null(Z)) {
      Z_tr <- Z[-current,, drop = F]
      Z_test <- Z[current,, drop = F]
      Z_tr <- scale(Z_tr)
      bool_cst_z <- is.nan(colSums(Z_tr))
      if (any(bool_cst_z)) {
        Z_tr <- Z_tr[, !bool_cst_z, drop = F]
        Z_test <- Z_test[, !bool_cst_z, drop = F]
        q <- ncol(Z_tr)
      }
    } else {
      Z_tr <- Z_test <- NULL
    }

    lb_vec <- vector(length = size_p_guess_grid)

    for(ind_pg in 1:size_p_guess_grid) {

      pg <-  p_guess_grid[ind_pg]

      if (verbose) cat(paste("Evaluating p_guess = ", pg, "... \n", sep=""))

      list_hyper_pg <- auto_set_hyperparam_(Y_tr, p, pg, q)
      list_init_pg <- auto_init_param_(Y_tr, p, pg, user_seed, q)

      if (is.null(q)) {
        vb_tr <- locus_core_(Y_tr, X_tr, d, n_tr, p, list_hyper_pg, list_init_pg,
                             tol_cv, maxit_cv, batch_cv, verbose = F, full_output = T)

        list2env(vb_tr, envir=environment())
        rm(vb_tr)

        lb_vec[ind_pg] <- lower_bound_(Y_test, X_test, d, n_test, p,
                                       mu_beta_vb, sig2_beta_vb, sig2_inv_vb,
                                       tau_vb, gam_vb, om_vb, eta, kappa,
                                       lambda, nu, a, b, a_vb, b_vb,
                                       m1_beta, m2_beta, sum_gam)
      } else {
        vb_tr <- locus_z_core_(Y_tr, X_tr, Z_tr, d, n_tr, p, q, list_hyper_pg,
                               list_init_pg, tol_cv, maxit_cv, batch_cv, verbose = F,
                               full_output = T)
        list2env(vb_tr, envir=environment())
        rm(vb_tr)

        lb_vec[ind_pg] <- lower_bound_z_(Y_test, X_test, Z_test, d, n_test, p, q,
                                         mu_alpha_vb, sig2_alpha_vb, zeta2_inv_vb,
                                         mu_beta_vb, sig2_beta_vb, sig2_inv_vb,
                                         tau_vb, gam_vb, om_vb,
                                         eta, kappa, lambda, nu, a, b, a_vb, b_vb,
                                         phi, phi_vb, xi,
                                         m2_alpha, m1_beta, m2_beta, sum_gam)
      }

      if (verbose) { cat(paste("Lower bound on test set, fold ", k, ", p_guess ",
                               pg, ": ", lb_vec[ind_pg], ". \n", sep = ""))
        cat("-------------------------\n") }
    }
    lb_vec
  }

  RNGkind("L'Ecuyer-CMRG") # ensure reproducibility when using mclapply

  lb_mat <- parallel::mclapply(1:n_folds, function(k) evaluate_fold_(k), mc.cores = n_cpus)
  lb_mat <- do.call(rbind, lb_mat)

  rownames(lb_mat) <- paste("fold_", 1:n_folds, sep = "")
  colnames(lb_mat) <- paste("p_guess_", p_guess_grid, sep = "")

  p_guess_opt <- p_guess_grid[which.max(colMeans(lb_mat))]

  if (verbose) {
    cat("Lower bounds on test sets for each fold and each grid element: \n")
    print(lb_mat)
    cat(paste("===== ...end of cross-validation with selected p_guess = ",
              p_guess_opt, " ===== \n", sep=""))
  }

  p_guess_opt

}

