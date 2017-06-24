# This file is part of the `locus` R package:
#     https://github.com/hruffieux/locus
#
# Internal core function to call the variational algorithm for group selection
# with identity link, no fixed covariates and no external annotation variables.
# See help of `locus` function for details.
#
locus_group_core_ <- function(Y, list_X, list_hyper, gam_vb, list_mu_beta_vb,
                               sig2_inv_vb, tau_vb, tol, maxit, verbose,
                               batch = "y", full_output = FALSE, debug = FALSE) {


  # Y must have been centered, and X, standardized.

  d <- ncol(Y)
  n <- nrow(Y)
  G <- length(list_X)

  g_sizes <- sapply(list_X, ncol)

  with(list_hyper, { # list_init not used with the with() function to avoid
                     # copy-on-write for large objects

    list_sig2_beta_star_inv <- lapply(list_X, function(X_g) crossprod(X_g) + diag(sig2_inv_vb, nrow = ncol(X_g)))

    list_sig2_beta_star <- lapply(list_sig2_beta_star_inv, solve)

    list_m1_beta <- update_g_m1_beta_(list_mu_beta_vb, gam_vb)

    list_m1_btb <- update_g_m1_btb_(gam_vb, list_mu_beta_vb, list_sig2_beta_star, tau_vb)

    list_m1_btXtXb <- update_g_m1_btXtXb_(list_X, gam_vb, list_mu_beta_vb,
                                          list_sig2_beta_star, tau_vb)

    mat_x_m1 <- update_g_mat_x_m1_(list_X, list_m1_beta)


    log_tau_vb <- update_log_tau_vb_(eta, kappa)  # do not update tau_vb here as
                                                  # its current form was already used
                                                  # in list_m1_btb as part of the vb
                                                  # parameter sig2_beta = sig2_beta_star / tau_vb
    rs_gam <- rowSums(gam_vb)
    digam_sum <- digamma(a + b + d)

    converged <- FALSE
    lb_new <- -Inf
    it <- 0

    while ((!converged) & (it < maxit)) {

      lb_old <- lb_new
      it <- it + 1

      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("Iteration ", format(it), "... \n", sep = ""))

      # % #
      lambda_vb <- update_g_lambda_vb_(lambda, g_sizes, rs_gam)
      nu_vb <- update_g_nu_vb_(nu, list_m1_btb, tau_vb)

      list_sig2_beta_star_inv <- lapply(list_sig2_beta_star_inv, function(sig2_beta_star_inv)
        sig2_beta_star_inv - diag(sig2_inv_vb, nrow = nrow(sig2_beta_star_inv))) # to avoid recomputing X_g^TX_g each time

      sig2_inv_vb <- lambda_vb / nu_vb

      list_sig2_beta_star_inv <- lapply(list_sig2_beta_star_inv, function(sig2_beta_star_inv)
        sig2_beta_star_inv + diag(sig2_inv_vb, nrow = nrow(sig2_beta_star_inv)))

      list_sig2_beta_star <- lapply(list_sig2_beta_star_inv, solve)

      vec_log_det <- log_det(list_sig2_beta_star)
      # % #

      log_sig2_inv_vb <- update_log_sig2_inv_vb_(lambda_vb, nu_vb)


      # different possible batch-coordinate ascent schemes:

      if (batch == "y") { # optimal scheme

        log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
        log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

        for (g in sample(1:G)) {
          mat_x_m1 <- mat_x_m1 - list_X[[g]] %*% list_m1_beta[[g]]

          list_mu_beta_vb[[g]] <- list_sig2_beta_star[[g]] %*% crossprod(list_X[[g]], Y - mat_x_m1)

          gam_vb[g, ] <- exp(-log_one_plus_exp_(log_1_min_om_vb[g] - log_om_vb[g] -
                                                  g_sizes[g] * (log_sig2_inv_vb + log_tau_vb - log(tau_vb)) / 2 - # |g| * log(tau_vb) /2 came out of the determinant
                                                  colSums(list_mu_beta_vb[[g]] *
                                                            (list_sig2_beta_star_inv[[g]] %*% list_mu_beta_vb[[g]])) * tau_vb / 2 -
                                                  vec_log_det[g] / 2))

          list_m1_beta[[g]] <- sweep(list_mu_beta_vb[[g]], 2, gam_vb[g, ], `*`)

          mat_x_m1 <- mat_x_m1 + list_X[[g]] %*% list_m1_beta[[g]]
        }

        rs_gam <- rowSums(gam_vb)

      } else if (batch == "x") { # used only internally, convergence not ensured

        stop("Not implemented")

      } else if (batch == "x-y") { # used only internally, convergence not ensured

        stop("Not implemented")

      } else if (batch == "0") { # no batch, used only internally

        for (k in sample(1:d)) {

          log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
          log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

          for (g in sample(1:G)) {

            mat_x_m1[, k] <- mat_x_m1[, k] - list_X[[g]] %*% list_m1_beta[[g]][, k]

            list_mu_beta_vb[[g]][, k] <- list_sig2_beta_star[[g]] %*% crossprod(list_X[[g]], Y[, k] - mat_x_m1[, k])

            gam_vb[g, k] <- exp(-log_one_plus_exp_(log_1_min_om_vb[g] - log_om_vb[g] -
                                                     g_sizes[g] * (log_sig2_inv_vb + log_tau_vb[k] - log(tau_vb[k])) / 2 -
                                                     sum(list_mu_beta_vb[[g]][, k] * (list_sig2_beta_star_inv[[g]] %*% list_mu_beta_vb[[g]][, k])) * tau_vb[k] / 2 -
                                                     vec_log_det[g] / 2))

            list_m1_beta[[g]][, k] <- list_mu_beta_vb[[g]][, k] * gam_vb[g, k]

            mat_x_m1[, k] <- mat_x_m1[, k] + list_X[[g]] %*% list_m1_beta[[g]][, k]

          }

          rs_gam <- rowSums(gam_vb)

        }

      } else {

        stop ("Batch scheme not defined. Exit.")

      }


      list_m1_btb <- update_g_m1_btb_(gam_vb, list_mu_beta_vb, list_sig2_beta_star, tau_vb)

      list_m1_btXtXb <- update_g_m1_btXtXb_(list_X, gam_vb, list_mu_beta_vb,
                                            list_sig2_beta_star, tau_vb)


      a_vb <- update_a_vb(a, rs_gam)
      b_vb <- update_b_vb(b, d, rs_gam)
      om_vb <- a_vb / (a_vb + b_vb)


      # % #
      eta_vb <- update_g_eta_vb_(n, eta, g_sizes, gam_vb)
      kappa_vb <- update_g_kappa_vb_(Y, list_X, kappa, list_m1_beta, list_m1_btb,
                                     list_m1_btXtXb, mat_x_m1, sig2_inv_vb)


      lb_new <- elbo_group_(Y, list_X, a, a_vb, b, b_vb, eta, eta_vb, g_sizes,
                            gam_vb, kappa, kappa_vb, lambda, lambda_vb, nu, nu_vb,
                            rs_gam, list_sig2_beta_star, sig2_inv_vb, tau_vb,
                            vec_log_det, list_m1_beta, list_m1_btb,
                            list_m1_btXtXb, mat_x_m1)

      tau_vb <- eta_vb / kappa_vb # has to be updated after the elbo, as list_sig2_beta_star depends on it.
      # % #

      log_tau_vb <- update_log_tau_vb_(eta_vb, kappa_vb)

      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("ELBO = ", format(lb_new), "\n\n", sep = ""))

      if (debug && lb_new < lb_old)
        stop("ELBO not increasing monotonically. Exit. ")

      converged <- (abs(lb_new - lb_old) < tol)

    }


    if (verbose) {
      if (converged) {
        cat(paste("Convergence obtained after ", format(it), " iterations. \n",
                  "Optimal marginal log-likelihood variational lower bound ",
                  "(ELBO) = ", format(lb_new), ". \n\n", sep = ""))
      } else {
        warning("Maximal number of iterations reached before convergence. Exit.")
      }
    }

    lb_opt <- lb_new

    if (full_output) { # for internal use only
      create_named_list_(a, a_vb, b, b_vb, eta, eta_vb, g_sizes,
                         gam_vb, kappa, kappa_vb, lambda, lambda_vb, nu, nu_vb,
                         rs_gam, list_sig2_beta_star, sig2_inv_vb, tau_vb,
                         vec_log_det, list_m1_beta, list_m1_btb, list_m1_btXtXb)
    } else {
      names_y <- colnames(Y)

      names_G <- unlist(lapply(list_X,
                               function(X_g) paste(as.character(colnames(X_g)), collapse = "-")))

      rownames(gam_vb) <- names_G
      colnames(gam_vb) <- names_y
      names(om_vb) <- names_G

      diff_lb <- abs(lb_opt - lb_old)

      create_named_list_(gam_vb, om_vb, converged, it, lb_opt, diff_lb)
    }
  })

}


# Internal function which implements the marginal log-likelihood variational
# lower bound (ELBO) corresponding to the `locus_group_core` algorithm.
#
elbo_group_ <- function(Y, list_X, a, a_vb, b, b_vb, eta, eta_vb, g_sizes,
                        gam_vb, kappa, kappa_vb, lambda, lambda_vb, nu, nu_vb,
                        rs_gam, list_sig2_beta_star, sig2_inv_vb, tau_vb,
                        vec_log_det, list_m1_beta, list_m1_btb, list_m1_btXtXb,
                        mat_x_m1) {

  n <- nrow(Y)

  log_tau_vb <- digamma(eta_vb) - log(kappa_vb)
  log_sig2_inv_vb <- digamma(lambda_vb) - log(nu_vb)
  log_om_vb <- digamma(a_vb) - digamma(a_vb + b_vb)
  log_1_min_om_vb <- digamma(b_vb) - digamma(a_vb + b_vb)


  elbo_A <- e_g_y_(n, kappa, kappa_vb, list_m1_btb, log_tau_vb, sig2_inv_vb, tau_vb)

  elbo_B <- e_g_beta_gamma_(gam_vb, g_sizes, log_om_vb, log_1_min_om_vb,
                            log_sig2_inv_vb, log_tau_vb, list_m1_btb,
                            list_sig2_beta_star, sig2_inv_vb, tau_vb, vec_log_det)

  elbo_C <- e_tau_(eta, eta_vb, kappa, kappa_vb, log_tau_vb, tau_vb)

  elbo_D <- e_sig2_inv_(lambda, lambda_vb, log_sig2_inv_vb, nu, nu_vb, sig2_inv_vb)

  elbo_E <- e_omega_(a, a_vb, b, b_vb, log_om_vb, log_1_min_om_vb)


  elbo_A + elbo_B + elbo_C + elbo_D + elbo_E

}

