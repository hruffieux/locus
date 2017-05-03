# This file is part of the `locus` R package:
#     https://github.com/hruffieux/locus
#
# Internal core function to call the variational algorithm for identity link, no
# fixed covariates and no external annotation variables.
# See help of `locus` function for details.
#
locus_group_core_ <- function(Y, list_X, list_hyper, gam_vb, list_mu_beta_vb,
                              sig2_inv_vb, tau_vb, tol, maxit, verbose,
                              batch = "y", full_output = FALSE, debug = FALSE) {

  # Y must have been centered, and X, standardized.

  d <- ncol(Y)
  n <- nrow(Y)
  G <- length(list_X)
  g_size <- sapply(list_X, ncol)

  with(list_hyper, { # list_init not used with the with() function to avoid
                     # copy-on-write for large objects

    list_sig2_beta_star_inv <- lapply(list_X, function(X_g) crossprod(X_g) + sig2_beta_inv)
    list_sig2_beta_star <- lapply(list_sig2_beta_star_inv, solve)

    list_m1_beta <- update_g_m1_beta_(G, list_mu_beta_vb, gam_vb)
    list_m2_beta <- update_g_m2_beta_(gam_vb, list_mu_beta_vb,
                                      list_sig2_beta_star, tau_vb)
    list_m2_X_beta <- update_g_m2_X_beta_(list_X, gam_vb, list_mu_beta_vb,
                                          list_sig2_beta_star, tau_vb)

    mat_x_m1 <- update_g_mat_x_m1_(list_X, list_m1_beta)

    rs_gam <- rowSums(gam_vb)
    sum_gam <- sum(rs_gam)

    converged <- FALSE
    lb_new <- -Inf
    it <- 0

    while ((!converged) & (it <= maxit)) {

      lb_old <- lb_new
      it <- it + 1

      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("Iteration ", format(it), "... \n", sep = ""))

      # % #
      lambda_vb <- update_g_lambda_vb_(lambda, g_size, rs_gam)
      nu_vb <- update_g_nu_vb_(nu, list_m2_beta, tau_vb)

      list_sig2_beta_star_inv <- mapply(`-`, list_sig2_beta_star_inv,
                                        sig2_inv_vb, SIMPLIFY = F) # to avoid recomputing X_g^TX_g each time

      sig2_inv_vb <- lambda_vb / nu_vb

      list_sig2_beta_star_inv <- mapply(`+`, list_sig2_beta_star_inv,
                                        sig2_inv_vb, SIMPLIFY = F)
      list_sig2_beta_star <- lapply(list_sig2_beta_star_inv, solve)
      # % #

      # % #
      eta_vb <- update_g_eta_vb_(n, eta, g_size, gam_vb)
      kappa_vb <- update_g_kappa_vb_(Y, kappa, list_m1_beta, list_m2_beta,
                                     list_m2_X_beta, sig2_inv_vb)

      tau_vb <- eta_vb / kappa_vb
      # % #

      log_tau_vb <- update_log_tau_vb_(eta_vb, kappa_vb)
      log_sig2_inv_vb <- update_log_sig2_inv_vb_(lambda_vb, nu_vb)

      digam_sum <- digamma(a + b + d)


      # different possible batch-coordinate ascent schemes:

      if (batch == "y") { # optimal scheme

        log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
        log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)


        for (g in 1:G) {
          mat_x_m1 <- mat_x_m1 - list_X[[g]] %*% list_m1_beta[[g]]

          list_mu_beta_vb[[g]] <- list_sig2_beta_star[[g]] * crossprod(list_X[[g]], Y - mat_x_m1)

          gam_vb[g, ] <- exp(-log_one_plus_exp_(log_1_min_om_vb[g] - log_om_vb[g] -
                                                  g_size[g] * log_sig2_inv_vb / 2 -
                                                  sweep(sweep(crossprod(list_mu_beta_vb[[g]],
                                                                  list_sig2_beta_star_inv[[g]] %*% list_mu_beta_vb[[g]]) / 2, 2, tau_vb, `*`),
                                                        2, g_size[g] * (log(tau_vb) - log_tau_vb) / 2, `+`)  -
                                                  log(det(list_sig2_beta_star[[g]]) / 2)))

          list_m1_beta[[g]] <- sweep(list_mu_beta_vb[[g]], 2, gam_vb[g, ], `*`)

          mat_x_m1 <- mat_x_m1 + list_X[[g]] %*% list_m1_beta[[g]]
        }


        # C++ Eigen call for expensive updates
#        coreLoop(X, Y, gam_vb, log_om_vb, log_1_min_om_vb, log_sig2_inv_vb,
#                 log_tau_vb, m1_beta, mat_x_m1, mu_beta_vb, sig2_beta_vb, tau_vb)


        rs_gam <- rowSums(gam_vb)

      } else if (batch == "x") { # used only internally, convergence not ensured


        stop("Not implemented")

        # log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
        # log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)
        #
        # for (k in 1:d) {
        #
        #   mu_beta_vb[, k] <- sig2_beta_vb[k] * tau_vb[k] *
        #     (crossprod(Y[, k] - mat_x_m1[, k], X) + (n - 1) * m1_beta[, k])
        #
        #
        #   gam_vb[, k] <- exp(-log_one_plus_exp_(log_1_min_om_vb - log_om_vb -
        #                                           log_tau_vb[k] / 2 - log_sig2_inv_vb / 2 -
        #                                           mu_beta_vb[, k] ^ 2 / (2 * sig2_beta_vb[k]) -
        #                                           log(sig2_beta_vb[k]) / 2))
        #
        #   m1_beta[, k] <- mu_beta_vb[, k] * gam_vb[, k]
        #
        #   mat_x_m1[, k] <- X %*% m1_beta[, k]
        #
        # }
        #
        # rs_gam <- rowSums(gam_vb)

      } else if (batch == "x-y") { # used only internally, convergence not ensured

        stop("Not implemented")

        # log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
        # log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)
        #
        # # C++ Eigen call for expensive updates
        # coreBatch(X, Y, gam_vb, log_om_vb, log_1_min_om_vb, log_sig2_inv_vb,
        #           log_tau_vb, m1_beta, mat_x_m1, mu_beta_vb, sig2_beta_vb, tau_vb)
        #
        # rs_gam <- rowSums(gam_vb)

      } else if (batch == "0") { # no batch, used only internally


        stop("Not implemented")

        # for (k in 1:d) {
        #
        #   log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
        #   log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)
        #
        #   for (j in 1:p) {
        #
        #     mat_x_m1[, k] <- mat_x_m1[, k] - X[, j] * m1_beta[j, k]
        #
        #     mu_beta_vb[j, k] <- sig2_beta_vb[k] * tau_vb[k] * crossprod(Y[, k] - mat_x_m1[, k], X[, j])
        #
        #     gam_vb[j, k] <- exp(-log_one_plus_exp_(log_1_min_om_vb[j] - log_om_vb[j] -
        #                                              log_tau_vb[k] / 2 - log_sig2_inv_vb / 2 -
        #                                              mu_beta_vb[j, k] ^ 2 / (2 * sig2_beta_vb[k]) -
        #                                              log(sig2_beta_vb[k]) / 2))
        #
        #     m1_beta[j, k] <- mu_beta_vb[j, k] * gam_vb[j, k]
        #
        #     mat_x_m1[, k] <- mat_x_m1[, k] + X[, j] * m1_beta[j, k]
        #
        #   }
        #
        #   rs_gam <- rowSums(gam_vb)
        #
        # }

      } else {

        stop ("Batch scheme not defined. Exit.")

      }


      list_m2_beta <- update_g_m2_beta_(gam_vb, list_mu_beta_vb,
                                        list_sig2_beta_star, tau_vb)
      list_m2_X_beta <- update_g_m2_X_beta_(list_X, gam_vb, list_mu_beta_vb,
                                            list_sig2_beta_star, tau_vb)


      a_vb <- update_a_vb(a, rs_gam)
      b_vb <- update_b_vb(b, d, rs_gam)
      om_vb <- a_vb / (a_vb + b_vb)

      sum_gam <- sum(rs_gam)

      lb_new <- elbo_(Y, list_X, a, a_vb, b, b_vb, eta, gam_vb, kappa, lambda, nu,
                      list_sig2_beta_star, list_sig2_beta_star_inv, sig2_inv_vb,
                      tau_vb, list_m1_beta, list_m2_beta, mat_x_m1, sum_gam)

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
      create_named_list_(a, a_vb, b, b_vb, eta, gam_vb, kappa, lambda, nu,
                         list_sig2_beta_star, list_sig2_beta_star_inv, sig2_inv_vb,
                         tau_vb, list_m1_beta, list_m2_beta, sum_gam)
    } else {
      names_x <- colnames(X)
      names_y <- colnames(Y)

      rownames(gam_vb) <- names_x
      colnames(gam_vb) <- names_y
      names(om_vb) <- names_x

      diff_lb <- abs(lb_opt - lb_old)

      create_named_list_(gam_vb, om_vb, converged, it, lb_opt, diff_lb)
    }
  })

}


# Internal function which implements the marginal log-likelihood variational
# lower bound (ELBO) corresponding to the `locus_core` algorithm.
#
elbo_ <- function(Y, X, a, a_vb, b, b_vb, eta, gam_vb, kappa, lambda, nu,
                  sig2_beta_vb, sig2_inv_vb, tau_vb, m1_beta, m2_beta, mat_x_m1,
                  sum_gam) {

  n <- nrow(Y)

  eta_vb <- update_g_eta_vb_(n, eta, g_size, gam_vb)
  kappa_vb <- update_g_kappa_vb_(Y, kappa, list_m1_beta, list_m2_beta,
                                 list_m2_X_beta, sig2_inv_vb)

  lambda_vb <- update_g_lambda_vb_(lambda, g_size, rs_gam)
  nu_vb <- update_g_nu_vb_(nu, list_m2_beta, tau_vb)

  log_tau_vb <- digamma(eta_vb) - log(kappa_vb)
  log_sig2_inv_vb <- digamma(lambda_vb) - log(nu_vb)
  log_om_vb <- digamma(a_vb) - digamma(a_vb + b_vb)
  log_1_min_om_vb <- digamma(b_vb) - digamma(a_vb + b_vb)


  elbo_A <- e_y_(n, kappa, kappa_vb, log_tau_vb, do.call(rbind, list_m2_beta),
                 sig2_inv_vb, tau_vb)

  elbo_B <- e_g_beta_gamma_(gam_vb, g_sizes, log_om_vb, log_1_min_om_vb,
                            log_sig2_inv_vb, log_tau_vb, list_m2_beta,
                            sig2_beta_vb, sig2_inv_vb, tau_vb)

  elbo_C <- e_tau_(eta, eta_vb, kappa, kappa_vb, log_tau_vb, tau_vb)

  elbo_D <- e_sig2_inv_(lambda, lambda_vb, log_sig2_inv_vb, nu, nu_vb, sig2_inv_vb)

  elbo_E <- e_omega_(a, a_vb, b, b_vb, log_om_vb, log_1_min_om_vb)


  elbo_A + elbo_B + elbo_C + elbo_D + elbo_E

}

