# This file is part of the `locus` R package:
#     https://github.com/hruffieux/locus
#
# Internal core function to call the variational algorithm for structured
# sparse regression with identity link, no fixed covariates.
# See help of `locus` function for details.
#
locus_struct_core_ <- function(Y, X, list_hyper, gam_vb, mu_beta_vb, sig2_beta_vb,
                               tau_vb, vec_fac_st, tol, maxit, verbose, batch = "y",
                               full_output = FALSE, debug = FALSE) {

  # Y must have been centered, and X standardized.

  d <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)

  with(list_hyper, { # list_init not used with the with() function to avoid
                     # copy-on-write for large objects

    mu_theta_vb <- m0

    list_S0_inv <- lapply(unique(vec_fac_st), function(bl) {
      corX <- cor(X[, vec_fac_st == bl, drop = FALSE])
      corX <- as.matrix(Matrix::nearPD(corX, corr = TRUE, do2eigen = TRUE)$mat) # regularization in case of non-positive definiteness.
      as.matrix(solve(corX) / s02)})

    list_sig2_theta_vb <- update_sig2_theta_vb_(d, list_S0_inv)

    m1_beta <- update_m1_beta_(gam_vb, mu_beta_vb)
    m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE)

    mat_x_m1 <- update_mat_x_m1_(X, m1_beta)


    converged <- FALSE
    lb_new <- -Inf
    it <- 0

    while ((!converged) & (it < maxit)) {

      lb_old <- lb_new
      it <- it + 1

      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("Iteration ", format(it), "... \n", sep = ""))

      # % #
      lambda_vb <- update_lambda_vb_(lambda, sum(gam_vb))
      nu_vb <- update_nu_vb_(nu, m2_beta, tau_vb)

      sig2_inv_vb <- lambda_vb / nu_vb
      # % #

      # % #
      eta_vb <- update_eta_vb_(n, eta, gam_vb)
      kappa_vb <- update_kappa_vb_(Y, kappa, mat_x_m1, m1_beta, m2_beta, sig2_inv_vb)

      tau_vb <- eta_vb / kappa_vb
      # % #

      sig2_beta_vb <- update_sig2_beta_vb_(n, sig2_inv_vb, tau_vb)

      log_tau_vb <- update_log_tau_vb_(eta_vb, kappa_vb)
      log_sig2_inv_vb <- update_log_sig2_inv_vb_(lambda_vb, nu_vb)


      # different possible batch-coordinate ascent schemes:

      if (batch == "y") { # optimal scheme

        for (j in 1:p) {

          mat_x_m1 <- mat_x_m1 - tcrossprod(X[, j], m1_beta[j, ])

          mu_beta_vb[j, ] <- sig2_beta_vb * (tau_vb * crossprod(Y - mat_x_m1, X[, j]))

          gam_vb[j, ] <- exp(-log_one_plus_exp_(pnorm(mu_theta_vb[j], lower.tail = FALSE, log.p = TRUE) -
                                                  pnorm(mu_theta_vb[j], log.p = TRUE) -
                                                  log_tau_vb / 2 - log_sig2_inv_vb / 2 -
                                                  mu_beta_vb[j, ] ^ 2 / (2 * sig2_beta_vb) -
                                                  log(sig2_beta_vb) / 2))

          m1_beta[j, ] <- gam_vb[j, ] * mu_beta_vb[j, ]

          mat_x_m1 <- mat_x_m1 + tcrossprod(X[, j], m1_beta[j, ])

        }


      } else if (batch == "0"){ # no batch, used only internally
        # schemes "x" of "x-y" are not batch concave
        # hence not implemented as they may diverge

        for (k in 1:d) {

          for (j in 1:p) {

            mat_x_m1[, k] <- mat_x_m1[, k] - X[, j] * m1_beta[j, k]

            mu_beta_vb[j, k] <- sig2_beta_vb[k] * tau_vb[k] * crossprod(Y[, k] - mat_x_m1[, k], X[, j])

            gam_vb[j, k] <- exp(-log_one_plus_exp_(pnorm(mu_theta_vb[j], lower.tail = FALSE, log.p = TRUE) -
                                                     pnorm(mu_theta_vb[j], log.p = TRUE) -
                                                     log_tau_vb[k] / 2 - log_sig2_inv_vb / 2 -
                                                     mu_beta_vb[j, k] ^ 2 / (2 * sig2_beta_vb[k]) -
                                                     log(sig2_beta_vb[k]) / 2))

            m1_beta[j, k] <- gam_vb[j, k] * mu_beta_vb[j, k]

            mat_x_m1[, k] <- mat_x_m1[, k] + X[, j] * m1_beta[j, k]

          }
        }

      } else {

        stop ("Batch scheme not defined. Exit.")

      }

      m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE)


      W <- update_W_struct_(gam_vb, mu_theta_vb)

      mu_theta_vb <- update_mu_theta_vb_(W, m0, list_S0_inv, list_sig2_theta_vb, vec_fac_st)


      lb_new <- elbo_struct_(Y, eta, eta_vb, gam_vb, kappa, kappa_vb, lambda,
                             lambda_vb, m0, mu_theta_vb, nu, nu_vb, sig2_beta_vb,
                             list_S0_inv, list_sig2_theta_vb, sig2_inv_vb, tau_vb,
                             m1_beta, m2_beta, mat_x_m1, vec_fac_st)

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

      create_named_list_(eta, eta_vb, gam_vb, kappa, kappa_vb, lambda,
                         lambda_vb, m0, mu_theta_vb, nu, nu_vb, sig2_beta_vb,
                         list_S0_inv, list_sig2_theta_vb, sig2_inv_vb, tau_vb,
                         m1_beta, m2_beta, vec_fac_st)

    } else {

      names_x <- colnames(X)
      names_y <- colnames(Y)

      rownames(gam_vb) <- names_x
      colnames(gam_vb) <- names_y
      names(mu_theta_vb) <- names_x

      diff_lb <- abs(lb_opt - lb_old)

      create_named_list_(gam_vb, mu_theta_vb, converged, it, lb_opt, diff_lb)

    }
  })

}



# Internal function which implements the marginal log-likelihood variational
# lower bound (ELBO) corresponding to the `locus_struct_core` algorithm.
#
elbo_struct_ <- function(Y, eta, eta_vb, gam_vb, kappa, kappa_vb, lambda,
                         lambda_vb, m0, mu_theta_vb, nu, nu_vb, sig2_beta_vb,
                         list_S0_inv, list_sig2_theta_vb, sig2_inv_vb, tau_vb,
                         m1_beta, m2_beta, mat_x_m1, vec_fac_st) {


  n <- nrow(Y)

  # needed for monotonically increasing elbo.
  #
  eta_vb <- update_eta_vb_(n, eta, gam_vb)
  kappa_vb <- update_kappa_vb_(Y, kappa, mat_x_m1, m1_beta, m2_beta, sig2_inv_vb)

  lambda_vb <- update_lambda_vb_(lambda, sum(gam_vb))
  nu_vb <- update_nu_vb_(nu, m2_beta, tau_vb)

  log_tau_vb <- update_log_tau_vb_(eta_vb, kappa_vb)
  log_sig2_inv_vb <- update_log_sig2_inv_vb_(lambda_vb, nu_vb)


  elbo_A <- e_y_(n, kappa, kappa_vb, log_tau_vb, m2_beta, sig2_inv_vb, tau_vb)

  elbo_B <- e_beta_gamma_struct_(gam_vb, log_sig2_inv_vb, log_tau_vb,
                                 mu_theta_vb, m2_beta, sig2_beta_vb,
                                 list_sig2_theta_vb, sig2_inv_vb, tau_vb)

  elbo_C <- e_theta_(m0, mu_theta_vb, list_S0_inv, list_sig2_theta_vb, vec_fac_st)

  elbo_D <- e_tau_(eta, eta_vb, kappa, kappa_vb, log_tau_vb, tau_vb)

  elbo_E <- e_sig2_inv_(lambda, lambda_vb, log_sig2_inv_vb, nu, nu_vb, sig2_inv_vb)

  elbo_A + elbo_B + elbo_C + elbo_D + elbo_E

}

