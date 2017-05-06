# This file is part of the `locus` R package:
#     https://github.com/hruffieux/locus
#
# Internal core function to call the variational algorithm for logit link,
# optional fixed covariates and external annotation variables.
# See help of `locus` function for details.
#
locus_logit_info_core_ <- function(Y, X, Z, V, list_hyper, chi_vb, gam_vb,
                                   mu_alpha_vb, mu_beta_vb, mu_c0_vb, mu_c_vb,
                                   sig2_alpha_vb, sig2_beta_vb, tol, maxit,
                                   verbose, batch = "y", full_output = FALSE,
                                   debug = FALSE) {

  # 1/2 must have been substracted from Y, and X, Z and V must have been standardized (except intercept in Z).

  d <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  q <- ncol(Z)
  r <- ncol(V)

  with(list_hyper, { # list_init not used with the with() function to avoid
    # copy-on-write for large objects

    m2_alpha <- update_m2_alpha_(mu_alpha_vb, sig2_alpha_vb)
    m1_beta <- update_m1_beta_(gam_vb, mu_beta_vb)
    m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb)

    mat_x_m1 <- update_mat_x_m1_(X, m1_beta)
    mat_z_mu <- update_mat_z_mu_(Z, mu_alpha_vb)
    mat_v_mu <- update_mat_v_mu_(V, mu_c0_vb, mu_c_vb)

    sig2_c0_vb <- update_sig2_c0_vb_(d, s02)
    sig2_c_vb <-  update_sig2_c_vb_(p, s2)

    phi_vb <- update_phi_z_vb_(phi, d)

    converged <- FALSE
    lb_new <- -Inf
    it <- 0

    while ((!converged) & (it < maxit)) {

      lb_old <- lb_new
      it <- it + 1

      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("Iteration ", format(it), "... \n", sep = ""))

      # % #
      chi_vb <- update_chi_vb_(X, Z, m1_beta, m2_beta, mat_x_m1, mat_z_mu, sig2_alpha_vb)

      psi_vb <- update_psi_logit_vb_(chi_vb)
      # % #

      # % #
      xi_vb <- update_xi_bin_vb_(xi, m2_alpha)

      zeta2_inv_vb <- phi_vb / xi_vb
      # % #

      # % #
      lambda_vb <- update_lambda_vb_(lambda, sum(gam_vb))
      nu_vb <- update_nu_bin_vb_(nu, m2_beta)

      sig2_inv_vb <- lambda_vb / nu_vb
      # % #

      sig2_alpha_vb <- update_sig2_alpha_logit_vb_(Z, psi_vb, zeta2_inv_vb)
      sig2_beta_vb <- update_sig2_beta_logit_vb_(X, psi_vb, sig2_inv_vb)

      log_sig2_inv_vb <- update_log_sig2_inv_vb_(lambda_vb, nu_vb)

      W <- update_W_info_(gam_vb, mat_v_mu)


      # different possible batch-coordinate ascent schemes:

      if (batch == "y") { # some updates are made batch-wise

        for (i in 1:q) {
          mat_z_mu <- mat_z_mu - tcrossprod(Z[, i], mu_alpha_vb[i, ])

          mu_alpha_vb[i, ] <- sig2_alpha_vb[i, ] * crossprod(Y - 2 * psi_vb * (mat_z_mu + mat_x_m1), Z[, i])

          mat_z_mu <- mat_z_mu + tcrossprod(Z[, i], mu_alpha_vb[i, ])
        }

        log_Phi_mat_v_mu <- pnorm(mat_v_mu, log.p = TRUE)
        log_1_min_Phi_mat_v_mu <- pnorm(mat_v_mu, lower.tail = FALSE, log.p = TRUE)

        # C++ Eigen call for expensive updates
        coreLogitInfoLoop(X, Y, gam_vb, log_Phi_mat_v_mu, log_1_min_Phi_mat_v_mu,
                          log_sig2_inv_vb, m1_beta, mat_x_m1, mat_z_mu, mu_beta_vb,
                          psi_vb, sig2_beta_vb)


        mat_v_mu <- sweep(mat_v_mu, 1, mu_c0_vb, `-`)

        mu_c0_vb <- update_mu_c0_vb_(W, mat_v_mu, m0, s02, sig2_c0_vb)

        mat_v_mu <- sweep(mat_v_mu, 1, mu_c0_vb, `+`)


        for (l in 1:r) {

          mat_v_mu <- mat_v_mu - tcrossprod(V[, l], mu_c_vb[l, ])

          mu_c_vb[l, ] <- sig2_c_vb * crossprod(W - mat_v_mu, V[, l])

          mat_v_mu <- mat_v_mu + tcrossprod(V[, l], mu_c_vb[l, ])

        }

      } else if (batch == "0"){

        for (k in 1:d) {

          for (i in 1:q) {

            mat_z_mu[, k] <- mat_z_mu[, k] - Z[, i] * mu_alpha_vb[i, k]

            mu_alpha_vb[i, k] <- sig2_alpha_vb[i, k] * crossprod(Z[, i], Y[, k] - 2 * psi_vb[, k] * (mat_z_mu[, k] + mat_x_m1[, k]))

            mat_z_mu[, k] <- mat_z_mu[, k] + Z[, i] * mu_alpha_vb[i, k]

          }

          for (j in 1:p) {

            mat_x_m1[, k] <- mat_x_m1[, k] - X[, j] * m1_beta[j, k]

            mu_beta_vb[j, k] <- sig2_beta_vb[j, k] *
              crossprod(X[, j], Y[, k] - 2 * psi_vb[, k] * (mat_z_mu[, k] + mat_x_m1[, k]))

            gam_vb[j, k] <- exp(-log_one_plus_exp_(pnorm(mat_v_mu[j, k], lower.tail = FALSE, log.p = TRUE) -
                                                     pnorm(mat_v_mu[j, k], log.p = TRUE) -
                                                     log_sig2_inv_vb / 2 -
                                                     mu_beta_vb[j, k] ^ 2 / (2 * sig2_beta_vb[j, k]) -
                                                     log(sig2_beta_vb[j, k]) / 2))

            m1_beta[j, k] <- mu_beta_vb[j, k] * gam_vb[j, k]

            mat_x_m1[, k] <- mat_x_m1[, k] + X[, j] * m1_beta[j, k]

          }

          mat_v_mu <- sweep(mat_v_mu, 1, mu_c0_vb, `-`)

          mu_c0_vb <- update_mu_c0_vb_(W, mat_v_mu, m0, s02, sig2_c0_vb)

          mat_v_mu <- sweep(mat_v_mu, 1, mu_c0_vb, `+`)

          for (l in 1:r) {

            mat_v_mu[, k] <- mat_v_mu[, k] - V[, l] * mu_c_vb[l, k]

            mu_c_vb[l, k] <- sig2_c_vb * crossprod(W[, k] - mat_v_mu[, k], V[, l])

            mat_v_mu[, k] <- mat_v_mu[, k] + V[, l] * mu_c_vb[l, k]

          }

        }

      } else {

        stop ("Batch scheme not defined. Exit.")

      }

      m2_alpha <- update_m2_alpha_(mu_alpha_vb, sig2_alpha_vb)
      m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb)


      lb_new <- elbo_logit_info_(Y, X, Z, V, chi_vb, gam_vb, m0,  mu_c0_vb,
                                 mu_c_vb, lambda, nu, phi, phi_vb, psi_vb,
                                 sig2_alpha_vb, sig2_beta_vb, sig2_c0_vb,
                                 sig2_c_vb, sig2_inv_vb, s02, s2, xi,
                                 zeta2_inv_vb, mu_alpha_vb, m2_alpha, m1_beta,
                                 m2_beta, mat_x_m1, mat_v_mu, mat_z_mu)

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
      create_named_list_(chi_vb, gam_vb, m0,  mu_c0_vb, mu_c_vb, lambda, nu,
                         phi, phi_vb, psi_vb, sig2_alpha_vb, sig2_beta_vb,
                         sig2_c0_vb, sig2_c_vb, sig2_inv_vb, s02, s2, xi,
                         zeta2_inv_vb, mu_alpha_vb, m2_alpha, m1_beta, m2_beta,
                         mat_x_m1, mat_v_mu, mat_z_mu)
    } else {

      names_x <- colnames(X)
      names_y <- colnames(Y)
      names_v <- colnames(V)
      names_z <- colnames(Z)

      rownames(mu_alpha_vb) <- names_z
      colnames(mu_alpha_vb) <- names_y

      rownames(gam_vb) <- names_x
      colnames(gam_vb) <- names_y

      names(mu_c0_vb) <- names_x

      rownames(mu_c_vb) <- names_v
      colnames(mu_c_vb) <- names_y

      diff_lb <- abs(lb_opt - lb_old)

      create_named_list_(gam_vb, mu_c0_vb, mu_c_vb, mu_alpha_vb, converged, it,
                         lb_opt, diff_lb)

    }
  })

}


# Internal function which implements the marginal log-likelihood variational
# lower bound (ELBO) corresponding to the `locus_logit_info_core` algorithm.
#
elbo_logit_info_ <- function(Y, X, Z, V, chi_vb, gam_vb, m0,  mu_c0_vb,
                             mu_c_vb, lambda, nu, phi, phi_vb, psi_vb,
                             sig2_alpha_vb, sig2_beta_vb, sig2_c0_vb,
                             sig2_c_vb, sig2_inv_vb, s02, s2, xi,
                             zeta2_inv_vb, mu_alpha_vb, m2_alpha,
                             m1_beta, m2_beta, mat_x_m1, mat_v_mu, mat_z_mu) {

  lambda_vb <- update_lambda_vb_(lambda, sum(gam_vb))
  nu_vb <- update_nu_bin_vb_(nu, m2_beta)

  xi_vb <- update_xi_bin_vb_(xi, m2_alpha)

  log_sig2_inv_vb <- update_log_sig2_inv_vb_(lambda_vb, nu_vb)
  log_zeta2_inv_vb <- update_log_zeta2_inv_vb_(phi_vb, xi_vb)


  elbo_A <- e_y_logit_(X, Y, Z, chi_vb, m1_beta, m2_alpha, m2_beta, mat_x_m1,
                       mat_z_mu, mu_alpha_vb, psi_vb)

  elbo_B <- e_beta_gamma_info_bin_(V, gam_vb, log_sig2_inv_vb, mat_v_mu, m2_beta,
                                   sig2_beta_vb, sig2_c0_vb, sig2_c_vb, sig2_inv_vb)

  elbo_C <- e_c0_(m0, mu_c0_vb, s02, sig2_c0_vb)

  elbo_D <- e_c_(mu_c_vb, s2, sig2_c_vb)

  elbo_E <- e_sig2_inv_(lambda, lambda_vb, log_sig2_inv_vb, nu, nu_vb, sig2_inv_vb)

  elbo_F <- e_alpha_logit_(m2_alpha, log_zeta2_inv_vb, sig2_alpha_vb, zeta2_inv_vb)

  elbo_G <- e_zeta2_inv_(log_zeta2_inv_vb, phi, phi_vb, xi, xi_vb, zeta2_inv_vb)

  elbo_A + elbo_B + elbo_C + elbo_D + elbo_E + elbo_F + elbo_G

}

