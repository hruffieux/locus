locus_z_info_core_ <- function(Y, X, Z, V, list_hyper, gam_vb, mu_alpha_vb,
                               mu_beta_vb, mu_c0_vb, mu_c_vb, sig2_alpha_vb,
                               sig2_beta_vb, tau_vb, tol, maxit, batch,
                               verbose, full_output = FALSE) {

  # Y must have been centered, and X, Z and V standardized (except intercept in Z).

  d <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  q <- ncol(Z)
  r <- ncol(V)

  with(list_hyper, { # list_init not used with the with() function to avoid
    # copy-on-write for large objects

    m2_alpha <- update_m2_alpha_(mu_alpha_vb, sig2_alpha_vb)

    m1_beta <- update_m1_beta_(gam_vb, mu_beta_vb)
    m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE)

    mat_x_m1 <- update_mat_x_m1_(X, m1_beta)
    mat_z_mu <- update_mat_z_mu_(Z, mu_alpha_vb)
    mat_v_mu <- update_mat_v_mu_(V, mu_c0_vb, mu_c_vb)

    sig2_c0_vb <- update_sig2_c0_vb_(d, s02)
    sig2_c_vb <-  update_sig2_c_vb_(p, s2)

    phi_vb <- update_phi_z_vb_(phi, d)

    converged <- FALSE
    lb_old <- -Inf
    it <- 1

    while ((!converged) & (it <= maxit)) {

      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("Iteration ", format(it), "... \n", sep = ""))

      # % #
      xi_vb <- update_xi_z_vb_(xi, tau_vb, m2_alpha)

      zeta2_inv_vb <- phi_vb / xi_vb
      # % #

      # % #
      lambda_vb <- update_lambda_vb_(lambda, sum(gam_vb))
      nu_vb <- update_nu_vb_(nu, m2_beta, tau_vb)

      sig2_inv_vb <- lambda_vb / nu_vb
      # % #

      # % #
      eta_vb <- update_eta_z_vb_(n, q, eta, gam_vb)
      kappa_vb <- update_kappa_z_vb_(Y, X, Z, kappa, mu_alpha_vb, m1_beta,
                                     m2_alpha, m2_beta, mat_x_m1, mat_z_mu,
                                     sig2_inv_vb, zeta2_inv_vb)

      tau_vb <- eta_vb / kappa_vb
      # % #

      sig2_alpha_vb <- update_sig2_alpha_vb_(n, zeta2_inv_vb, tau_vb)
      sig2_beta_vb <- update_sig2_beta_vb_(n, sig2_inv_vb, tau_vb)

      log_tau_vb <- update_log_tau_vb_(eta_vb, kappa_vb)
      log_sig2_inv_vb <- update_log_sig2_inv_vb_(lambda_vb, nu_vb)

      W <- update_W_info_(gam_vb, mat_v_mu)

      if (batch) { # some updates are made batch-wise

        for (i in 1:q) {
          mat_z_mu <- mat_z_mu - tcrossprod(Z[, i], mu_alpha_vb[i, ])

          mu_alpha_vb[i, ] <- sig2_alpha_vb[i, ] * (tau_vb *
                                                      crossprod(Y  - mat_z_mu - mat_x_m1, Z[, i]))

          mat_z_mu <- mat_z_mu + tcrossprod(Z[, i], mu_alpha_vb[i, ])
        }

        log_Phi_mat_v_mu <- pnorm(mat_v_mu, log.p = TRUE)
        log_1_min_Phi_mat_v_mu <- pnorm(mat_v_mu, lower.tail = FALSE, log.p = TRUE)

        coreZInfoLoop(X, Y, gam_vb, log_Phi_mat_v_mu, log_1_min_Phi_mat_v_mu,
                      log_sig2_inv_vb, log_tau_vb, m1_beta, mat_x_m1, mat_z_mu,
                      mu_beta_vb, sig2_beta_vb, tau_vb)


        mat_v_mu <- sweep(mat_v_mu, 1, mu_c0_vb, `-`)

        mu_c0_vb <- update_mu_c0_vb_(W, mat_v_mu, m0, s02, sig2_c0_vb)

        mat_v_mu <- sweep(mat_v_mu, 1, mu_c0_vb, `+`)


        for (l in 1:r) {

          mat_v_mu <- mat_v_mu - tcrossprod(V[, l], mu_c_vb[l, ])

          mu_c_vb[l, ] <- sig2_c_vb * crossprod(W - mat_v_mu, V[, l])

          mat_v_mu <- mat_v_mu + tcrossprod(V[, l], mu_c_vb[l, ])

        }

      } else {

        for (k in 1:d) {

          for (i in 1:q) {

            mat_z_mu[, k] <- mat_z_mu[, k] - Z[, i] * mu_alpha_vb[i, k]

            mu_alpha_vb[i, k] <- sig2_alpha_vb[i, k] * tau_vb[k] *
              crossprod(Z[, i], Y[,k]  - mat_z_mu[, k] - mat_x_m1[, k])

            mat_z_mu[, k] <- mat_z_mu[, k] + Z[, i] * mu_alpha_vb[i, k]
          }


          for (j in 1:p) {

            mat_x_m1[, k] <- mat_x_m1[, k] - X[, j] * m1_beta[j, k]

            mu_beta_vb[j, k] <- sig2_beta_vb[k] * tau_vb[k] * crossprod(Y[, k] - mat_x_m1[, k] - mat_z_mu[, k], X[, j])

            gam_vb[j, k] <- exp(-log_one_plus_exp_(pnorm(mat_v_mu[j, k], lower.tail = FALSE, log.p = TRUE) -
                                                     pnorm(mat_v_mu[j, k], log.p = TRUE) -
                                                     log_tau_vb[k] / 2 - log_sig2_inv_vb / 2 -
                                                     mu_beta_vb[j, k] ^ 2 / (2 * sig2_beta_vb[k]) -
                                                     log(sig2_beta_vb[k]) / 2))

            m1_beta[j, k] <- gam_vb[j, k] * mu_beta_vb[j, k]

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

      }

      m2_alpha <- update_m2_alpha_(mu_alpha_vb, sig2_alpha_vb)
      m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE)


      lb_new <- lower_bound_z_info_(Y, X, Z, V, eta, gam_vb, kappa, lambda,
                                    m0, mu_alpha_vb, mu_c0_vb, mu_c_vb, nu,
                                    phi, phi_vb, sig2_alpha_vb, sig2_beta_vb,
                                    sig2_c0_vb, sig2_c_vb, sig2_inv_vb, s02, s2,
                                    tau_vb,  xi, zeta2_inv_vb, m2_alpha, m1_beta,
                                    m2_beta, mat_x_m1, mat_v_mu, mat_z_mu)


      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("Lower bound = ", format(lb_new), "\n\n", sep = ""))

      converged <- (abs(lb_new - lb_old) < tol)

      lb_old <- lb_new
      it <- it + 1
    }



    if (verbose) {
      if (converged) {
        cat(paste("Convergence obtained after ", format(it),
                  " iterations with variational lower bound = ",
                  format(lb_new), ". \n\n", sep = ""))
      } else {
        warning("Maximal number of iterations reached before convergence. Exit.")
      }
    }

    lb_opt <- lb_new

    if (full_output) { # for internal use only

      create_named_list_(eta, gam_vb, kappa, lambda,
                         m0, mu_alpha_vb, mu_c0_vb, mu_c_vb, nu,
                         phi, phi_vb, sig2_alpha_vb, sig2_beta_vb,
                         sig2_c0_vb, sig2_c_vb, sig2_inv_vb, s02, s2,
                         tau_vb,  xi, zeta2_inv_vb, m2_alpha, m1_beta,
                         m2_beta, mat_x_m1, mat_v_mu, mat_z_mu)

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

      create_named_list_(gam_vb, mu_c0_vb, mu_c_vb, mu_alpha_vb, converged, it, lb_opt)

    }
  })

}




lower_bound_z_info_ <- function(Y, X, Z, V, eta, gam_vb, kappa, lambda,
                                m0, mu_alpha_vb, mu_c0_vb, mu_c_vb, nu,
                                phi, phi_vb, sig2_alpha_vb, sig2_beta_vb,
                                sig2_c0_vb, sig2_c_vb, sig2_inv_vb, s02, s2,
                                tau_vb,  xi, zeta2_inv_vb, m2_alpha, m1_beta,
                                m2_beta, mat_x_m1, mat_v_mu, mat_z_mu) {

  n <- nrow(Y)
  q <- ncol(Z)

  xi_vb <- update_xi_z_vb_(xi, tau_vb, m2_alpha)

  eta_vb <- update_eta_z_vb_(n, q, eta, gam_vb)

  kappa_vb <- update_kappa_z_vb_(Y, X, Z, kappa, mu_alpha_vb, m1_beta, m2_alpha,
                                 m2_beta, mat_x_m1, mat_z_mu, sig2_inv_vb,
                                 zeta2_inv_vb)

  lambda_vb <- update_lambda_vb_(lambda, sum(gam_vb))
  nu_vb <- update_nu_vb_(nu, m2_beta, tau_vb)

  log_tau_vb <- update_log_tau_vb_(eta_vb, kappa_vb)
  log_zeta2_inv_vb <- update_log_zeta2_inv_vb_(phi_vb, xi_vb)
  log_sig2_inv_vb <- update_log_sig2_inv_vb_(lambda_vb, nu_vb)

  A <- sum(-n / 2 * log(2 * pi) + n / 2 * log_tau_vb -
             tau_vb * (kappa_vb - colSums(m2_beta) * sig2_inv_vb / 2 -
                         crossprod(m2_alpha, zeta2_inv_vb) / 2  - kappa))

  eps <- .Machine$double.eps^0.75 # to control the argument of the log when gamma is very small

  B <- sum(log_sig2_inv_vb * gam_vb / 2 +
             sweep(gam_vb, 2, log_tau_vb, `*`) / 2 -
             sweep(m2_beta, 2, tau_vb, `*`) * sig2_inv_vb / 2 +
             gam_vb * pnorm(mat_v_mu, log.p = TRUE) +
             sweep((1 - gam_vb) * pnorm(mat_v_mu, lower.tail = FALSE, log.p = TRUE), 1, sig2_c_vb * rowSums(V^2) / 2, `-`) -
             sig2_c0_vb / 2 + 1 / 2 * sweep(gam_vb, 2, log(sig2_beta_vb) + 1, `*`) -
             gam_vb * log(gam_vb + eps) - (1 - gam_vb) * log(1 - gam_vb + eps))

  G <- sum(log(sig2_c0_vb) + 1 - log(s02) - (mu_c0_vb^2 + sig2_c0_vb - 2*mu_c0_vb * m0 + m0^2) / s02) / 2

  H <- sum(log(sig2_c_vb) + 1 - log(s2) - (mu_c_vb^2 + sig2_c_vb) / s2) / 2

  J <- sum((eta - eta_vb) * log_tau_vb -
             (kappa - kappa_vb) * tau_vb + eta * log(kappa) -
             eta_vb * log(kappa_vb) - lgamma(eta) + lgamma(eta_vb))

  K <- (lambda - lambda_vb) * log_sig2_inv_vb - (nu - nu_vb) * sig2_inv_vb +
    lambda * log(nu) - lambda_vb * log(nu_vb) - lgamma(lambda) +
    lgamma(lambda_vb)

  L <- sum(sweep( sweep( sweep( sweep(m2_alpha, MARGIN = 2, tau_vb, `*`),
                                MARGIN = 1, - zeta2_inv_vb / 2, `*`),
                         MARGIN = 2, log_tau_vb / 2, `+`),
                  MARGIN = 1, log_zeta2_inv_vb / 2, `+`) +
             log(sig2_alpha_vb) / 2 + 1 / 2)

  M <- sum((phi - phi_vb) * log_zeta2_inv_vb - (xi - xi_vb) * zeta2_inv_vb +
             phi * log(xi) - phi_vb * log(xi_vb) - lgamma(phi) + lgamma(phi_vb))

  A + B + G + H + J + K + L + M

}

