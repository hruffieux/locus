locus_logit_info_core_ <- function(Y, X, Z, V, list_hyper, chi_vb, gam_vb,
                                   mu_alpha_vb, mu_beta_vb, mu_c0_vb, mu_c_vb,
                                   sig2_alpha_vb, sig2_beta_vb, tol, maxit,
                                   batch, verbose, full_output = FALSE) {

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
    lb_old <- -Inf
    it <- 1

    while ((!converged) & (it <= maxit)) {

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

      if (batch == "y") { # some updates are made batch-wise

        for (i in 1:q) {
          mat_z_mu <- mat_z_mu - tcrossprod(Z[, i], mu_alpha_vb[i, ])

          mu_alpha_vb[i, ] <- sig2_alpha_vb[i, ] * crossprod(Y - 2 * psi_vb * (mat_z_mu + mat_x_m1), Z[, i])

          mat_z_mu <- mat_z_mu + tcrossprod(Z[, i], mu_alpha_vb[i, ])
        }

        log_Phi_mat_v_mu <- pnorm(mat_v_mu, log.p = TRUE)
        log_1_min_Phi_mat_v_mu <- pnorm(mat_v_mu, lower.tail = FALSE, log.p = TRUE)

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

      } else {

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

      }

      m2_alpha <- update_m2_alpha_(mu_alpha_vb, sig2_alpha_vb)
      m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb)


      lb_new <- lower_bound_logit_info_(Y, X, Z, V, chi_vb, gam_vb, m0,  mu_c0_vb,
                                        mu_c_vb, lambda, nu, phi, phi_vb, psi_vb,
                                        sig2_alpha_vb, sig2_beta_vb, sig2_c0_vb,
                                        sig2_c_vb, sig2_inv_vb, s02, s2, xi,
                                        zeta2_inv_vb, mu_alpha_vb, m2_alpha,
                                        m1_beta, m2_beta, mat_x_m1, mat_v_mu,
                                        mat_z_mu)

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

      create_named_list_(gam_vb, mu_c0_vb, mu_c_vb, mu_alpha_vb, converged, it, lb_opt)

    }
  })

}


lower_bound_logit_info_ <- function(Y, X, Z, V, chi_vb, gam_vb, m0,  mu_c0_vb,
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

  A <- sum(log_sigmoid(chi_vb)  + Y * (mat_x_m1 + mat_z_mu)  -  chi_vb / 2 -
             psi_vb * (X^2 %*% m2_beta + mat_x_m1^2 - X^2 %*% m1_beta^2 +
                         Z^2 %*% m2_alpha + mat_z_mu^2 - Z^2 %*% mu_alpha_vb^2 +
                         2 * mat_x_m1 * mat_z_mu - chi_vb^2))

  eps <- .Machine$double.eps^0.75 # to control the argument of the log when gamma is very small
  B <- sum(gam_vb * log_sig2_inv_vb / 2 - m2_beta * sig2_inv_vb / 2 +
             gam_vb * pnorm(mat_v_mu, log.p = TRUE) +
             sweep((1 - gam_vb) * pnorm(mat_v_mu, lower.tail = FALSE, log.p = TRUE), 1, sig2_c_vb * rowSums(V^2) / 2, `-`) -
             sig2_c0_vb / 2 + gam_vb * (log(sig2_beta_vb) + 1) / 2 -
             gam_vb * log(gam_vb + eps) - (1 - gam_vb) * log(1 - gam_vb + eps))

  G <- sum(log(sig2_c0_vb) + 1 - log(s02) - (mu_c0_vb^2 + sig2_c0_vb - 2*mu_c0_vb * m0 + m0^2) / s02) / 2

  H <- sum(log(sig2_c_vb) + 1 - log(s2) - (mu_c_vb^2 + sig2_c_vb) / s2) / 2


  J <- sum((lambda - lambda_vb) * log_sig2_inv_vb - (nu - nu_vb) * sig2_inv_vb +
             lambda * log(nu) - lambda_vb * log(nu_vb) - lgamma(lambda) +
             lgamma(lambda_vb))

  K <- sum( sweep(-sweep(m2_alpha, 1, zeta2_inv_vb, `*`), 1,
                  log_zeta2_inv_vb, `+`) + log(sig2_alpha_vb) + 1) / 2

  L <- sum((phi - phi_vb) * log_zeta2_inv_vb - (xi - xi_vb) * zeta2_inv_vb +
             phi * log(xi) - phi_vb * log(xi_vb) - lgamma(phi) +
             lgamma(phi_vb))

  A + B + G + H + J + K + L

}

