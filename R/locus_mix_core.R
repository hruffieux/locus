## for continuous and binary responses using a probit link for the latter.
locus_mix_core_ <- function(Y, X, Z, ind_bin, list_hyper, gam_vb, mu_alpha_vb,
                            mu_beta_vb, sig2_alpha_vb, sig2_beta_vb, tau_vb,
                            tol, maxit, batch, verbose, full_output = FALSE) {

  # Y must have its continuous variables centered,
  # and X and Z must have been standardized (except intercept in Z).

  d <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  q <- ncol(Z)

  W <- Y
  Y_bin <- Y[, ind_bin, drop = FALSE]
  Y_cont <- Y[, -ind_bin, drop = FALSE]
  rm(Y)

  with(list_hyper, {  # list_init not used with the with() function to avoid
                      # copy-on-write for large objects
    m2_alpha <- update_m2_alpha_(mu_alpha_vb, sig2_alpha_vb)

    m1_beta <- update_m1_beta_(gam_vb, mu_beta_vb)
    m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE)

    mat_x_m1 <- update_mat_x_m1_(X, m1_beta)
    mat_z_mu <- update_mat_z_mu_(Z, mu_alpha_vb)

    # no drop = FALSE for W, as replacement not allowed in this case
    W[, ind_bin] <- update_W_probit_(Y_bin,
                                     mat_z_mu[, ind_bin, drop = FALSE],
                                     mat_x_m1[, ind_bin, drop = FALSE])

    phi_vb <- update_phi_z_vb_(phi, d)

    rs_gam <- rowSums(gam_vb)
    sum_gam <- sum(rs_gam)

    log_tau_vb <- rep(0, d)

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
      lambda_vb <- update_lambda_vb_(lambda, sum_gam)
      nu_vb <- update_nu_vb_(nu, m2_beta, tau_vb)

      sig2_inv_vb <- lambda_vb / nu_vb
      # % #

      # % #
      eta_vb <- update_eta_z_vb_(n, q, eta, gam_vb[, -ind_bin, drop = FALSE])

      kappa_vb <- update_kappa_z_vb_(Y_cont, X, Z,
                                     kappa, mu_alpha_vb[, -ind_bin, drop = FALSE],
                                     m1_beta[, -ind_bin, drop = FALSE],
                                     m2_alpha[, -ind_bin, drop = FALSE],
                                     m2_beta[, -ind_bin, drop = FALSE],
                                     mat_x_m1[, -ind_bin, drop = FALSE],
                                     mat_z_mu[, -ind_bin, drop = FALSE],
                                     sig2_inv_vb, zeta2_inv_vb, intercept = TRUE)

      tau_vb[-ind_bin] <- eta_vb / kappa_vb
      # % #

      sig2_alpha_vb <- update_sig2_alpha_vb_(n, zeta2_inv_vb, tau_vb, intercept = TRUE)
      sig2_beta_vb <- update_sig2_beta_vb_(n, sig2_inv_vb, tau_vb)

      log_tau_vb[-ind_bin] <- update_log_tau_vb_(eta_vb, kappa_vb)
      log_sig2_inv_vb <- update_log_sig2_inv_vb_(lambda_vb, nu_vb)

      digam_sum <- digamma(a + b + d)

      if (batch == "y") {

        log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
        log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

        for (i in 1:q) {
          mat_z_mu <- mat_z_mu - tcrossprod(Z[, i], mu_alpha_vb[i, ])

          mu_alpha_vb[i, ] <- sig2_alpha_vb[i, ] * (tau_vb *
                                                      crossprod(W  - mat_z_mu - mat_x_m1, Z[, i]))

          mat_z_mu <- mat_z_mu + tcrossprod(Z[, i], mu_alpha_vb[i, ])
        }

        coreZLoop(X, W, gam_vb, log_om_vb, log_1_min_om_vb, log_sig2_inv_vb,
                  log_tau_vb, m1_beta, mat_x_m1, mat_z_mu, mu_beta_vb,
                  sig2_beta_vb, tau_vb)

        rs_gam <- rowSums(gam_vb)

      } else {

        for (k in 1:d) {

          log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
          log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

          for (i in 1:q) {

            mat_z_mu[, k] <- mat_z_mu[, k] - Z[, i] * mu_alpha_vb[i, k]

            mu_alpha_vb[i, k] <- sig2_alpha_vb[i, k] * tau_vb[k] *
              crossprod(W[, k] - mat_z_mu[, k] - mat_x_m1[, k], Z[, i])

            mat_z_mu[, k] <- mat_z_mu[, k] + Z[, i] * mu_alpha_vb[i, k]
          }

          for (j in 1:p) {

            mat_x_m1[, k] <- mat_x_m1[, k] - X[, j] * m1_beta[j, k]

            mu_beta_vb[j, k] <- sig2_beta_vb[k] * tau_vb[k] *
              crossprod(W[, k] - mat_x_m1[, k] - mat_z_mu[, k], X[, j])

            gam_vb[j, k] <- exp(-log_one_plus_exp_(log_1_min_om_vb[j] - log_om_vb[j] -
                                                    log_tau_vb[k] / 2 - log_sig2_inv_vb / 2 -
                                                    mu_beta_vb[j, k] ^ 2 / (2 * sig2_beta_vb[k]) -
                                                    log(sig2_beta_vb[k]) / 2))

            m1_beta[j, k] <- mu_beta_vb[j, k] * gam_vb[j, k]

            mat_x_m1[, k] <- mat_x_m1[, k] + X[, j] * m1_beta[j, k]

          }

          rs_gam <- rowSums(gam_vb)

        }

      }

      m2_alpha <- update_m2_alpha_(mu_alpha_vb, sig2_alpha_vb)
      m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE)

      W[, ind_bin] <- update_W_probit_(Y_bin,
                                       mat_z_mu[, ind_bin, drop = FALSE],
                                       mat_x_m1[, ind_bin, drop = FALSE])

      a_vb <- update_a_vb(a, rs_gam)
      b_vb <- update_b_vb(b, d, rs_gam)
      om_vb <- a_vb / (a_vb + b_vb)

      sum_gam <- sum(rs_gam)

      lb_new <- lower_bound_mix_(Y_bin, Y_cont, ind_bin, X, Z, a, a_vb, b, b_vb,
                                 eta, gam_vb, kappa, lambda, mu_alpha_vb, nu,
                                 phi, phi_vb, sig2_alpha_vb, sig2_beta_vb,
                                 sig2_inv_vb, tau_vb, log_tau_vb, xi,
                                 zeta2_inv_vb, m2_alpha, m1_beta, m2_beta,
                                 mat_x_m1, mat_z_mu, sum_gam)

      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("Lower bound = ", format(lb_new), "\n\n", sep = ""))

      converged <- (abs(lb_new-lb_old) < tol)

      lb_old <- lb_new
      it <- it + 1

    }

    if (verbose) {
      if (converged) {
        cat(paste("Convergence obtained after ", format(it),
                  " iterations with variational lower bound = ",
                  format(lb_new), ". \n\n", sep = ""))
      } else {
        cat("Maximal number of iterations reached before convergence. Exit.")
      }
    }

    lb_opt <- lb_new

    if (full_output) { # for internal use only
      create_named_list_(ind_bin, a, a_vb, b, b_vb, eta, gam_vb, kappa, lambda,
                         mu_alpha_vb, nu, phi, phi_vb, sig2_alpha_vb,
                         sig2_beta_vb, sig2_inv_vb, tau_vb, log_tau_vb, xi,
                         zeta2_inv_vb, m2_alpha, m1_beta, m2_beta, mat_x_m1,
                         mat_z_mu, sum_gam)
    } else {
      names_x <- colnames(X)
      names_y <- colnames(W)
      names_z <- colnames(Z)

      rownames(gam_vb) <- names_x
      colnames(gam_vb) <- names_y
      names(om_vb) <- names_x
      rownames(mu_alpha_vb) <- names_z
      colnames(mu_alpha_vb) <- names_y

      create_named_list_(lb_opt, gam_vb, om_vb, mu_alpha_vb)
    }
  })
}


lower_bound_mix_ <- function(Y_bin, Y_cont, ind_bin, X, Z, a, a_vb, b, b_vb, eta,
                             gam_vb, kappa, lambda, mu_alpha_vb, nu, phi,
                             phi_vb, sig2_alpha_vb, sig2_beta_vb, sig2_inv_vb,
                             tau_vb, log_tau_vb, xi, zeta2_inv_vb, m2_alpha,
                             m1_beta, m2_beta, mat_x_m1, mat_z_mu, sum_gam) {

  n <- nrow(Z)
  q <- ncol(Z)

  xi_vb <- update_xi_z_vb_(xi, tau_vb, m2_alpha)

  eta_vb <- update_eta_z_vb_(n, q, eta, gam_vb[, -ind_bin, drop = FALSE])


  kappa_vb <- update_kappa_z_vb_(Y_cont, X, Z,
                                 kappa, mu_alpha_vb[, -ind_bin, drop = FALSE],
                                 m1_beta[, -ind_bin, drop = FALSE],
                                 m2_alpha[, -ind_bin, drop = FALSE],
                                 m2_beta[, -ind_bin, drop = FALSE],
                                 mat_x_m1[, -ind_bin, drop = FALSE],
                                 mat_z_mu[, -ind_bin, drop = FALSE],
                                 sig2_inv_vb, zeta2_inv_vb, intercept = TRUE)

  lambda_vb <- update_lambda_vb_(lambda, sum_gam)
  nu_vb <- update_nu_vb_(nu, m2_beta, tau_vb)

  log_tau_vb[-ind_bin] <- digamma(eta_vb) - log(kappa_vb)
  log_zeta2_inv_vb <- digamma(phi_vb) - log(xi_vb)
  log_sig2_inv_vb <- digamma(lambda_vb) - log(nu_vb)
  log_om_vb <- digamma(a_vb) - digamma(a_vb + b_vb)
  log_1_min_om_vb <- digamma(b_vb) - digamma(a_vb + b_vb)

  A_cont <- sum(-n / 2 * log(2 * pi) + n / 2 * log_tau_vb[-ind_bin] -
             tau_vb[-ind_bin] * (kappa_vb -
                                   colSums(m2_beta[, -ind_bin, drop = FALSE]) * sig2_inv_vb / 2 -
                                   crossprod(m2_alpha[, -ind_bin, drop = FALSE], zeta2_inv_vb) / 2 - kappa))

  U <- mat_x_m1[, ind_bin, drop = FALSE] + mat_z_mu[, ind_bin, drop = FALSE]

  A_bin <- sum(Y_bin * pnorm(U, log.p = TRUE) +
                 (1 - Y_bin) * pnorm(U, lower.tail = FALSE, log.p = TRUE) -
                 Z^2 %*% sig2_alpha_vb[, ind_bin, drop = FALSE]  / 2 -
                 X^2 %*% (m2_beta[, ind_bin, drop = FALSE] - m1_beta[, ind_bin, drop = FALSE]^2) / 2)


  eps <- .Machine$double.eps # to control the argument of the log when gamma is very small
  B <- sum(log_sig2_inv_vb * gam_vb / 2 +
             sweep(gam_vb, MARGIN = 2, log_tau_vb, `*`) / 2 -
             sweep(m2_beta, MARGIN = 2, tau_vb, `*`) * sig2_inv_vb / 2 +
             sweep(gam_vb, MARGIN = 1, log_om_vb, `*`) +
             sweep(1 - gam_vb, MARGIN = 1, log_1_min_om_vb, `*`) +
             sweep(gam_vb, 2, log(sig2_beta_vb) + 1, `*`) / 2 -
             gam_vb * log(gam_vb + eps) - (1 - gam_vb) * log(1 - gam_vb + eps))

  G <- sum((eta - eta_vb) * log_tau_vb[-ind_bin] - (kappa - kappa_vb) * tau_vb[-ind_bin] +
             eta * log(kappa) - eta_vb * log(kappa_vb) - lgamma(eta) +
             lgamma(eta_vb))

  H <- (lambda - lambda_vb) * log_sig2_inv_vb - (nu - nu_vb) * sig2_inv_vb +
    lambda * log(nu) - lambda_vb * log(nu_vb) - lgamma(lambda) + lgamma(lambda_vb)

  J <- sum((a - a_vb) * log_om_vb + (b - b_vb) * log_1_min_om_vb - lbeta(a, b) +
             lbeta(a_vb, b_vb))

  K <- sum(sweep( sweep( sweep( sweep(m2_alpha, MARGIN = 2, tau_vb, `*`),
                                MARGIN = 1, - zeta2_inv_vb / 2, `*`),
                         MARGIN = 2, log_tau_vb / 2, `+`),
                  MARGIN = 1, log_zeta2_inv_vb / 2, `+`) +
             log(sig2_alpha_vb) / 2 + 1 / 2)

  L <- sum((phi - phi_vb) * log_zeta2_inv_vb - (xi - xi_vb) * zeta2_inv_vb +
             phi * log(xi) - phi_vb * log(xi_vb) - lgamma(phi) + lgamma(phi_vb))


  A_cont + A_bin + B + G + H + J + K + L

}

