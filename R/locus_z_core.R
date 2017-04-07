## with covariates
locus_z_core_ <- function(Y, X, Z, list_hyper, gam_vb, mu_alpha_vb, mu_beta_vb,
                          sig2_alpha_vb, sig2_beta_vb, tau_vb, tol, maxit,
                          verbose, batch = "x-y", full_output = FALSE, debug = FALSE) {

  # Y must have been centered, and X and Z, standardized (except the intercept in Z).

  d <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  q <- ncol(Z)

  with(list_hyper, {  # list_init not used with the with() function to avoid
                      # copy-on-write for large objects

    m2_alpha <- update_m2_alpha_(mu_alpha_vb, sig2_alpha_vb)

    m1_beta <- update_m1_beta_(gam_vb, mu_beta_vb)
    m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE)

    mat_x_m1 <- update_mat_x_m1_(X, m1_beta)
    mat_z_mu <- update_mat_z_mu_(Z, mu_alpha_vb)

    phi_vb <- update_phi_z_vb_(phi, d)

    rs_gam <- rowSums(gam_vb)
    sum_gam <- sum(rs_gam)

    converged <- FALSE
    lb_old <- -Inf
    it <- 1

    while ((!converged) & (it <= maxit)) {

      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("Iteration ", format(it), "... \n", sep = ""))

      # % #
      xi_vb <- update_xi_z_vb_(xi, tau_vb, m2_alpha) ###

      zeta2_inv_vb <- phi_vb / xi_vb
      # % #

      # % #
      lambda_vb <- update_lambda_vb_(lambda, sum_gam)
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

      digam_sum <- digamma(a + b + d)


      # different possible batch-coordinate ascent schemes:

      if (batch == "y") { # used only internally

        log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
        log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

        for (i in 1:q) {
          mat_z_mu <- mat_z_mu - tcrossprod(Z[, i], mu_alpha_vb[i, ])

          mu_alpha_vb[i, ] <- sig2_alpha_vb[i, ] * (tau_vb *
                                                      crossprod(Y  - mat_z_mu - mat_x_m1, Z[, i]))

          mat_z_mu <- mat_z_mu + tcrossprod(Z[, i], mu_alpha_vb[i, ])
        }

        # C++ Eigen call for expensive updates
        coreZLoop(X, Y, gam_vb, log_om_vb, log_1_min_om_vb, log_sig2_inv_vb,
                  log_tau_vb, m1_beta, mat_x_m1, mat_z_mu, mu_beta_vb,
                  sig2_beta_vb, tau_vb)

        rs_gam <- rowSums(gam_vb)

      } else if (batch == "x") { # used only internally

        log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
        log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

        for (k in 1:d) {

          mu_alpha_vb[, k] <- sig2_alpha_vb[, k] * tau_vb[k] *
            (crossprod(Y[, k]  - mat_z_mu[, k] - mat_x_m1[, k], Z) +  (n - 1) * mu_alpha_vb[, k])

          mat_z_mu[, k] <- Z %*% mu_alpha_vb[, k]

          mu_beta_vb[, k] <- sig2_beta_vb[k] * tau_vb[k] *
            (crossprod(Y[, k] -  mat_z_mu[, k] - mat_x_m1[, k], X) + (n - 1) * m1_beta[, k])


          gam_vb[, k] <- exp(-log_one_plus_exp_(log_1_min_om_vb - log_om_vb -
                                                  log_tau_vb[k] / 2 - log_sig2_inv_vb / 2 -
                                                  mu_beta_vb[, k] ^ 2 / (2 * sig2_beta_vb[k]) -
                                                  log(sig2_beta_vb[k]) / 2))

          m1_beta[, k] <- mu_beta_vb[, k] * gam_vb[, k]

          mat_x_m1[, k] <- X %*% m1_beta[, k]

        }

        rs_gam <- rowSums(gam_vb)

      } else if (batch == "x-y") { # optimal scheme

        log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
        log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

        mu_alpha_vb <- sweep(sig2_alpha_vb * (crossprod(Z, Y - mat_z_mu - mat_x_m1) + (n - 1) * mu_alpha_vb), 2, tau_vb, `*`)

        mat_z_mu <- Z %*% mu_alpha_vb

        # C++ Eigen call for expensive updates
        coreZBatch(X, Y, gam_vb, log_om_vb, log_1_min_om_vb, log_sig2_inv_vb,
                      log_tau_vb, m1_beta, mat_x_m1, mat_z_mu, mu_beta_vb, sig2_beta_vb, tau_vb)

        rs_gam <- rowSums(gam_vb)

      } else if (batch == "0"){ # no batch, used only internally

        for (k in 1:d) {

          log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
          log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

          for (i in 1:q) {

            mat_z_mu[, k] <- mat_z_mu[, k] - Z[, i] * mu_alpha_vb[i, k]

            mu_alpha_vb[i, k] <- sig2_alpha_vb[i, k] * tau_vb[k] *
              crossprod(Z[, i], Y[,k]  - mat_z_mu[, k] - mat_x_m1[, k])

            mat_z_mu[, k] <- mat_z_mu[, k] + Z[, i] * mu_alpha_vb[i, k]
          }

          for (j in 1:p) {

            mat_x_m1[, k] <- mat_x_m1[, k] - X[, j] * m1_beta[j, k]

            mu_beta_vb[j, k] <- sig2_beta_vb[k] * tau_vb[k] *
              crossprod(Y[,k] - mat_x_m1[, k] - mat_z_mu[, k], X[, j])

            gam_vb[j, k] <- exp(-log_one_plus_exp_(log_1_min_om_vb[j] - log_om_vb[j] -
                                                     log_tau_vb[k] / 2 - log_sig2_inv_vb / 2 -
                                                     mu_beta_vb[j, k] ^ 2 / (2 * sig2_beta_vb[k]) -
                                                     log(sig2_beta_vb[k]) / 2))

            m1_beta[j, k] <- mu_beta_vb[j, k] * gam_vb[j, k]

            mat_x_m1[, k] <- mat_x_m1[, k] + X[, j] * m1_beta[j, k]

          }

          rs_gam <- rowSums(gam_vb)

        }

      } else {

        stop ("Batch scheme not defined. Exit.")

      }


      m2_alpha <- update_m2_alpha_(mu_alpha_vb, sig2_alpha_vb)
      m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE)

      a_vb <- update_a_vb(a, rs_gam)
      b_vb <- update_b_vb(b, d, rs_gam)
      om_vb <- a_vb / (a_vb + b_vb)

      sum_gam <- sum(rs_gam)

      lb_new <- lower_bound_z_(Y, X, Z, a, a_vb, b, b_vb, eta, gam_vb, kappa,
                               lambda, mu_alpha_vb, nu, phi, phi_vb,
                               sig2_alpha_vb, sig2_beta_vb, sig2_inv_vb, tau_vb,
                               xi, zeta2_inv_vb, m2_alpha, m1_beta, m2_beta,
                               mat_x_m1, mat_z_mu, sum_gam)


      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("ELBO = ", format(lb_new), "\n\n", sep = ""))

      if (debug && lb_new < lb_old)
        stop("ELBO not increasing monotonically. Exit. ")

      converged <- (abs(lb_new-lb_old) < tol)

      lb_old <- lb_new
      it <- it + 1

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
      create_named_list_(a, a_vb, b, b_vb, eta, gam_vb, kappa, lambda,
                         mu_alpha_vb, nu, phi, phi_vb, sig2_alpha_vb,
                         sig2_beta_vb, sig2_inv_vb, tau_vb, xi, zeta2_inv_vb,
                         m2_alpha, m1_beta, m2_beta, sum_gam)
    } else {
      names_x <- colnames(X)
      names_y <- colnames(Y)
      names_z <- colnames(Z)

      rownames(gam_vb) <- names_x
      colnames(gam_vb) <- names_y
      names(om_vb) <- names_x
      rownames(mu_alpha_vb) <- names_z
      colnames(mu_alpha_vb) <- names_y

      create_named_list_(gam_vb, om_vb, mu_alpha_vb, converged, it, lb_opt)
    }
  })
}


lower_bound_z_ <- function(Y, X, Z, a, a_vb, b, b_vb, eta, gam_vb, kappa, lambda,
                           mu_alpha_vb, nu, phi, phi_vb, sig2_alpha_vb,
                           sig2_beta_vb, sig2_inv_vb, tau_vb, xi, zeta2_inv_vb,
                           m2_alpha, m1_beta, m2_beta, mat_x_m1, mat_z_mu, sum_gam) {


  n <- nrow(Y)
  q <- ncol(Z)

  xi_vb <- update_xi_z_vb_(xi, tau_vb, m2_alpha)

  eta_vb <- update_eta_z_vb_(n, q, eta, gam_vb)

  kappa_vb <- update_kappa_z_vb_(Y, X, Z, kappa, mu_alpha_vb, m1_beta, m2_alpha,
                                 m2_beta, mat_x_m1, mat_z_mu, sig2_inv_vb,
                                 zeta2_inv_vb)

  lambda_vb <- update_lambda_vb_(lambda, sum_gam)
  nu_vb <- update_nu_vb_(nu, m2_beta, tau_vb)

  log_tau_vb <- digamma(eta_vb) - log(kappa_vb)
  log_zeta2_inv_vb <- digamma(phi_vb) - log(xi_vb)
  log_sig2_inv_vb <- digamma(lambda_vb) - log(nu_vb)
  log_om_vb <- digamma(a_vb) - digamma(a_vb + b_vb)
  log_1_min_om_vb <- digamma(b_vb) - digamma(a_vb + b_vb)

  A <- sum(-n / 2 * log(2 * pi) + n / 2 * log_tau_vb -
             tau_vb * (kappa_vb - colSums(m2_beta) * sig2_inv_vb / 2 -
                         crossprod(m2_alpha, zeta2_inv_vb) / 2 - kappa))

  eps <- .Machine$double.eps # to control the argument of the log when gamma is very small
  B <- sum(log_sig2_inv_vb * gam_vb / 2 +
             sweep(gam_vb, MARGIN = 2, log_tau_vb, `*`) / 2 -
             sweep(m2_beta, MARGIN = 2, tau_vb, `*`) * sig2_inv_vb / 2 +
             sweep(gam_vb, MARGIN = 1, log_om_vb, `*`) +
             sweep(1 - gam_vb, MARGIN = 1, log_1_min_om_vb, `*`) +
             1 / 2 * sweep(gam_vb, 2, log(sig2_beta_vb) + 1, `*`) -
             gam_vb * log(gam_vb + eps) - (1 - gam_vb) * log(1 - gam_vb + eps))

  G <- sum((eta - eta_vb) * log_tau_vb - (kappa - kappa_vb) * tau_vb +
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

  A + B + G + H + J + K + L

}

