# This file is part of the `locus` R package:
#     https://github.com/hruffieux/locus
#
# Internal core function to call the variational algorithm for identity-probit link,
# optional fixed covariates and no external annotation variables.
# See help of `locus` function for details.
#
locus_mix_core_ <- function(Y, X, Z, ind_bin, list_hyper, gam_vb, mu_alpha_vb,
                            mu_beta_vb, sig2_alpha_vb, sig2_beta_vb, tau_vb,
                            tol, maxit, verbose, batch = "y",
                            full_output = FALSE, debug = FALSE) {

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
    digam_sum <- digamma(a + b + d)

    log_tau_vb <- rep(0, d)

    converged <- FALSE
    lb_new <- -Inf
    it <- 0

    while ((!converged) & (it < maxit)) {

      lb_old <- lb_new
      it <- it + 1

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

      kappa_vb <- update_kappa_z_vb_(Y_cont, Z, kappa,
                                     mu_alpha_vb[, -ind_bin, drop = FALSE],
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


      # different possible batch-coordinate ascent schemes:

      if (batch == "y") { # optimal scheme

        log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
        log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

        for (i in sample(1:q)) {
          mat_z_mu <- mat_z_mu - tcrossprod(Z[, i], mu_alpha_vb[i, ])

          mu_alpha_vb[i, ] <- sig2_alpha_vb[i, ] * (tau_vb *
                                                      crossprod(W  - mat_z_mu - mat_x_m1, Z[, i]))

          mat_z_mu <- mat_z_mu + tcrossprod(Z[, i], mu_alpha_vb[i, ])
        }

        # C++ Eigen call for expensive updates
        shuffled_ind <- as.numeric(sample(0:(p-1))) # Zero-based index in C++

        coreZLoop(X, W, gam_vb, log_om_vb, log_1_min_om_vb, log_sig2_inv_vb,
                  log_tau_vb, m1_beta, mat_x_m1, mat_z_mu, mu_beta_vb,
                  sig2_beta_vb, tau_vb, shuffled_ind)

        rs_gam <- rowSums(gam_vb)

      } else if (batch == "x") { # used internally for testing purposes,
        # convergence not ensured as ELBO not batch-concave

        log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
        log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

        for (k in sample(1:d)) {

          mu_alpha_vb[, k] <- sig2_alpha_vb[, k] * tau_vb[k] *
            (crossprod(W[, k]  - mat_z_mu[, k] - mat_x_m1[, k], Z) +  (n - 1) * mu_alpha_vb[, k])

          mu_alpha_vb[1, k] <- mu_alpha_vb[1, k] + sig2_alpha_vb[1, k] * tau_vb[k] * mu_alpha_vb[1, k]

          mat_z_mu[, k] <- Z %*% mu_alpha_vb[, k]

          mu_beta_vb[, k] <- sig2_beta_vb[k] * tau_vb[k] * (crossprod(W[, k] -  mat_z_mu[, k] - mat_x_m1[, k], X) +
                                                              (n - 1) * m1_beta[, k])


          gam_vb[, k] <- exp(-log_one_plus_exp_(log_1_min_om_vb - log_om_vb -
                                                  log_tau_vb[k] / 2 - log_sig2_inv_vb / 2 -
                                                  mu_beta_vb[, k] ^ 2 / (2 * sig2_beta_vb[k]) -
                                                  log(sig2_beta_vb[k]) / 2))

          m1_beta[, k] <- mu_beta_vb[, k] * gam_vb[, k]

          mat_x_m1[, k] <- X %*% m1_beta[, k]

        }

        rs_gam <- rowSums(gam_vb)

      } else if (batch == "x-y") { # used only internally, convergence not ensured

        log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
        log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

        mu_alpha_vb <- sweep(sig2_alpha_vb * (crossprod(Z, W - mat_z_mu - mat_x_m1) + (n - 1) * mu_alpha_vb), 2, tau_vb, `*`)

        mu_alpha_vb[1, ] <- mu_alpha_vb[1, ] + sig2_alpha_vb[1, ] * tau_vb * mu_alpha_vb[1, ]

        mat_z_mu <- Z %*% mu_alpha_vb

        # C++ Eigen call for expensive updates
        coreZBatch(X, W, gam_vb, log_om_vb, log_1_min_om_vb, log_sig2_inv_vb,
                   log_tau_vb, m1_beta, mat_x_m1, mat_z_mu, mu_beta_vb,
                   sig2_beta_vb, tau_vb)

        rs_gam <- rowSums(gam_vb)

      } else if (batch == "0"){ # no batch, used only internally

        for (k in sample(1:d)) {

          log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
          log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

          for (i in sample(1:q)) {

            mat_z_mu[, k] <- mat_z_mu[, k] - Z[, i] * mu_alpha_vb[i, k]

            mu_alpha_vb[i, k] <- sig2_alpha_vb[i, k] * tau_vb[k] *
              crossprod(W[, k] - mat_z_mu[, k] - mat_x_m1[, k], Z[, i])

            mat_z_mu[, k] <- mat_z_mu[, k] + Z[, i] * mu_alpha_vb[i, k]
          }

          for (j in sample(1:p)) {

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

      } else {

        stop ("Batch scheme not defined. Exit.")

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

      lb_new <- elbo_mix_(Y_bin, Y_cont, ind_bin, X, Z, a, a_vb, b, b_vb, eta,
                          gam_vb, kappa, lambda, mu_alpha_vb, nu, phi, phi_vb,
                          sig2_alpha_vb, sig2_beta_vb, sig2_inv_vb, tau_vb,
                          log_tau_vb, xi, zeta2_inv_vb, m2_alpha, m1_beta,
                          m2_beta, mat_x_m1, mat_z_mu, sum_gam)

      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("ELBO = ", format(lb_new), "\n\n", sep = ""))

      if (debug && lb_new < lb_old)
        stop("ELBO not increasing monotonically. Exit. ")

      converged <- (abs(lb_new-lb_old) < tol)

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
      create_named_list_(ind_bin, a, a_vb, b, b_vb, eta, gam_vb, kappa, lambda,
                         mu_alpha_vb, nu, phi, phi_vb, sig2_alpha_vb,
                         sig2_beta_vb, sig2_inv_vb, tau_vb, log_tau_vb, xi,
                         zeta2_inv_vb, m2_alpha, m1_beta, m2_beta, sum_gam)
    } else {
      names_x <- colnames(X)
      names_y <- colnames(W)
      names_z <- colnames(Z)

      rownames(gam_vb) <- names_x
      colnames(gam_vb) <- names_y
      names(om_vb) <- names_x
      rownames(mu_alpha_vb) <- names_z
      colnames(mu_alpha_vb) <- names_y

      diff_lb <- abs(lb_opt - lb_old)

      create_named_list_(gam_vb, om_vb, mu_alpha_vb, converged, it, lb_opt, diff_lb)
    }
  })
}


# Internal function which implements the marginal log-likelihood variational
# lower bound (ELBO) corresponding to the `locus_mix_core` algorithm.
#
elbo_mix_ <- function(Y_bin, Y_cont, ind_bin, X, Z, a, a_vb, b, b_vb, eta,
                      gam_vb, kappa, lambda, mu_alpha_vb, nu, phi, phi_vb,
                      sig2_alpha_vb, sig2_beta_vb, sig2_inv_vb, tau_vb,
                      log_tau_vb, xi, zeta2_inv_vb, m2_alpha, m1_beta, m2_beta,
                      mat_x_m1, mat_z_mu, sum_gam) {

  n <- nrow(Z)
  q <- ncol(Z)

  xi_vb <- update_xi_z_vb_(xi, tau_vb, m2_alpha)

  eta_vb <- update_eta_z_vb_(n, q, eta, gam_vb[, -ind_bin, drop = FALSE])


  kappa_vb <- update_kappa_z_vb_(Y_cont, Z, kappa,
                                 mu_alpha_vb[, -ind_bin, drop = FALSE],
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


  elbo_A_cont <- e_y_(n, kappa, kappa_vb, log_tau_vb[-ind_bin],
                      m2_beta[, -ind_bin, drop = FALSE], sig2_inv_vb,
                      tau_vb[-ind_bin], m2_alpha[, -ind_bin, drop = FALSE],
                      zeta2_inv_vb)

  elbo_A_bin <- e_y_probit_(X, Y_bin, Z, m1_beta[, ind_bin, drop = FALSE],
                            m2_beta[, ind_bin, drop = FALSE],
                            mat_x_m1[, ind_bin, drop = FALSE],
                            mat_z_mu[, ind_bin, drop = FALSE],
                            sig2_alpha_vb[, ind_bin, drop = FALSE], sweep = FALSE)


  elbo_B <- e_beta_gamma_(gam_vb, log_om_vb, log_1_min_om_vb, log_sig2_inv_vb,
                          log_tau_vb, m2_beta, sig2_beta_vb, sig2_inv_vb, tau_vb)

  elbo_C <- e_tau_(eta, eta_vb, kappa, kappa_vb, log_tau_vb[-ind_bin], tau_vb[-ind_bin])

  elbo_D <- e_sig2_inv_(lambda, lambda_vb, log_sig2_inv_vb, nu, nu_vb, sig2_inv_vb)

  elbo_E <- e_omega_(a, a_vb, b, b_vb, log_om_vb, log_1_min_om_vb)

  elbo_F <- e_alpha_(m2_alpha, log_tau_vb, log_zeta2_inv_vb, sig2_alpha_vb,
                     tau_vb, zeta2_inv_vb)

  elbo_G <- e_zeta2_inv_(log_zeta2_inv_vb, phi, phi_vb, xi, xi_vb, zeta2_inv_vb)

  elbo_A_cont + elbo_A_bin + elbo_B + elbo_C + elbo_D + elbo_E + elbo_F + elbo_G

}

