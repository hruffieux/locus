# This file is part of the `locus` R package:
#     https://github.com/hruffieux/locus
#
# Internal core function to call the variational algorithm for probit link, optional
# fixed covariates and no external annotation variables.
# See help of `locus` function for details.
#
locus_probit_core_ <- function(Y, X, Z, list_hyper, gam_vb, mu_alpha_vb,
                               mu_beta_vb, sig2_alpha_vb, sig2_beta_vb, tol,
                               maxit, verbose, batch = "y", full_output = FALSE,
                               debug = FALSE) {

  # an intercept must be present in Z (column of ones), and X and Z must be standardized (except intercept in Z)

  d <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  q <- ncol(Z)

  with(list_hyper, { # list_init not used with the with() function to avoid
    # copy-on-write for large objects

    m2_alpha <- update_m2_alpha_(mu_alpha_vb, sig2_alpha_vb, sweep = TRUE)
    m1_beta <- update_m1_beta_(gam_vb, mu_beta_vb)
    m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb)

    mat_x_m1 <- update_mat_x_m1_(X, m1_beta)
    mat_z_mu <- update_mat_z_mu_(Z, mu_alpha_vb)

    W <- update_W_probit_(Y, mat_z_mu, mat_x_m1)

    rs_gam <- rowSums(gam_vb)
    sum_gam <- sum(rs_gam)

    phi_vb <- update_phi_z_vb_(phi, d)

    converged <- FALSE
    lb_new <- -Inf
    it <- 0

    while ((!converged) & (it <= maxit)) {

      lb_old <- lb_new
      it <- it + 1

      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("Iteration ", format(it), "... \n", sep = ""))

      # % #
      xi_vb <- update_xi_bin_vb_(xi, m2_alpha)

      zeta2_inv_vb <- phi_vb / xi_vb
      # % #

      # % #
      lambda_vb <- update_lambda_vb_(lambda, sum_gam)
      nu_vb <- update_nu_bin_vb_(nu, m2_beta)

      sig2_inv_vb <- lambda_vb / nu_vb
      # % #

      sig2_alpha_vb <- update_sig2_alpha_vb_(n, zeta2_inv_vb, intercept = TRUE)
      sig2_beta_vb <- update_sig2_beta_vb_(n, sig2_inv_vb)

      log_sig2_inv_vb <- update_log_sig2_inv_vb_(lambda_vb, nu_vb)

      digam_sum <- digamma(a + b + d)


      # different possible batch-coordinate ascent schemes:

      if (batch == "y") { # optimal scheme

        log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
        log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

        for (i in 1:q) {

          mat_z_mu <- mat_z_mu - tcrossprod(Z[, i], mu_alpha_vb[i, ])

          mu_alpha_vb[i, ] <- sig2_alpha_vb[i] * crossprod(W  - mat_z_mu - mat_x_m1, Z[, i])

          mat_z_mu <- mat_z_mu + tcrossprod(Z[, i], mu_alpha_vb[i, ])

        }

        # C++ Eigen call for expensive updates
        coreProbitLoop(X, W, gam_vb, log_om_vb, log_1_min_om_vb, log_sig2_inv_vb,
                       m1_beta, mat_x_m1, mat_z_mu, mu_beta_vb, sig2_beta_vb)

        rs_gam <- rowSums(gam_vb)

      } else if (batch == "x") {  # used internally for testing purposes,
                                  # convergence not ensured as ELBO not batch-concave

        log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
        log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

        for (k in 1:d) {

          mu_alpha_vb[, k] <- sig2_alpha_vb * (crossprod(W[, k]  - mat_z_mu[, k] - mat_x_m1[, k], Z) +  (n - 1) * mu_alpha_vb[, k])
          mu_alpha_vb[1, k] <- mu_alpha_vb[1, k] + sig2_alpha_vb[1] * mu_alpha_vb[1, k] # correction for the intercept (sums to 1)

          mat_z_mu[, k] <- Z %*% mu_alpha_vb[, k]

          mu_beta_vb[, k] <- sig2_beta_vb * (crossprod(W[, k] -  mat_z_mu[, k] - mat_x_m1[, k], X) + (n - 1) * m1_beta[, k])

          gam_vb[, k] <- exp(-log_one_plus_exp_(log_1_min_om_vb - log_om_vb -
                                                  log_sig2_inv_vb / 2 -
                                                  mu_beta_vb[, k] ^ 2 / (2 * sig2_beta_vb) -
                                                  log(sig2_beta_vb) / 2))

          m1_beta[, k] <- mu_beta_vb[, k] * gam_vb[, k]

          mat_x_m1[, k] <- X %*% m1_beta[, k]

        }

        rs_gam <- rowSums(gam_vb)

      } else if (batch == "x-y") { # used internally for testing purposes,
                                   # convergence not ensured as ELBO not batch-concave

        log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
        log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

        mu_alpha_vb <- sweep(crossprod(Z, W  - mat_z_mu - mat_x_m1) +  (n - 1) * mu_alpha_vb, 1, sig2_alpha_vb, `*`)
        mu_alpha_vb[1, ] <- mu_alpha_vb[1, ] + sig2_alpha_vb[1] * mu_alpha_vb[1, ] # correction for the intercept (sums to 1)

        mat_z_mu <- Z %*% mu_alpha_vb

        # C++ Eigen call for expensive updates
        coreProbitBatch(X, W, gam_vb, log_om_vb, log_1_min_om_vb, log_sig2_inv_vb,
                        m1_beta, mat_x_m1, mat_z_mu, mu_beta_vb, sig2_beta_vb)

        rs_gam <- rowSums(gam_vb)

      } else if (batch == "0") { # no batch, used only internally

        for (k in 1:d) {

          log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
          log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

          for (i in 1:q) {

            mat_z_mu[, k] <- mat_z_mu[, k] - Z[, i] * mu_alpha_vb[i, k]

            mu_alpha_vb[i, k] <- sig2_alpha_vb[i] * crossprod(Z[, i], W[, k] - mat_z_mu[, k] - mat_x_m1[, k])

            mat_z_mu[, k] <- mat_z_mu[, k] + Z[, i] * mu_alpha_vb[i, k]

          }

          for (j in 1:p) {

            mat_x_m1[, k] <- mat_x_m1[, k] - X[, j] * m1_beta[j, k]

            mu_beta_vb[j, k] <- sig2_beta_vb * crossprod(W[, k] - mat_x_m1[, k] - mat_z_mu[, k], X[, j])

            gam_vb[j, k] <- exp(-log_one_plus_exp_(log_1_min_om_vb[j] - log_om_vb[j] -
                                                    log_sig2_inv_vb / 2 -
                                                    mu_beta_vb[j, k] ^ 2 / (2 * sig2_beta_vb) -
                                                    log(sig2_beta_vb) / 2))


            m1_beta[j, k] <- mu_beta_vb[j, k] * gam_vb[j, k]

            mat_x_m1[, k] <- mat_x_m1[, k] + X[, j] * m1_beta[j, k]

          }

          rs_gam <- rowSums(gam_vb)

        }

      } else {

        stop ("Batch scheme not defined. Exit.")

      }

      sum_gam <- sum(rs_gam)

      m2_alpha <- update_m2_alpha_(mu_alpha_vb, sig2_alpha_vb, sweep = TRUE)
      m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb)

      W <- update_W_probit_(Y, mat_z_mu, mat_x_m1)

      a_vb <- update_a_vb(a, rs_gam)
      b_vb <- update_b_vb(b, d, rs_gam)
      om_vb <- a_vb / (a_vb + b_vb)


      lb_new <- lower_bound_probit_(Y, X, Z, a, a_vb, b, b_vb, gam_vb, lambda,
                                    nu, phi, phi_vb, sig2_alpha_vb, sig2_beta_vb,
                                    sig2_inv_vb, xi, zeta2_inv_vb, mu_alpha_vb,
                                    m1_beta, m2_alpha, m2_beta, mat_x_m1,
                                    mat_z_mu, sum_gam)

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
      create_named_list_(a, a_vb, b, b_vb, gam_vb, lambda, nu, phi, phi_vb,
                         sig2_alpha_vb, sig2_beta_vb, sig2_inv_vb, xi,
                         zeta2_inv_vb, mu_alpha_vb, m1_beta, m2_alpha, m2_beta,
                         sum_gam)
    } else {

      names_x <- colnames(X)
      names_y <- colnames(Y)
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
# lower bound (ELBO) corresponding to the `locus_probit_core` algorithm.
#
lower_bound_probit_ <- function(Y, X, Z, a, a_vb, b, b_vb, gam_vb, lambda, nu,
                                phi, phi_vb, sig2_alpha_vb, sig2_beta_vb,
                                sig2_inv_vb, xi, zeta2_inv_vb, mu_alpha_vb,
                                m1_beta, m2_alpha, m2_beta, mat_x_m1, mat_z_mu,
                                sum_gam) {

  lambda_vb <- update_lambda_vb_(lambda, sum_gam)
  nu_vb <- update_nu_bin_vb_(nu, m2_beta)

  xi_vb <- update_xi_bin_vb_(xi, m2_alpha)

  log_sig2_inv_vb <- digamma(lambda_vb) - log(nu_vb)
  log_zeta2_inv_vb <- digamma(phi_vb) - log(xi_vb)
  log_om_vb <- digamma(a_vb) - digamma(a_vb + b_vb)
  log_1_min_om_vb <- digamma(b_vb) - digamma(a_vb + b_vb)

  eps <- .Machine$double.eps # to control the argument of the log when gamma is very small

  U <- mat_x_m1 + mat_z_mu

  A <- sum(Y * pnorm(U, log.p = TRUE) +
             sweep((1 - Y) * pnorm(U, lower.tail = FALSE, log.p = TRUE), 1, Z^2 %*% sig2_alpha_vb / 2, `-`) -
             X^2 %*% (m2_beta - m1_beta^2) / 2)

  B <- sum(gam_vb * log_sig2_inv_vb / 2 - m2_beta * sig2_inv_vb / 2 +
             sweep(gam_vb, 1, log_om_vb, `*`) +
             sweep(1 - gam_vb, 1, log_1_min_om_vb, `*`) +
             gam_vb * (log(sig2_beta_vb) + 1) / 2 -
             gam_vb * log(gam_vb + eps) - (1 - gam_vb) * log(1 - gam_vb + eps))

  G <- sum((lambda - lambda_vb) * log_sig2_inv_vb - (nu - nu_vb) * sig2_inv_vb +
             lambda * log(nu) - lambda_vb * log(nu_vb) - lgamma(lambda) +
             lgamma(lambda_vb))

  H <- sum((a - a_vb) * log_om_vb + (b - b_vb) * log_1_min_om_vb - lbeta(a, b) +
             lbeta(a_vb, b_vb))

  J <- sum(sweep(-sweep(m2_alpha, 1, zeta2_inv_vb, `*`), 1,
                 log_zeta2_inv_vb + log(sig2_alpha_vb), `+`) + 1) / 2

  K <- sum((phi - phi_vb) * log_zeta2_inv_vb - (xi - xi_vb) * zeta2_inv_vb +
             phi * log(xi) - phi_vb * log(xi_vb) - lgamma(phi) + lgamma(phi_vb))

  A + B + G + H + J + K

}

