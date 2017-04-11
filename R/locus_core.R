# This file is part of the `locus` R package:
#     https://github.com/hruffieux/locus
#
# Internal core function to call the variational algorithm for identity link, no
# fixed covariates and no external annotation variables.
# See help of `locus` function for details.
#
locus_core_ <- function(Y, X, list_hyper, gam_vb, mu_beta_vb, sig2_beta_vb,
                        tau_vb, tol, maxit, verbose, batch = "y",
                        full_output = FALSE, debug = FALSE) {

  # Y must have been centered, and X, standardized.

  d <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)

  with(list_hyper, { # list_init not used with the with() function to avoid
                     # copy-on-write for large objects

    m1_beta <- update_m1_beta_(gam_vb, mu_beta_vb)
    m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE)

    mat_x_m1 <- update_mat_x_m1_(X, m1_beta)

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
      lambda_vb <- update_lambda_vb_(lambda, sum_gam)
      nu_vb <- update_nu_vb_(nu, m2_beta, tau_vb)

      sig2_inv_vb <- lambda_vb / nu_vb
      # % #

      # % #
      eta_vb <- update_eta_vb_(n, eta, gam_vb)
      kappa_vb <- update_kappa_vb_(Y, X, kappa, mat_x_m1, m1_beta, m2_beta, sig2_inv_vb)

      tau_vb <- eta_vb / kappa_vb
      # % #

      sig2_beta_vb <- update_sig2_beta_vb_(n, sig2_inv_vb, tau_vb)

      log_tau_vb <- update_log_tau_vb_(eta_vb, kappa_vb)
      log_sig2_inv_vb <- update_log_sig2_inv_vb_(lambda_vb, nu_vb)

      digam_sum <- digamma(a + b + d)


      # different possible batch-coordinate ascent schemes:

      if (batch == "y") { # optimal scheme

        log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
        log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

        # C++ Eigen call for expensive updates
        coreLoop(X, Y, gam_vb, log_om_vb, log_1_min_om_vb, log_sig2_inv_vb,
                 log_tau_vb, m1_beta, mat_x_m1, mu_beta_vb, sig2_beta_vb, tau_vb)


        rs_gam <- rowSums(gam_vb)

      } else if (batch == "x") { # used only internally, convergence not ensured

        log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
        log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

        for (k in 1:d) {

          mu_beta_vb[, k] <- sig2_beta_vb[k] * tau_vb[k] *
            (crossprod(Y[, k] - mat_x_m1[, k], X) + (n - 1) * m1_beta[, k])


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

        # C++ Eigen call for expensive updates
        coreBatch(X, Y, gam_vb, log_om_vb, log_1_min_om_vb, log_sig2_inv_vb,
                  log_tau_vb, m1_beta, mat_x_m1, mu_beta_vb, sig2_beta_vb, tau_vb)

        rs_gam <- rowSums(gam_vb)

      } else if (batch == "0") { # no batch, used only internally

        for (k in 1:d) {

          log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
          log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

          for (j in 1:p) {

            mat_x_m1[, k] <- mat_x_m1[, k] - X[, j] * m1_beta[j, k]

            mu_beta_vb[j, k] <- sig2_beta_vb[k] * tau_vb[k] * crossprod(Y[, k] - mat_x_m1[, k], X[, j])

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

      m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE)

      a_vb <- update_a_vb(a, rs_gam)
      b_vb <- update_b_vb(b, d, rs_gam)
      om_vb <- a_vb / (a_vb + b_vb)

      sum_gam <- sum(rs_gam)

      lb_new <- lower_bound_(Y, X, a, a_vb, b, b_vb, eta, gam_vb, kappa, lambda,
                             nu, sig2_beta_vb, sig2_inv_vb, tau_vb, m1_beta,
                             m2_beta, mat_x_m1, sum_gam)

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
      create_named_list_(a, a_vb, b, b_vb, eta, gam_vb, kappa, lambda,
                         nu, sig2_beta_vb, sig2_inv_vb, tau_vb, m1_beta,
                         m2_beta, mat_x_m1, sum_gam)
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
lower_bound_ <- function(Y, X, a, a_vb, b, b_vb, eta, gam_vb, kappa, lambda, nu,
                         sig2_beta_vb, sig2_inv_vb, tau_vb, m1_beta, m2_beta,
                         mat_x_m1, sum_gam) {

  n <- nrow(Y)

  eta_vb <- update_eta_vb_(n, eta, gam_vb)
  kappa_vb <- update_kappa_vb_(Y, X, kappa, mat_x_m1, m1_beta, m2_beta, sig2_inv_vb)

  lambda_vb <- update_lambda_vb_(lambda, sum_gam)
  nu_vb <- update_nu_vb_(nu, m2_beta, tau_vb)

  log_tau_vb <- digamma(eta_vb) - log(kappa_vb)
  log_sig2_inv_vb <- digamma(lambda_vb) - log(nu_vb)
  log_om_vb <- digamma(a_vb) - digamma(a_vb + b_vb)
  log_1_min_om_vb <- digamma(b_vb) - digamma(a_vb + b_vb)

  A <- sum(-n / 2 * log(2 * pi) + n / 2 * log_tau_vb -
             tau_vb * (kappa_vb - colSums(m2_beta) * sig2_inv_vb / 2 - kappa))

  eps <- .Machine$double.eps # to control the argument of the log when gamma is very small
  B <- sum(log_sig2_inv_vb * gam_vb / 2 +
             sweep(gam_vb, 2, log_tau_vb, `*`) / 2 -
             sweep(m2_beta, 2, tau_vb, `*`) * sig2_inv_vb / 2 +
             sweep(gam_vb, 1, log_om_vb, `*`) +
             sweep(1 - gam_vb, 1, log_1_min_om_vb, `*`) +
             1 / 2 * sweep(gam_vb, 2, log(sig2_beta_vb) + 1, `*`) -
             gam_vb * log(gam_vb + eps) - (1 - gam_vb) * log(1 - gam_vb + eps))

  G <- sum((eta - eta_vb) * log_tau_vb -
             (kappa - kappa_vb) * tau_vb + eta * log(kappa) -
             eta_vb * log(kappa_vb) - lgamma(eta) + lgamma(eta_vb))

  H <- (lambda - lambda_vb) * log_sig2_inv_vb - (nu - nu_vb) * sig2_inv_vb +
    lambda * log(nu) - lambda_vb * log(nu_vb) - lgamma(lambda) +
    lgamma(lambda_vb)

  J <- sum((a - a_vb) * log_om_vb + (b - b_vb) * log_1_min_om_vb - lbeta(a, b) +
             lbeta(a_vb, b_vb))

  A + B + G + H + J

}

