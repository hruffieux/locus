locus_core_ <- function(Y, X, list_hyper, gam_vb, mu_beta_vb, sig2_beta_vb,
                        tau_vb, tol, maxit, batch, verbose, full_output = FALSE) {

  # Y must have been centered, and X, standardized.

  d <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)

  with(list_hyper, { # list_init not used with the with() function to avoid
                     # copy-on-write for large objects
    m1_beta <- mu_beta_vb * gam_vb
    m2_beta <- sweep(mu_beta_vb ^ 2, 2, sig2_beta_vb, `+`) * gam_vb

    mat_x_m1 <-  X %*% m1_beta

    rowsums_gam <- rowSums(gam_vb)
    sum_gam <- sum(rowsums_gam)

    lambda_vb <- nu_vb <- eta_vb <- kappa_vb <- a_vb <- b_vb <- NULL

    converged <- FALSE
    lb_old <- -Inf
    it <- 1

    while ((!converged) & (it <= maxit)) {

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

      sig2_beta_vb <- 1 / ((n - 1 + sig2_inv_vb) * tau_vb)

      log_tau_vb <- digamma(eta_vb) - log(kappa_vb)
      log_sig2_inv_vb <- digamma(lambda_vb) - log(nu_vb)

      vec_part_digam <- digamma(a + b + d)

      if (batch) { # some updates are made batch-wise

        log_om_vb <- digamma(a + rowsums_gam) - vec_part_digam
        log_1_min_om_vb <- digamma(b - rowsums_gam + d) - vec_part_digam

        for (j in 1:p) {
          mat_x_m1 <- mat_x_m1 - tcrossprod(X[, j], m1_beta[j, ])

          mu_beta_vb[j, ] <- sig2_beta_vb * (tau_vb *
                                               crossprod(Y - mat_x_m1, X[, j]))

          log_part_gam_vb <- log_om_vb[j] + log(sig2_beta_vb) / 2 +
            mu_beta_vb[j, ] ^ 2 / (2 * sig2_beta_vb)

          log_part2_gam_vb <- log_1_min_om_vb[j] - log_tau_vb / 2 -
            log_sig2_inv_vb / 2

          gam_vb[j, ] <- exp(log_part_gam_vb -
                               log_sum_exp_mat_(list(log_part_gam_vb, log_part2_gam_vb)))

          m1_beta[j, ] <- mu_beta_vb[j, ] * gam_vb[j, ]

          mat_x_m1 <- mat_x_m1 + tcrossprod(X[, j], m1_beta[j, ])
        }

        rowsums_gam <- rowSums(gam_vb)

      } else {

        for (k in 1:d) {

          log_om_vb <- digamma(a + rowsums_gam) - vec_part_digam
          log_1_min_om_vb <- digamma(b - rowsums_gam + d) - vec_part_digam

          vec_x_j_k <-  X %*% m1_beta[, k]
          for (j in 1:p) {

            vec_x_j_k <- vec_x_j_k - X[, j] * m1_beta[j, k]

            mu_beta_vb[j, k] <- sig2_beta_vb[k] * tau_vb[k] *
              crossprod(X[, j], Y[,k] - vec_x_j_k)

            log_part_gam_vb <- log_om_vb[j] + log(sig2_beta_vb[k]) / 2 +
              mu_beta_vb[j, k] ^ 2 / (2 * sig2_beta_vb[k])

            log_part2_gam_vb <- log_1_min_om_vb[j] - log_tau_vb[k] / 2 -
              log_sig2_inv_vb / 2

            gam_vb[j, k] <- exp(log_part_gam_vb -
                                  log_sum_exp_(c(log_part_gam_vb, log_part2_gam_vb)))

            m1_beta[j, k] <- mu_beta_vb[j, k] * gam_vb[j, k]

            vec_x_j_k <- vec_x_j_k + X[, j] * m1_beta[j, k]

          }

          rowsums_gam <- rowSums(gam_vb)

        }

      }

      m2_beta <- sweep(mu_beta_vb ^ 2, 2, sig2_beta_vb, `+`) * gam_vb

      if (!batch)
        mat_x_m1 <-  X %*% m1_beta

      a_vb <- a + rowsums_gam
      b_vb <- b - rowsums_gam + d
      om_vb <- a_vb / (a_vb + b_vb)

      sum_gam <- sum(rowsums_gam)

      lb_new <- lower_bound_(Y, X, a, a_vb, b, b_vb, eta, gam_vb, kappa, lambda,
                             nu, sig2_beta_vb, sig2_inv_vb, tau_vb, m1_beta,
                             m2_beta, mat_x_m1, sum_gam)


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
        cat("Maximal number of iterations reached before convergence. Exit.")
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

      create_named_list_(lb_opt, gam_vb, om_vb)
    }
  })

}



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

