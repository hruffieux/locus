locus_bin_core_ <- function(Y, X, list_hyper, chi_vb, gam_vb, mu_beta_vb,
                            sig2_beta_vb, tol, maxit, batch, verbose,
                            full_output = FALSE) {

  d <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)

  # 1/2 must have been substracted from Y, and X must be standardized.

  with(list_hyper, { # list_init not used with the with() function to avoid
                     # copy-on-write for large objects

    m1_beta <- mu_beta_vb * gam_vb
    m2_beta <- (sig2_beta_vb + mu_beta_vb ^ 2) * gam_vb

    psi_vb <- update_psi_bin_vb_(chi_vb)

    rowsums_gam <- rowSums(gam_vb)

    lambda_vb <- nu_vb <- a_vb <- b_vb <- NULL

    converged <- FALSE
    lb_old <- -Inf
    it <- 1

    while ((!converged) & (it <= maxit)) {

      #if (verbose & (it == 1 | it %% 5 == 0)
      if (verbose)
        cat(paste("Iteration ", format(it), "... \n", sep = ""))

      # % #
      lambda_vb <- update_lambda_bin_vb_(lambda, gam_vb)
      nu_vb <- update_nu_bin_vb_(nu, m2_beta)

      sig2_inv_vb <- lambda_vb / nu_vb
      # % #

      sig2_beta_vb <- 1 / sweep(2 * crossprod(X * X, psi_vb), 2, sig2_inv_vb, `+`)

      log_sig2_inv_vb <- digamma(lambda_vb) - log(nu_vb)

      vec_part_digam <- digamma(a + b + d)

      if (batch) { # some updates are made batch-wise

        log_om_vb <- digamma(a + rowsums_gam) - vec_part_digam
        log_1_min_om_vb <- digamma(b - rowsums_gam + d) - vec_part_digam

        mat_x_m1_j <-  X %*% m1_beta

        for (j in 1:p) {
          mat_x_m1_j <- mat_x_m1_j - tcrossprod(X[, j], m1_beta[j, ])

          mu_beta_vb[j, ] <- sig2_beta_vb[j, ] * crossprod(Y - 2 * psi_vb * mat_x_m1_j, X[, j])

          log_part_gam_vb <- log_om_vb[j] + log(sig2_beta_vb[j, ]) / 2 +
            mu_beta_vb[j, ] ^ 2 / (2 * sig2_beta_vb[j, ])

          log_part2_gam_vb <- log_1_min_om_vb[j] - log_sig2_inv_vb / 2

          gam_vb[j, ] <- exp(log_part_gam_vb -
                               log_sum_exp_mat_(list(log_part_gam_vb, log_part2_gam_vb)))

          m1_beta[j, ] <- mu_beta_vb[j, ] * gam_vb[j, ]

          mat_x_m1_j <- mat_x_m1_j + tcrossprod(X[, j], m1_beta[j, ])
        }

        rowsums_gam <- rowSums(gam_vb)

      } else {

        for (k in 1:d) {

          log_om_vb <- digamma(a + rowsums_gam) - vec_part_digam
          log_1_min_om_vb <- digamma(b - rowsums_gam + d) - vec_part_digam

          vec_x_j_k <-  X %*% m1_beta[, k]
          for (j in 1:p) {

            vec_x_j_k <- vec_x_j_k - X[, j] * m1_beta[j, k]

            mu_beta_vb[j, k] <- sig2_beta_vb[, k] * crossprod(X[, j], Y[, k] - 2 * psi_vb[, k] * vec_x_j_k)

            log_part_gam_vb <- log_om_vb[j] + log(sig2_beta_vb[j, k]) / 2 +
              mu_beta_vb[j, k] ^ 2 / (2 * sig2_beta_vb[j, k])

            log_part2_gam_vb <- log_1_min_om_vb[j] - log_sig2_inv_vb[k] / 2

            gam_vb[j, k] <- exp(log_part_gam_vb -
                                  log_sum_exp_(c(log_part_gam_vb, log_part2_gam_vb)))

            m1_beta[j, k] <- mu_beta_vb[j, k] * gam_vb[j, k]

            vec_x_j_k <- vec_x_j_k + X[, j] * m1_beta[j, k]

          }

          rowsums_gam <- rowSums(gam_vb)

        }

      }

      m2_beta <- (sig2_beta_vb + mu_beta_vb ^ 2) * gam_vb
      m3_beta <- sig2_beta_vb * gam_vb

      a_vb <- a + rowsums_gam
      b_vb <- b - rowsums_gam + d
      om_vb <- a_vb / (a_vb + b_vb)


      chi_vb <- sqrt(apply(m3_beta, 2, function(m3_beta_k) rowSums(sweep(X, 2, sqrt(m3_beta_k), `*`)^2) ) +
                       (X %*% (m1_beta))^2)

      psi_vb <- update_psi_bin_vb_(chi_vb)

      lb_new <- lower_bound_bin_(Y, X, a, a_vb, b, b_vb, chi_vb, gam_vb,
                                 lambda, nu, psi_vb, sig2_beta_vb, sig2_inv_vb,
                                 m1_beta, m2_beta)


      #if (verbose & (it == 1 | it %% 5 == 0))
      #  cat(paste("Lower bound = ", format(lb_new), "\n\n", sep = ""))

      if (verbose)
        cat(paste("Lower bound = ", lb_new, "\n\n", sep = ""))

      converged <- (abs(lb_new - lb_old) < tol)

      lb_old <- lb_new
      it <- it + 1
    }



    if (verbose) {
      if (converged) {
        cat(paste("Convergence obtained after ", format(it),
                  " iterations with variational lower bound = ",
                  format(lb_new), ". \n\n",
                  sep = ""))
      } else {
        cat("Maximal number of iterations reached before convergence. Exit.")
      }
    }

    lb_opt <- lb_new
    x_prpnst <- rowSums(gam_vb)
    y_prpnst <- colSums(gam_vb)


    if (full_output) { # for internal use only
      create_named_list_(a, a_vb, b, b_vb, chi_vb, gam_vb,
                         lambda, nu, psi_vb, sig2_beta_vb, sig2_inv_vb,
                         m1_beta, m2_beta)
    } else {
      names_x <- colnames(X)
      names_y <- colnames(Y)

      rownames(gam_vb) <- names_x
      colnames(gam_vb) <- names_y
      names(om_vb) <- names_x
      names(x_prpnst) <- names_x
      names(y_prpnst) <- names_y

      create_named_list_(lb_opt, gam_vb, om_vb, x_prpnst, y_prpnst)
    }
  })

}


update_psi_bin_vb_ <- function(chi_vb) {

  sig <- function(chi) {
    1 / (1 + exp(-chi))
  }

  (sig(chi_vb) - 1 / 2) / (2 * chi_vb)

}

update_lambda_bin_vb_ <- function(lambda, gam_vb) {

  lambda + colSums(gam_vb) / 2

}

update_nu_bin_vb_ <- function(nu, m2_beta) {

  nu + colSums(m2_beta) / 2

}

lower_bound_bin_ <- function(Y, X, a, a_vb, b, b_vb, chi_vb, gam_vb,
                             lambda, nu, psi_vb, sig2_beta_vb, sig2_inv_vb,
                             m1_beta, m2_beta) {

  lambda_vb <- update_lambda_bin_vb_(lambda, gam_vb)
  nu_vb <- update_nu_bin_vb_(nu, m2_beta)

  log_sig2_inv_vb <- digamma(lambda_vb) - log(nu_vb)
  log_om_vb <- digamma(a_vb) - digamma(a_vb + b_vb)
  log_1_min_om_vb <- digamma(b_vb) - digamma(a_vb + b_vb)

  X_sq <- X^2
  A <- sum(chi_vb / 2 - log_sum_exp_mat_(list(-chi_vb/2, chi_vb/2)) +
             (X %*% m1_beta) * Y -  chi_vb / 2 -
             psi_vb * (X_sq %*% m2_beta + (X %*% m1_beta)^2 -
                         X_sq %*% m1_beta^2 - chi_vb^2))

  eps <- .Machine$double.eps # to control the argument of the log when gamma is very small
  B <- sum( sweep(gam_vb, 2, log_sig2_inv_vb, `*`) / 2 -
             sweep(m2_beta, 2, sig2_inv_vb, `*`) / 2 +
             sweep(gam_vb, 1, log_om_vb, `*`) +
             sweep(1 - gam_vb, 1, log_1_min_om_vb, `*`) +
             gam_vb * (log(sig2_beta_vb) + 1) / 2 -
             gam_vb * log(gam_vb + eps) - (1 - gam_vb) * log(1 - gam_vb + eps))

  H <- sum((lambda - lambda_vb) * log_sig2_inv_vb - (nu - nu_vb) * sig2_inv_vb +
             lambda * log(nu) - lambda_vb * log(nu_vb) - lgamma(lambda) +
             lgamma(lambda_vb))

  J <- sum((a - a_vb) * log_om_vb + (b - b_vb) * log_1_min_om_vb - lbeta(a, b) +
             lbeta(a_vb, b_vb))

  A + B + H + J

}

