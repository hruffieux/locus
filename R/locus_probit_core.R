locus_probit_core_ <- function(Y, X, Z, list_hyper, gam_vb, mu_alpha_vb,
                               mu_beta_vb, sig2_alpha_vb, sig2_beta_vb, tol,
                               maxit, batch, verbose, full_output = FALSE) {

  # an intercept must be present in Z (column of ones), and X and Z must be standardized (except intercept in Z)

  d <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  q <- ncol(Z)

  with(list_hyper, { # list_init not used with the with() function to avoid
    # copy-on-write for large objects

    m2_alpha <- sweep(mu_alpha_vb ^ 2, 1, sig2_alpha_vb, `+`)
    m1_beta <- mu_beta_vb * gam_vb
    m2_beta <- (mu_beta_vb ^ 2 + sig2_beta_vb) * gam_vb

    mat_z_mu <-  Z %*% mu_alpha_vb
    mat_x_m1 <-  X %*% m1_beta

    W <- update_W_probit_(Y, mat_z_mu, mat_x_m1)

    rowsums_gam <- rowSums(gam_vb)
    sum_gam <- sum(rowsums_gam)

    phi_vb <- update_phi_z_vb_(phi, d)

    converged <- FALSE
    lb_old <- -Inf
    it <- 1

    while ((!converged) & (it <= maxit)) {

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

      sig2_alpha_vb <- n - 1 + zeta2_inv_vb
      sig2_alpha_vb[1] <- sig2_alpha_vb[1] + 1 # the first column of Z was not scaled, it is the intercept.
      sig2_alpha_vb <- 1 / sig2_alpha_vb
      sig2_beta_vb <- 1 / (n - 1 + sig2_inv_vb)

      log_sig2_inv_vb <- digamma(lambda_vb) - log(nu_vb)

      vec_part_digam <- digamma(a + b + d)

      if (batch) { # some updates are made batch-wise

        log_om_vb <- digamma(a + rowsums_gam) - vec_part_digam
        log_1_min_om_vb <- digamma(b - rowsums_gam + d) - vec_part_digam

        for (i in 1:q) {

          mat_z_mu <- mat_z_mu - tcrossprod(Z[, i], mu_alpha_vb[i, ])

          mu_alpha_vb[i, ] <- sig2_alpha_vb[i] * crossprod(W  - mat_z_mu - mat_x_m1, Z[, i])

          mat_z_mu <- mat_z_mu + tcrossprod(Z[, i], mu_alpha_vb[i, ])

        }

        for (j in 1:p) {

          mat_x_m1 <- mat_x_m1 - tcrossprod(X[, j], m1_beta[j, ])

          mu_beta_vb[j,] <- sig2_beta_vb * crossprod(W - mat_x_m1 - mat_z_mu, X[, j])

          log_part_gam_vb <- log_om_vb[j] + log(sig2_beta_vb) / 2 +
            mu_beta_vb[j, ] ^ 2 / (2 * sig2_beta_vb)

          log_part2_gam_vb <- log_1_min_om_vb[j] - rep(log_sig2_inv_vb, d) / 2

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

          vec_z_r_k <-  Z %*% mu_alpha_vb[, k]
          vec_x_j_k <-  X %*% m1_beta[, k]

          for (i in 1:q) {

            vec_z_r_k <- vec_z_r_k - Z[, i] * mu_alpha_vb[i, k]

            mu_alpha_vb[i, k] <- sig2_alpha_vb[i] * crossprod(Z[, i], W[, k] - vec_z_r_k - vec_x_j_k)

            vec_z_r_k <- vec_z_r_k + Z[, i] * mu_alpha_vb[i, k]
          }

          for (j in 1:p) {

            vec_x_j_k <- vec_x_j_k - X[, j] * m1_beta[j, k]

            mu_beta_vb[j, k] <- sig2_beta_vb * crossprod(W[, k] - vec_z_r_k - vec_x_j_k, X[, j])

            log_part_gam_vb <- log_om_vb[j] + log(sig2_beta_vb) / 2 +
              mu_beta_vb[j, k] ^ 2 / (2 * sig2_beta_vb)

            log_part2_gam_vb <- log_1_min_om_vb[j] - log_sig2_inv_vb / 2

            gam_vb[j, k] <- exp(log_part_gam_vb -
                                  log_sum_exp_(c(log_part_gam_vb, log_part2_gam_vb)))

            m1_beta[j, k] <- mu_beta_vb[j, k] * gam_vb[j, k]

            vec_x_j_k <- vec_x_j_k + X[, j] * m1_beta[j, k]

          }

          rowsums_gam <- rowSums(gam_vb)

        }

      }

      sum_gam <- sum(rowsums_gam)

      m2_alpha <- sweep(mu_alpha_vb ^ 2, 1, sig2_alpha_vb, `+`)
      m2_beta <- (mu_beta_vb ^ 2 + sig2_beta_vb) * gam_vb

      if (!batch) {
        mat_z_mu <-  Z %*% mu_alpha_vb
        mat_x_m1 <-  X %*% m1_beta
      }

      W <- update_W_probit_(Y, mat_z_mu, mat_x_m1)

      a_vb <- a + rowsums_gam
      b_vb <- b - rowsums_gam + d
      om_vb <- a_vb / (a_vb + b_vb)


      lb_new <- lower_bound_probit_(Y, X, Z, a, a_vb, b, b_vb, gam_vb, lambda,
                                    nu, phi, phi_vb, sig2_alpha_vb, sig2_beta_vb,
                                    sig2_inv_vb, xi, zeta2_inv_vb, mu_alpha_vb,
                                    m1_beta, m2_alpha, m2_beta, mat_x_m1,
                                    mat_z_mu, sum_gam)

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
      create_named_list_(a, a_vb, b, b_vb, gam_vb, lambda, nu, phi, phi_vb,
                         sig2_alpha_vb, sig2_beta_vb, sig2_inv_vb, xi,
                         zeta2_inv_vb, mu_alpha_vb, m1_beta, m2_alpha, m2_beta,
                         mat_x_m1, mat_z_mu, sum_gam)
    } else {

      names_x <- colnames(X)
      names_y <- colnames(Y)
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

