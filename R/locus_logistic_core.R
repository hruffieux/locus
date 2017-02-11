locus_logit_core_ <- function(Y, X, Z, list_hyper, chi_vb, gam_vb, mu_alpha_vb,
                            mu_beta_vb, sig2_alpha_vb, sig2_beta_vb, tol, maxit,
                            batch, verbose, full_output = FALSE) {

  d <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  q <- ncol(Z)

  # 1/2 must have been substracted from Y, and X must be standardized.

  with(list_hyper, { # list_init not used with the with() function to avoid
                     # copy-on-write for large objects

    m2_alpha <- (sig2_alpha_vb + mu_alpha_vb ^ 2)
    m1_beta <- mu_beta_vb * gam_vb
    m2_beta <- (sig2_beta_vb + mu_beta_vb ^ 2) * gam_vb

    psi_vb <- update_psi_logit_vb_(chi_vb)

    rowsums_gam <- rowSums(gam_vb)

    lambda_vb <- nu_vb <- a_vb <- b_vb <- phi_vb <- xi_vb <- NULL

    converged <- FALSE
    lb_old <- -Inf
    it <- 1

    phi_vb <- update_phi_logit_vb_(phi)

    while ((!converged) & (it <= maxit)) {

      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("Iteration ", format(it), "... \n", sep = ""))

      # % #
      xi_vb <- update_xi_logit_vb_(xi, m2_alpha)

      zeta2_inv_vb <- phi_vb / xi_vb
      # % #

      # % #
      lambda_vb <- update_lambda_logit_vb_(lambda, gam_vb)
      nu_vb <- update_nu_logit_vb_(nu, m2_beta)

      sig2_inv_vb <- lambda_vb / nu_vb
      # % #

      sig2_alpha_vb <- 1 / (2 * crossprod(Z ^ 2, psi_vb) + zeta2_inv_vb)
      sig2_beta_vb <- 1 / sweep(2 * crossprod(X ^ 2, psi_vb), 2, sig2_inv_vb, `+`)

      log_sig2_inv_vb <- digamma(lambda_vb) - log(nu_vb)

      vec_part_digam <- digamma(a + b + d)

      if (batch) { # some updates are made batch-wise

        log_om_vb <- digamma(a + rowsums_gam) - vec_part_digam
        log_1_min_om_vb <- digamma(b - rowsums_gam + d) - vec_part_digam

        mat_z_mu <-  Z %*% mu_alpha_vb
        mat_x_m1 <-  X %*% m1_beta

        for (i in 1:q) {
          mat_z_mu <- mat_z_mu - tcrossprod(Z[, i], mu_alpha_vb[i, ])

          mu_alpha_vb[i, ] <- sig2_alpha_vb[i, ] * crossprod(Y - 2 * psi_vb * (mat_z_mu + mat_x_m1), Z[, i])

          mat_z_mu <- mat_z_mu + tcrossprod(Z[, i], mu_alpha_vb[i, ])
        }

        for (j in 1:p) {
          mat_x_m1 <- mat_x_m1 - tcrossprod(X[, j], m1_beta[j, ])

          mu_beta_vb[j, ] <- sig2_beta_vb[j, ] * crossprod(Y - 2 * psi_vb * (mat_x_m1 + mat_z_mu), X[, j])

          log_part_gam_vb <- log_om_vb[j] + log(sig2_beta_vb[j, ]) / 2 +
            mu_beta_vb[j, ] ^ 2 / (2 * sig2_beta_vb[j, ])

          log_part2_gam_vb <- log_1_min_om_vb[j] - log_sig2_inv_vb / 2

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

          vec_z_i_k <-  Z %*% mu_alpha_vb[, k]
          vec_x_j_k <-  X %*% m1_beta[, k]

          for (i in 1:q) {
            vec_z_i_k <- vec_z_i_k - Z[, i] * mu_alpha_vb[i, k]

            mu_alpha_vb[i, k] <- sig2_alpha_vb[i, k] * crossprod(Z[, i], Y[, k] - 2 * psi_vb[, k] * (vec_z_i_k + vec_x_j_k))

            vec_z_i_k <- vec_z_i_k + Z[, i] * mu_alpha_vb[i, k]
          }

          for (j in 1:p) {

            vec_x_j_k <- vec_x_j_k - X[, j] * m1_beta[j, k]

            mu_beta_vb[j, k] <- sig2_beta_vb[j, k] * crossprod(X[, j], Y[, k] - 2 * psi_vb[, k] * (vec_z_i_k + vec_x_j_k))

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

      m2_alpha <- (sig2_alpha_vb + mu_alpha_vb ^ 2)
      m2_beta <- (sig2_beta_vb + mu_beta_vb ^ 2) * gam_vb

      a_vb <- a + rowsums_gam
      b_vb <- b - rowsums_gam + d
      om_vb <- a_vb / (a_vb + b_vb)

      # % #
      if (!batch) {
        mat_x_m1 <- X %*% m1_beta
        mat_z_mu <- Z %*% mu_alpha_vb
      }

      chi_vb <- sqrt(X^2 %*% m2_beta + mat_x_m1^2 - X^2 %*% m1_beta^2 +
                       Z^2 %*% sig2_alpha_vb + mat_z_mu^2 +
                       2 * mat_x_m1 * mat_z_mu)

      psi_vb <- update_psi_logit_vb_(chi_vb)


      # % #

      lb_new <- lower_bound_logit_(Y, X, Z, a, a_vb, b, b_vb, chi_vb, gam_vb,
                                 lambda, nu, phi, phi_vb, psi_vb, sig2_alpha_vb,
                                 sig2_beta_vb, sig2_inv_vb, xi, zeta2_inv_vb,
                                 mu_alpha_vb, m1_beta, m2_alpha, m2_beta, mat_x_m1,
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
                  format(lb_new), ". \n\n",
                  sep = ""))
      } else {
        cat("Maximal number of iterations reached before convergence. Exit.")
      }
    }

    lb_opt <- lb_new

    if (full_output) { # for internal use only
      create_named_list_(a, a_vb, b, b_vb, chi_vb, gam_vb, lambda, nu, phi,
                         phi_vb, psi_vb, sig2_alpha_vb, sig2_beta_vb,
                         sig2_inv_vb, xi, zeta2_inv_vb, mu_alpha_vb, m1_beta,
                         m2_alpha, m2_beta, mat_x_m1, mat_z_mu)
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


update_phi_logit_vb_ <- function(phi) {

  phi + 1 / 2

}

update_xi_logit_vb_ <- function(xi, m2_alpha) {

  xi + m2_alpha / 2

}

update_psi_logit_vb_ <- function(chi_vb) {

  sig <- function(chi) {
    1 / (1 + exp(-chi))
  }

  (sig(chi_vb) - 1 / 2) / (2 * chi_vb)

}

update_lambda_logit_vb_ <- function(lambda, gam_vb) {

  lambda + colSums(gam_vb) / 2

}

update_nu_logit_vb_ <- function(nu, m2_beta) {

  nu + colSums(m2_beta) / 2

}

lower_bound_logit_ <- function(Y, X, Z, a, a_vb, b, b_vb, chi_vb, gam_vb,
                             lambda, nu, phi, phi_vb, psi_vb, sig2_alpha_vb,
                             sig2_beta_vb, sig2_inv_vb, xi, zeta2_inv_vb,
                             mu_alpha_vb, m1_beta, m2_alpha, m2_beta, mat_x_m1,
                             mat_z_mu) {

  lambda_vb <- update_lambda_logit_vb_(lambda, gam_vb)
  nu_vb <- update_nu_logit_vb_(nu, m2_beta)

  xi_vb <- update_xi_logit_vb_(xi, m2_alpha)

  log_sig2_inv_vb <- digamma(lambda_vb) - log(nu_vb)
  log_zeta2_inv_vb <- digamma(phi_vb) - log(xi_vb)
  log_om_vb <- digamma(a_vb) - digamma(a_vb + b_vb)
  log_1_min_om_vb <- digamma(b_vb) - digamma(a_vb + b_vb)


  A <- sum(chi_vb / 2 - log_sum_exp_mat_(list(-chi_vb/2, chi_vb/2)) +
             Y * (mat_x_m1 + mat_z_mu)  -  chi_vb / 2 -
             psi_vb * (X^2 %*% m2_beta + mat_x_m1^2 - X^2 %*% m1_beta^2 +
                         Z^2 %*% m2_alpha + mat_z_mu^2 - Z^2 %*% mu_alpha_vb^2 +
                         2 * mat_x_m1 * mat_z_mu - chi_vb^2))

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

  K <- 1 / 2 * sum(log_zeta2_inv_vb - zeta2_inv_vb * m2_alpha + log(sig2_alpha_vb) + 1)

  L <- sum((phi - phi_vb) * log_zeta2_inv_vb - (xi - xi_vb) * zeta2_inv_vb +
             phi * log(xi) - phi_vb * log(xi_vb) - lgamma(phi) +
             lgamma(phi_vb))

  A + B + H + J + K + L

}

