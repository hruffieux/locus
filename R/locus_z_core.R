## with covariates.
locus_z_core_ <- function(Y, X, Z, d, n, p, q, list_hyper, gam_vb, mu_beta_vb,
                          sig2_beta_vb, tau_vb, mu_alpha_vb, sig2_alpha_vb, tol,
                          maxit, batch, verbose, full_output = FALSE) {

  # Y must have been centered, and X, standardized.

  with(list_hyper, {  # list_init not used with the with() function to avoid
                      # copy-on-write for large objects
    m2_alpha <- (sig2_alpha_vb + mu_alpha_vb ^ 2)

    m1_beta <- mu_beta_vb * gam_vb
    m2_beta <- sweep(mu_beta_vb ^ 2, 2, sig2_beta_vb, `+`) * gam_vb

    rowsums_gam <- rowSums(gam_vb)
    sum_gam <- sum(rowsums_gam)

    lambda_vb <- nu_vb <- eta_vb <- kappa_vb <- a_vb <- b_vb <- phi_vb <- xi_vb <- NULL ###

    converged <- FALSE
    lb_old <- -Inf
    it <- 1

    phi_vb <- update_phi_z_vb_(phi, d) ###

    while ((!converged) & (it <= maxit)) {

      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("Iteration ", format(it), "... \n", sep = ""))

      # % #
      xi_vb <- update_xi_z_vb_(xi, tau_vb, m2_alpha) ###

      zeta2_inv_vb <- phi_vb / xi_vb
      # % #

      # % #
      lambda_vb <- update_lambda_vb_(sum_gam, lambda)
      nu_vb <- update_nu_vb_(tau_vb, m2_beta, nu)

      sig2_inv_vb <- lambda_vb / nu_vb
      # % #

      # % #
      eta_vb <- update_eta_z_vb_(gam_vb, eta, q, n)
      kappa_vb <- update_kappa_z_vb_(Y, X, Z, d, n, p, q, sig2_inv_vb, zeta2_inv_vb,
                                     mu_alpha_vb, m2_alpha, m1_beta, m2_beta, kappa)

      tau_vb <- eta_vb / kappa_vb
      # % #

      sig2_alpha_vb <- 1 / tcrossprod(n - 1 + zeta2_inv_vb, tau_vb)
      sig2_beta_vb <- 1 / ((n - 1 + sig2_inv_vb) * tau_vb)

      log_tau_vb <- digamma(eta_vb) - log(kappa_vb)
      log_sig2_inv_vb <- digamma(lambda_vb) - log(nu_vb)

      vec_part_digam <- digamma(a + b + d)


      if (batch) {
        log_om_vb <- digamma(a + rowsums_gam) - vec_part_digam
        log_1_min_om_vb <- digamma(b - rowsums_gam + d) - vec_part_digam

        mat_z_mu_i <-  Z %*% mu_alpha_vb
        mat_x_m1_j <-  X %*% m1_beta

        for (i in 1:q) {
          mat_z_mu_i <- mat_z_mu_i - tcrossprod(Z[, i], mu_alpha_vb[i, ])

          mu_alpha_vb[i, ] <- sig2_alpha_vb[i, ] * (tau_vb *
                                                      crossprod(Y  - mat_z_mu_i - mat_x_m1_j, Z[, i]))

          mat_z_mu_i <- mat_z_mu_i + tcrossprod(Z[, i], mu_alpha_vb[i, ])
        }

        for (j in 1:p) {

          mat_x_m1_j <- mat_x_m1_j -  tcrossprod(X[, j], m1_beta[j, ])

          mu_beta_vb[j,] <- sig2_beta_vb * (tau_vb *
                                              crossprod(Y - mat_x_m1_j - mat_z_mu_i, X[, j]))

          log_part_gam_vb <- log_om_vb[j] + log(sig2_beta_vb) / 2 +
            mu_beta_vb[j, ] ^ 2 / (2 * sig2_beta_vb)

          log_part2_gam_vb <- log_1_min_om_vb[j] - log_tau_vb / 2 -
            log_sig2_inv_vb / 2

          gam_vb[j, ] <- exp(log_part_gam_vb - log_sum_exp_vec_(list(log_part_gam_vb, log_part2_gam_vb)))
          m1_beta[j, ] <- mu_beta_vb[j, ] * gam_vb[j, ]

          mat_x_m1_j <- mat_x_m1_j +  tcrossprod(X[, j], m1_beta[j, ])

        }

        rowsums_gam <- rowSums(gam_vb)

      } else {

        for (k in 1:d) {

          vec_z_i_k <-  Z %*% mu_alpha_vb[, k]
          vec_x_j_k <-  X %*% m1_beta[, k]

          for (i in 1:q) {

            vec_z_i_k <- vec_z_i_k - Z[, i] * mu_alpha_vb[i, k]

            mu_alpha_vb[i, k] <- sig2_alpha_vb[i, k] * tau_vb[k] *
              crossprod(Z[, i], Y[,k]  - vec_z_i_k - vec_x_j_k)

            vec_z_i_k <- vec_z_i_k + Z[, i] * mu_alpha_vb[i, k]
          }

          log_om_vb <- digamma(a + rowsums_gam) - vec_part_digam
          log_1_min_om_vb <- digamma(b - rowsums_gam + d) - vec_part_digam


          for (j in 1:p) {

            vec_x_j_k <- vec_x_j_k - X[, j] * m1_beta[j, k]

            mu_beta_vb[j, k] <- sig2_beta_vb[k] * tau_vb[k] *
              crossprod(X[, j], Y[,k] - vec_x_j_k - vec_z_i_k)

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

      m2_alpha <- (sig2_alpha_vb + mu_alpha_vb ^ 2) ####
      m2_beta <- sweep(mu_beta_vb ^ 2, 2, sig2_beta_vb, `+`) * gam_vb

      a_vb <- a + rowsums_gam
      b_vb <- b - rowsums_gam + d
      om_vb <- a_vb / (a_vb + b_vb)

      sum_gam <- sum(rowsums_gam)

      lb_new <- lower_bound_z_(Y, X, Z, d, n, p, q, mu_alpha_vb, sig2_alpha_vb,
                               zeta2_inv_vb, sig2_beta_vb, sig2_inv_vb, tau_vb,
                               gam_vb, eta, kappa, lambda, nu, a, b, a_vb, b_vb,
                               phi, phi_vb, xi, m2_alpha, m1_beta, m2_beta, sum_gam)


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
                  format(lb_new), ". \n\n",
                  sep = ""))
      } else {
        cat("Maximal number of iterations reached before convergence. Exit.")
      }
    }

    lb_opt <- lb_new

    x_prpnst <- rowSums(gam_vb)
    names(x_prpnst) <- colnames(X)

    y_prpnst <- colSums(gam_vb)
    names(y_prpnst) <- colnames(Y)


    if (full_output) { # for internal use only
      create_named_list_(mu_alpha_vb, sig2_alpha_vb, zeta2_inv_vb, mu_beta_vb,
                         sig2_beta_vb, sig2_inv_vb, tau_vb, gam_vb, om_vb, eta,
                         kappa, lambda, nu, a, b, a_vb, b_vb, phi, phi_vb, xi,
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
      names(x_prpnst) <- names_x
      names(y_prpnst) <- names_y

      create_named_list_(lb_opt, gam_vb, om_vb, mu_alpha_vb, x_prpnst, y_prpnst)
    }
  })
}

update_phi_z_vb_ <- function(phi, d) {

  phi + d / 2

}

update_xi_z_vb_ <- function(xi, tau_vb, m2_alpha) {

  xi + m2_alpha %*% tau_vb / 2

}

update_eta_z_vb_ <- function(gam_vb, eta, q, n) {

  eta + n / 2 + colSums(gam_vb) / 2 + q / 2

}

update_kappa_z_vb_ <- function(Y_mat, X_mat, Z_mat, d, n, p, q, sig2_inv_vb,
                               zeta2_inv_vb, mu_alpha_vb, m2_alpha, m1_beta,
                               m2_beta, kappa) {
  # put X_mat and Y_mat instead of X and Y to avoid conflicts with the function sapply,
  # which has also an "X" argument with different meaning...

  prod_alpha <- Z_mat %*% mu_alpha_vb
  prod_beta <- X_mat %*% m1_beta

  # Z must be standardized as we use (Z \hadamard Z)^T \one_n = (n-1)\one_q
  kappa_vb <- kappa + colSums(Y_mat ^ 2) / 2 - colSums(Y_mat * prod_beta) -
    colSums(Y_mat * prod_alpha) +
    colSums(prod_alpha * prod_beta) +
    (n - 1 + sig2_inv_vb) * colSums(m2_beta) / 2 +
    (n - 1) * colSums(m2_alpha) / 2 + crossprod(m2_alpha, zeta2_inv_vb) / 2

  if (q > 1) {

    mat_z_list <- lapply(q:1, function(i) {
      Z_mat[, i, drop = FALSE] %*% mu_alpha_vb[i,, drop = FALSE]
    }
    )

    cum_mat_z_list <- list(0)
    for (i in 1:(q-1)) {
      cum_mat_z_list[[i+1]] <- cum_mat_z_list[[i]] + mat_z_list[[i]]
    }
    cum_mat_z_list[[1]] <- NULL

    mix_z_sum <- sapply(q:2, function(i) {
      colSums(mat_z_list[[i]] * cum_mat_z_list[[i-1]])
    })

    if(d == 1) mix_z_sum <- t(mix_z_sum)

    kappa_vb <- kappa_vb + rowSums(mix_z_sum)
  }


  if (p > 1) {

    mat_x_list <- lapply(p:1, function(j) {
      X_mat[, j, drop = FALSE] %*% m1_beta[j,, drop = FALSE]
    }
    )

    cum_mat_x_list <- list(0)
    for (j in 1:(p-1)) {
      cum_mat_x_list[[j+1]] <- cum_mat_x_list[[j]] + mat_x_list[[j]]
    }
    cum_mat_x_list[[1]] <- NULL

    mix_x_sum <- sapply(p:2, function(j) {
      colSums(mat_x_list[[j]] * cum_mat_x_list[[j-1]])
    })

    if(d == 1) mix_x_sum <- t(mix_x_sum)

    kappa_vb <- kappa_vb + rowSums(mix_x_sum)

  }

  kappa_vb
}

lower_bound_z_ <- function(Y, X, Z, d, n, p, q, mu_alpha_vb, sig2_alpha_vb,
                           zeta2_inv_vb, sig2_beta_vb, sig2_inv_vb, tau_vb,
                           gam_vb, eta, kappa, lambda, nu, a, b, a_vb, b_vb,
                           phi, phi_vb, xi, m2_alpha, m1_beta, m2_beta, sum_gam) {

  xi_vb <- update_xi_z_vb_(xi, tau_vb, m2_alpha)

  eta_vb <- update_eta_z_vb_(gam_vb, eta, q, n)
  kappa_vb <- update_kappa_z_vb_(Y, X, Z, d, n, p, q, sig2_inv_vb, zeta2_inv_vb,
                                 mu_alpha_vb, m2_alpha, m1_beta, m2_beta, kappa)

  lambda_vb <- update_lambda_vb_(sum_gam, lambda)
  nu_vb <- update_nu_vb_(tau_vb, m2_beta, nu)

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

