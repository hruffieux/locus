# This file is part of the `locus` R package:
#     https://github.com/hruffieux/locus
#
# Internal core function to call the variational algorithm for logit link,
# optional fixed covariates and no external annotation variables.
# See help of `locus` function for details.
#
locus_logit_core_ <- function(Y, X, Z, list_hyper, chi_vb, gam_vb, alpha_vb,
                              mu_beta_vb, sig2_alpha_vb, sig2_beta_vb, tol,
                              maxit, verbose, batch = "y", full_output = FALSE,
                              debug = FALSE) {

  # 1/2 must have been substracted from Y and X must have been standardized (except intercept in Z)

  d <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  q <- ncol(Z)

  with(list_hyper, { # list_init not used with the with() function to avoid
    # copy-on-write for large objects

    m2_alpha <- update_m2_alpha_(alpha_vb, sig2_alpha_vb)
    beta_vb <- update_beta_vb_(gam_vb, mu_beta_vb)
    m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb)

    mat_x_m1 <- update_mat_x_m1_(X, beta_vb)
    mat_z_mu <- update_mat_z_mu_(Z, alpha_vb)

    phi_vb <- update_phi_z_vb_(phi, d)

    rs_gam <- rowSums(gam_vb)
    sum_gam <- sum(rs_gam)
    digam_sum <- digamma(a + b + d)

    converged <- FALSE
    lb_new <- -Inf
    it <- 0

    while ((!converged) & (it < maxit)) {

      lb_old <- lb_new
      it <- it + 1

      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste0("Iteration ", format(it), "... \n"))

      # % #
      chi_vb <- update_chi_vb_(X, Z, beta_vb, m2_beta, mat_x_m1, mat_z_mu, sig2_alpha_vb)

      psi_vb <- update_psi_logit_vb_(chi_vb)
      # % #

      # % #
      xi_vb <- update_xi_bin_vb_(xi, m2_alpha)

      zeta2_inv_vb <- phi_vb / xi_vb
      # % #

      # % #
      lambda_vb <- update_lambda_vb_(lambda, sum_gam)
      nu_vb <- update_nu_bin_vb_(nu, m2_beta)

      sig2_inv_vb <- lambda_vb / nu_vb
      # % #

      sig2_alpha_vb <- update_sig2_alpha_logit_vb_(Z, psi_vb, zeta2_inv_vb)
      sig2_beta_vb <- update_sig2_beta_logit_vb_(X, psi_vb, sig2_inv_vb)

      log_sig2_inv_vb <- update_log_sig2_inv_vb_(lambda_vb, nu_vb)

      # different possible batch-coordinate ascent schemes:

      if (batch == "y") { # optimal scheme

        log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
        log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

        for (i in sample(1:q)) {
          mat_z_mu <- mat_z_mu - tcrossprod(Z[, i], alpha_vb[i, ])

          alpha_vb[i, ] <- sig2_alpha_vb[i, ] * crossprod(Y - 2 * psi_vb * (mat_z_mu + mat_x_m1), Z[, i])

          mat_z_mu <- mat_z_mu + tcrossprod(Z[, i], alpha_vb[i, ])
        }

        # C++ Eigen call for expensive updates
        shuffled_ind <- as.numeric(sample(0:(p-1))) # Zero-based index in C++

        coreLogitLoop(X, Y, gam_vb, log_om_vb, log_1_min_om_vb, log_sig2_inv_vb,
                      beta_vb, mat_x_m1, mat_z_mu, mu_beta_vb, psi_vb,
                      sig2_beta_vb, shuffled_ind)

        rs_gam <- rowSums(gam_vb)

      } else if (batch == "x") {  # used internally for testing purposes,
        # convergence not ensured as ELBO not batch-concave

        log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
        log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

        for (k in sample(1:d)) {

          alpha_vb[, k] <- sig2_alpha_vb[, k] * (crossprod(Y[, k] - 2 * psi_vb[, k] * (mat_z_mu[, k] + mat_x_m1[, k]), Z) +
                                                      colSums(sweep(sweep(Z, 2, alpha_vb[, k], `*`), 1, 2 * psi_vb[, k], `*`) * Z))

          mat_z_mu[, k] <- Z %*% alpha_vb[, k]

          mu_beta_vb[, k] <- sig2_beta_vb[, k] * (crossprod(Y[, k] -  2 * psi_vb[, k] * (mat_z_mu[, k] + mat_x_m1[, k]), X) +
                                                    colSums(sweep(sweep(X, 2, beta_vb[, k], `*`), 1, 2 * psi_vb[, k], `*`) * X))


          gam_vb[, k] <- exp(-log_one_plus_exp_(log_1_min_om_vb - log_om_vb -
                                                  log_sig2_inv_vb / 2 -
                                                  mu_beta_vb[, k] ^ 2 / (2 * sig2_beta_vb[, k]) -
                                                  log(sig2_beta_vb[, k]) / 2))

          beta_vb[, k] <- mu_beta_vb[, k] * gam_vb[, k]

          mat_x_m1[, k] <- X %*% beta_vb[, k]

        }

        rs_gam <- rowSums(gam_vb)

      } else if (batch == "x-y") {  # used internally for testing purposes,
        # convergence not ensured as ELBO not batch-concave

        log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
        log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

        alpha_vb <- sig2_alpha_vb * (crossprod(Z, Y - 2 * psi_vb * (mat_z_mu + mat_x_m1)) +
                                          sapply(1:d, function(k)
                                            colSums(sweep(sweep(Z, 2, alpha_vb[, k], `*`), 1, 2 * psi_vb[, k], `*`) * Z)))

        mat_z_mu <- Z %*% alpha_vb

        mu_beta_vb <- sig2_beta_vb * (crossprod(X, Y -  2 * psi_vb * (mat_z_mu + mat_x_m1)) +
                                        sapply(1:d, function(k)
                                          colSums(sweep(sweep(X, 2, beta_vb[, k], `*`), 1, 2 * psi_vb[, k], `*`) * X)))


        gam_vb <- exp(-log_one_plus_exp_(log_1_min_om_vb - log_om_vb -
                                           log_sig2_inv_vb / 2 -
                                           mu_beta_vb ^ 2 / (2 * sig2_beta_vb) -
                                           log(sig2_beta_vb) / 2))

        beta_vb<- mu_beta_vb * gam_vb

        mat_x_m1 <- X %*% beta_vb


        rs_gam <- rowSums(gam_vb)

      } else if (batch == "0"){ # no batch, used only internally


        for (k in sample(1:d)) {

          log_om_vb <- update_log_om_vb(a, digam_sum, rs_gam)
          log_1_min_om_vb <- update_log_1_min_om_vb(b, d, digam_sum, rs_gam)

          for (i in sample(1:q)) {

            mat_z_mu[, k] <- mat_z_mu[, k] - Z[, i] * alpha_vb[i, k]

            alpha_vb[i, k] <- sig2_alpha_vb[i, k] *
              crossprod(Y[, k] - 2 * psi_vb[, k] * (mat_z_mu[, k] + mat_x_m1[, k]), Z[, i])

            mat_z_mu[, k] <- mat_z_mu[, k] + Z[, i] * alpha_vb[i, k]

          }

          for (j in sample(1:p)) {

            mat_x_m1[, k] <- mat_x_m1[, k] - X[, j] * beta_vb[j, k]

            mu_beta_vb[j, k] <- sig2_beta_vb[j, k] *
              crossprod(Y[, k] - 2 * psi_vb[, k] * (mat_z_mu[, k] + mat_x_m1[, k]), X[, j])

            gam_vb[j, k] <- exp(-log_one_plus_exp_(log_1_min_om_vb[j] - log_om_vb[j] -
                                                     log_sig2_inv_vb / 2 -
                                                     mu_beta_vb[j, k] ^ 2 / (2 * sig2_beta_vb[j, k]) -
                                                     log(sig2_beta_vb[j, k]) / 2))

            beta_vb[j, k] <- mu_beta_vb[j, k] * gam_vb[j, k]

            mat_x_m1[, k] <- mat_x_m1[, k] + X[, j] * beta_vb[j, k]

          }

          rs_gam <- rowSums(gam_vb)

        }

      } else {

        stop ("Batch scheme not defined. Exit.")

      }

      sum_gam <- sum(rs_gam)

      m2_alpha <- update_m2_alpha_(alpha_vb, sig2_alpha_vb)
      m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb)

      a_vb <- update_a_vb(a, rs_gam)
      b_vb <- update_b_vb(b, d, rs_gam)
      om_vb <- a_vb / (a_vb + b_vb)

      lb_new <- elbo_logit_(Y, X, Z, a, a_vb, b, b_vb, beta_vb, chi_vb, gam_vb,
                            lambda, nu, phi, phi_vb, psi_vb, sig2_alpha_vb,
                            sig2_beta_vb, sig2_inv_vb, xi, zeta2_inv_vb,
                            alpha_vb, m2_alpha, m2_beta, mat_x_m1,
                            mat_z_mu, sum_gam)

      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste0("ELBO = ", format(lb_new), "\n\n"))

      if (debug && lb_new < lb_old)
        stop("ELBO not increasing monotonically. Exit. ")

      converged <- (abs(lb_new - lb_old) < tol)

    }

    if (verbose) {
      if (converged) {
        cat(paste0("Convergence obtained after ", format(it), " iterations. \n",
                  "Optimal marginal log-likelihood variational lower bound ",
                  "(ELBO) = ", format(lb_new), ". \n\n"))
      } else {
        warning("Maximal number of iterations reached before convergence. Exit.")
      }
    }

    lb_opt <- lb_new

    names_x <- colnames(X)
    names_y <- colnames(Y)
    names_z <- colnames(Z)
    
    rownames(gam_vb) <- rownames(beta_vb) <- names_x
    colnames(gam_vb) <- colnames(beta_vb) <- names_y
    
    names(om_vb) <- names_x
    rownames(alpha_vb) <- names_z
    colnames(alpha_vb) <- names_y
    
    diff_lb <- abs(lb_opt - lb_old)
    
    if (full_output) { # for internal use only
      
      create_named_list_(a, a_vb, b, b_vb, beta_vb, chi_vb, gam_vb, lambda, 
                         mu_beta_vb, nu, om_vb, phi, phi_vb, psi_vb, sig2_alpha_vb, 
                         sig2_beta_vb, sig2_inv_vb, xi, zeta2_inv_vb, alpha_vb, 
                         m2_alpha, m2_beta, mat_x_m1, mat_z_mu, sum_gam, 
                         converged, it, lb_opt, diff_lb)
    } else {

      create_named_list_(beta_vb, gam_vb, om_vb, alpha_vb, converged, it, lb_opt, diff_lb)
    }
  })

}


# Internal function which implements the marginal log-likelihood variational
# lower bound (ELBO) corresponding to the `locus_logit_core` algorithm.
#
elbo_logit_ <- function(Y, X, Z, a, a_vb, b, b_vb, beta_vb, chi_vb, gam_vb, 
                        lambda, nu, phi, phi_vb, psi_vb, sig2_alpha_vb, 
                        sig2_beta_vb, sig2_inv_vb, xi, zeta2_inv_vb, alpha_vb,
                        m2_alpha, m2_beta, mat_x_m1, mat_z_mu, sum_gam) {

  lambda_vb <- update_lambda_vb_(lambda, sum_gam)
  nu_vb <- update_nu_bin_vb_(nu, m2_beta)

  xi_vb <- update_xi_bin_vb_(xi, m2_alpha)

  log_sig2_inv_vb <- digamma(lambda_vb) - log(nu_vb)
  log_zeta2_inv_vb <- digamma(phi_vb) - log(xi_vb)
  log_om_vb <- digamma(a_vb) - digamma(a_vb + b_vb)
  log_1_min_om_vb <- digamma(b_vb) - digamma(a_vb + b_vb)

  elbo_A <- e_y_logit_(X, Y, Z, chi_vb, beta_vb, m2_alpha, m2_beta, mat_x_m1,
                       mat_z_mu, alpha_vb, psi_vb)

  elbo_B <- e_beta_gamma_bin_(gam_vb, log_om_vb, log_1_min_om_vb,
                              log_sig2_inv_vb, m2_beta, sig2_beta_vb, sig2_inv_vb)

  elbo_C <- e_sig2_inv_(lambda, lambda_vb, log_sig2_inv_vb, nu, nu_vb, sig2_inv_vb)

  elbo_D <- e_omega_(a, a_vb, b, b_vb, log_om_vb, log_1_min_om_vb)

  elbo_E <- e_alpha_logit_(m2_alpha, log_zeta2_inv_vb, sig2_alpha_vb, zeta2_inv_vb)

  elbo_F <- e_zeta2_inv_(log_zeta2_inv_vb, phi, phi_vb, xi, xi_vb, zeta2_inv_vb)


  elbo_A + elbo_B + elbo_C + elbo_D + elbo_E + elbo_F

}

