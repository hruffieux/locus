# This file is part of the `locus` R package:
#     https://github.com/hruffieux/locus
#
# Internal functions gathering the variational updates for the core algorithms.
# Besides improving code readability via modular programming, the main purpose
# is to avoid copy-and-paste programming, as most of these updates (or slightly
# modified versions) are used more than once in the different core algorithms.
# For this reason, we choose to create functions for most variational updates,
# even for those consisting in very basic operations.
# Note that we don't modularize the body of the core for loops for performance
# reasons.


#####################
## alpha's updates ##
#####################

update_m2_alpha_ <- function(alpha_vb, sig2_alpha_vb, sweep = FALSE) {

  if(sweep) {

    sweep(alpha_vb ^ 2, 1, sig2_alpha_vb, `+`)

  } else {

    sig2_alpha_vb + alpha_vb ^ 2
  }

}


update_sig2_alpha_vb_ <- function(n, zeta2_inv_vb, tau_vb = NULL, intercept = FALSE, c = 1) {

  den <- n - 1 + zeta2_inv_vb

  if (intercept)
    den[1] <- den[1] + 1 # the first column of Z was not scaled, it is the intercept.

  if (is.null(tau_vb)) {

    1 / (c * den)

  } else {

    1 / (c * tcrossprod(den, as.matrix(tau_vb)))

  }

}


update_sig2_alpha_logit_vb_ <- function(Z, psi_vb, zeta2_inv_vb) {

  1 / sweep(2 * crossprod(Z ^ 2, psi_vb), 1, zeta2_inv_vb, `+`)

}


update_mat_z_mu_ <- function(Z, alpha_vb) Z %*% alpha_vb


####################
## beta's updates ##
####################

update_beta_vb_ <- function(gam_vb, mu_beta_vb) gam_vb * mu_beta_vb


update_g_beta_vb_ <- function(list_mu_beta_vb, gam_vb) {

  G <- length(list_mu_beta_vb)

  lapply(1:G, function(g) sweep(list_mu_beta_vb[[g]], 2, gam_vb[g, ], `*`))

}


update_m2_beta_ <- function(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = FALSE) {

  if(sweep) {

    sweep(mu_beta_vb ^ 2, 2, sig2_beta_vb, `+`) * gam_vb

  } else {

    (mu_beta_vb ^ 2 + sig2_beta_vb) * gam_vb

  }

}


update_sig2_beta_vb_ <- function(n, sig2_inv_vb, tau_vb = NULL, c = 1) {

  if(is.null(tau_vb)) {

    1 / (c * (n - 1 + sig2_inv_vb))

  } else {

    1 / (c * (n - 1 + sig2_inv_vb) * tau_vb)

  }
}


update_sig2_beta_logit_vb_ <- function(X, psi_vb, sig2_inv_vb) {

  1 / (2 * crossprod(X ^ 2, psi_vb) + sig2_inv_vb)

}


update_mat_x_m1_ <- function(X, beta_vb) X %*% beta_vb


update_g_mat_x_m1_ <- function(list_X, list_beta_vb) {

  G <- length(list_X)

  Reduce(`+`, lapply(1:G, function(g) list_X[[g]] %*% list_beta_vb[[g]]))
}


update_g_m1_btb_ <- function(gam_vb, list_mu_beta_vb, list_sig2_beta_star, tau_vb) { ## not list_sig2_beta_star_inv!

  d <- length(tau_vb)
  G <- length(list_mu_beta_vb)

  lapply(1:G, function(g) {
    gam_vb[g, ]^2 * colSums(list_mu_beta_vb[[g]]^2) + # colSums(A^2) = diag(crossprod(A))
      gam_vb[g, ] * (sum(diag(as.matrix(list_sig2_beta_star[[g]]))) / tau_vb +
                       sapply(1:d, function(k) (1-gam_vb[g, k]) * sum(list_mu_beta_vb[[g]][, k]^2))) # tr(mu_gt mu_gt^T) = sum(mu_gt^2)
  })

}


update_g_m1_btXtXb_ <- function(list_X, gam_vb, list_mu_beta_vb, list_sig2_beta_star, tau_vb) {

  d <- length(tau_vb)
  G <- length(list_mu_beta_vb)

  lapply(1:G, function(g) {
    gam_vb[g, ]^2 * colSums((list_X[[g]] %*% list_mu_beta_vb[[g]])^2) +
      gam_vb[g, ] * (sum(crossprod(list_X[[g]]) * list_sig2_beta_star[[g]]) / tau_vb +
                       sapply(1:d, function(k) (1-gam_vb[g, k]) * sum(crossprod(list_X[[g]]) * tcrossprod(list_mu_beta_vb[[g]][, k])))) # tr(AB^T) = sum_ij A_ij B_ij
  })
}


########################
## c0 and c's updates ##
########################

update_sig2_c0_vb_ <- function(d, s02, c = 1) 1 / (c * (d + (1/s02)))


###################
## chi's updates ##
###################

update_chi_vb_ <- function(X, Z, beta_vb, m2_beta, mat_x_m1, mat_z_mu, sig2_alpha_vb) {

  sqrt(X^2 %*% m2_beta + mat_x_m1^2 - X^2 %*% beta_vb^2 + Z^2 %*% sig2_alpha_vb +
         mat_z_mu^2 + 2 * mat_x_m1 * mat_z_mu)
}


update_psi_logit_vb_ <- function(chi_vb) {

  exp(log(exp(log_sigmoid_(chi_vb)) - 1 / 2) - log(2 * chi_vb))

}


#####################
## omega's updates ##
#####################

a_vb <- update_a_vb <- function(a, rs_gam, c = 1) c * (a + rs_gam) - c + 1


b_vb <- update_b_vb <- function(b, d, rs_gam, c = 1) c * (b - rs_gam + d) - c + 1


update_log_om_vb <- function(a, digam_sum, rs_gam, c = 1) digamma(c * (a + rs_gam) - c + 1) - digam_sum


update_log_1_min_om_vb <- function(b, d, digam_sum, rs_gam, c = 1) digamma(c * (b - rs_gam + d) - c + 1) - digam_sum



#####################
## sigma's updates ##
#####################

update_lambda_vb_ <- function(lambda, sum_gam, c = 1) c * (lambda + sum_gam / 2) - c + 1


update_g_lambda_vb_ <- function(lambda, g_sizes, rs_gam) lambda + sum(g_sizes * rs_gam) / 2


update_nu_vb_ <- function(nu, m2_beta, tau_vb, c = 1) c * as.numeric(nu + crossprod(tau_vb, colSums(m2_beta)) / 2)


update_g_nu_vb_ <- function(nu, list_m1_btb, tau_vb) nu + sum(tau_vb * Reduce(`+`, list_m1_btb))/2


update_nu_bin_vb_ <- function(nu, m2_beta) nu + sum(m2_beta) / 2


update_log_sig2_inv_vb_ <- function(lambda_vb, nu_vb) digamma(lambda_vb) - log(nu_vb)


###################
## tau's updates ##
###################

update_eta_vb_ <- function(n, eta, gam_vb, c = 1) c * (eta + n / 2 + colSums(gam_vb) / 2) - c + 1


update_g_eta_vb_ <- function(n, eta, g_sizes, gam_vb) eta + n / 2 + as.numeric(crossprod(gam_vb, g_sizes)) / 2


update_eta_z_vb_ <- function(n, q, eta, gam_vb, c = 1) c * (eta + n / 2 + colSums(gam_vb) / 2 + q / 2) - c + 1


update_kappa_vb_ <- function(Y, kappa, mat_x_m1, beta_vb, m2_beta, sig2_inv_vb, c = 1) {

  n <- nrow(Y)

  c * (kappa + (colSums(Y^2) - 2 * colSums(Y * mat_x_m1)  +
                  (n - 1 + sig2_inv_vb) * colSums(m2_beta) +
                  colSums(mat_x_m1^2) - (n - 1) * colSums(beta_vb^2))/ 2)

}


update_g_kappa_vb_ <- function(Y, list_X, kappa, list_beta_vb, list_m1_btb,
                               list_m1_btXtXb, mat_x_m1, sig2_inv_vb) {

  n <- nrow(Y)
  G <- length(list_beta_vb)

  # avoid using do.call() as can trigger node stack overflow
  kappa + (colSums(Y^2) - 2 * colSums(Y * mat_x_m1)  +
             Reduce(`+`, list_m1_btXtXb) +
             sig2_inv_vb * Reduce(`+`, list_m1_btb) +
             colSums(mat_x_m1^2) -
             Reduce(`+`, lapply(1:G, function(g) colSums((list_X[[g]] %*% list_beta_vb[[g]])^2) ))) / 2
}


update_kappa_z_vb_ <- function(Y, Z, kappa, alpha_vb, beta_vb, m2_alpha,
                               m2_beta, mat_x_m1, mat_z_mu, sig2_inv_vb,
                               zeta2_inv_vb, intercept = FALSE, c = 1) {
  n <- nrow(Y)

  kappa_vb <- c * (kappa + (colSums(Y^2) - 2 * colSums(Y * (mat_x_m1 + mat_z_mu))  +
                         (n - 1 + sig2_inv_vb) * colSums(m2_beta) +
                         colSums(mat_x_m1^2) - (n - 1) * colSums(beta_vb^2) +
                         (n - 1) * colSums(m2_alpha) +
                         crossprod(m2_alpha, zeta2_inv_vb) +
                         colSums(mat_z_mu^2) - (n - 1) * colSums(alpha_vb^2) +
                         2 * colSums(mat_x_m1 * mat_z_mu))/ 2)

  if (intercept)
    kappa_vb <- kappa_vb + c * (m2_alpha[1, ] - (alpha_vb[1, ])^2) / 2

  kappa_vb
}


update_log_tau_vb_ <- function(eta_vb, kappa_vb) digamma(eta_vb) - log(kappa_vb)


#####################
## theta's updates ##
#####################

update_theta_vb_ <- function(W, m0, S0_inv, sig2_theta_vb, vec_fac_st,
                                mat_add = 0, is_mat = FALSE, c = 1) {

  if (is.null(vec_fac_st)) {

    if (is_mat) {

      theta_vb <- c * sig2_theta_vb * (rowSums(W) + S0_inv * m0 - rowSums(mat_add)) 

    } else {

      theta_vb <- c * sig2_theta_vb * (rowSums(W) + S0_inv * m0 - sum(mat_add)) 

    }


  } else {

    if (c != 1)
      stop("Annealing not implemented when Sigma_0 is not the identity matrix.")

    bl_ids <- unique(vec_fac_st)
    n_bl <- length(bl_ids)

    if (is_mat) {

      theta_vb <- unlist(lapply(1:n_bl, function(bl) {
        sig2_theta_vb[[bl]] %*% (rowSums(W[vec_fac_st == bl_ids[bl], , drop = FALSE]) +
                                   S0_inv[[bl]] %*% m0[vec_fac_st == bl_ids[bl]] -
                                   rowSums(mat_add[vec_fac_st == bl_ids[bl], , drop = FALSE]))  # mat_add = sweep(mat_v_mu, 1, theta_vb, `-`)
      }))
    } else {

      theta_vb <- unlist(lapply(1:n_bl, function(bl) {
        sig2_theta_vb[[bl]] %*% (rowSums(W[vec_fac_st == bl_ids[bl], , drop = FALSE]) +
                                   S0_inv[[bl]] %*% m0[vec_fac_st == bl_ids[bl]] -
                                   sum(mat_add)) 
      }))
    }

  }

}


update_sig2_theta_vb_ <- function(d, p, list_struct, s02, X = NULL, c = 1) {

  if (is.null(list_struct)) {

    vec_fac_st <- NULL

    S0_inv <- 1 / s02 # stands for a diagonal matrix of size p with this value on the (constant) diagonal
    sig2_theta_vb <- as.numeric(update_sig2_c0_vb_(d, s02, c = c)) # idem

    vec_sum_log_det_theta <- - p * (log(s02) + log(d + S0_inv))

  } else {

    if (c != 1)
      stop("Annealing not implemented when Sigma_0 is not the identity matrix.")

    if (is.null(X))
      stop("X must be passed to the update_sig2_theta_function.")

    vec_fac_st <- list_struct$vec_fac_st
    n_cpus <- list_struct$n_cpus

    S0_inv <- parallel::mclapply(unique(vec_fac_st), function(bl) {

      corX <- cor(X[, vec_fac_st == bl, drop = FALSE])
      corX <- as.matrix(Matrix::nearPD(corX, corr = TRUE, do2eigen = TRUE)$mat) # regularization in case of non-positive definiteness.

      as.matrix(solve(corX) / s02)
    }, mc.cores = n_cpus)

    if (is.list(S0_inv)) {

      sig2_theta_vb <- parallel::mclapply(S0_inv, function(mat) {
        as.matrix(solve(mat + diag(d, nrow(mat))))
      }, mc.cores = n_cpus)

    } else {

      sig2_theta_vb <- 1 / (S0_inv + d)

    }

    vec_sum_log_det_theta <- log_det(S0_inv) + log_det(sig2_theta_vb) # vec_sum_log_det_theta[bl] = log(det(S0_inv_bl)) + log(det(sig2_theta_vb_bl))

  }

  create_named_list_(S0_inv, sig2_theta_vb, vec_sum_log_det_theta, vec_fac_st)
}



#################
## W's updates ##
#################

update_W_probit_ <- function(Y, mat_z_mu, mat_x_m1) {
  
  mat_z_mu + mat_x_m1 + inv_mills_ratio_matrix_(Y, mat_z_mu + mat_x_m1)

}


update_W_struct_ <- function(gam_vb, theta_vb) {

  log_pnorm <- pnorm(theta_vb, log.p = TRUE)
  log_1_pnorm <- pnorm(theta_vb, log.p = TRUE, lower.tail = FALSE)
  
  imr0 <- inv_mills_ratio_(0, theta_vb, log_1_pnorm, log_pnorm)
  
  sweep(sweep(gam_vb, 1, (inv_mills_ratio_(1, theta_vb, log_1_pnorm, log_pnorm) - imr0), `*`),
        1,  theta_vb + imr0, `+`)

}


####################
## zeta's updates ##
####################


update_phi_z_vb_ <- function(phi, d, c = 1) c * (phi + d / 2) - c + 1


update_xi_z_vb_ <- function(xi, tau_vb, m2_alpha, c = 1) c * (xi + m2_alpha %*% tau_vb / 2)


update_xi_bin_vb_ <- function(xi, m2_alpha) xi + rowSums(m2_alpha) / 2


update_log_zeta2_inv_vb_ <- function(phi_vb, xi_vb) digamma(phi_vb) - log(xi_vb)
