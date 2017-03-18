## These functions gather the variational updates for the core methods.
## Besides improving code readability via modular programming, the main purpose
## is to avoid copy-and-paste programming, as most of these updates (or slightly
## modified versions) are used more than once in the different core algorithms.
## For this reason, we choose to create functions for most variational updates,
## even for those consisting in very basic operations.
## Note that we don't modularize the body of the core for loops for performance
## reasons.


#####################
## alpha's updates ##
#####################

update_m2_alpha_ <- function(mu_alpha_vb, sig2_alpha_vb, sweep = FALSE) {

  if(sweep) {

    sweep(mu_alpha_vb ^ 2, 1, sig2_alpha_vb, `+`)

  } else {

    sig2_alpha_vb + mu_alpha_vb ^ 2
  }

}

update_sig2_alpha_vb_ <- function(n, zeta2_inv_vb, tau_vb = NULL, intercept = FALSE) {

  den <- n - 1 + zeta2_inv_vb

  if (intercept)
    den[1] <- den[1] + 1 # the first column of Z was not scaled, it is the intercept.

  if (is.null(tau_vb)) {

    1 / den

  } else {

    1 / tcrossprod(den, as.matrix(tau_vb))

  }

}


update_sig2_alpha_logit_vb_ <- function(Z, psi_vb, zeta2_inv_vb) {

  1 / sweep(2 * crossprod(Z ^ 2, psi_vb), 1, zeta2_inv_vb, `+`)

}

update_mat_z_mu_ <- function(Z, mu_alpha_vb) Z %*% mu_alpha_vb


####################
## beta's updates ##
####################

update_m1_beta_ <- function(gam_vb, mu_beta_vb) gam_vb * mu_beta_vb

update_m2_beta_ <- function(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = FALSE) {

  if(sweep) {

    sweep(mu_beta_vb ^ 2, 2, sig2_beta_vb, `+`) * gam_vb

  } else {

    (mu_beta_vb ^ 2 + sig2_beta_vb) * gam_vb

  }

}

update_sig2_beta_vb_ <- function(n, sig2_inv_vb, tau_vb = NULL) {

  if(is.null(tau_vb)) {

    1 / (n - 1 + sig2_inv_vb)

  } else {

    1 / ((n - 1 + sig2_inv_vb) * tau_vb)

  }
}

update_sig2_beta_logit_vb_ <- function(X, psi_vb, sig2_inv_vb) {

  1 / (2 * crossprod(X ^ 2, psi_vb) + sig2_inv_vb)

}


update_mat_x_m1_ <- function(X, m1_beta) X %*% m1_beta


########################
## c0 and c's updates ##
########################

update_mu_c0_vb_ <- function(W, mat_v_mu, m0, s02, sig2_c0_vb) sig2_c0_vb * (rowSums(W - mat_v_mu) + m0 / s02)

update_sig2_c0_vb_ <- function(d, s02) 1 / (d + (1/s02))

update_sig2_c_vb_ <- function(p, s2) 1 / (p - 1 + (1/s2))

update_mat_v_mu_ <- function(V, mu_c0_vb, mu_c_vb) sweep(V %*% mu_c_vb, 1, mu_c0_vb, `+`)


#####################
## omega's updates ##
#####################


a_vb <- update_a_vb <- function(a, rs_gam) a + rs_gam

b_vb <- update_b_vb <- function(b, d, rs_gam) b - rs_gam + d

update_log_om_vb <- function(a, digam_sum, rs_gam) digamma(a + rs_gam) - digam_sum

update_log_1_min_om_vb <- function(b, d, digam_sum, rs_gam) digamma(b - rs_gam + d) - digam_sum


###################
## chi's updates ##
###################

update_chi_vb_ <- function(X, Z, m1_beta, m2_beta, mat_x_m1, mat_z_mu, sig2_alpha_vb) {

  sqrt(X^2 %*% m2_beta + mat_x_m1^2 - X^2 %*% m1_beta^2 + Z^2 %*% sig2_alpha_vb +
         mat_z_mu^2 + 2 * mat_x_m1 * mat_z_mu)
}

log_sigmoid <- function(chi) {

  - log(1 + exp(- chi)) # chi is always positive so no overflow possible (underflow neither, thanks to the "+1")

}

update_psi_logit_vb_ <- function(chi_vb) {

  exp(log(exp(log_sigmoid(chi_vb)) - 1 / 2) - log(2 * chi_vb))

}


#####################
## sigma's updates ##
#####################

update_lambda_vb_ <- function(lambda, sum_gam) lambda + sum_gam / 2

update_nu_vb_ <- function(nu, m2_beta, tau_vb) as.numeric(nu + crossprod(tau_vb, colSums(m2_beta)) / 2)

update_nu_bin_vb_ <- function(nu, m2_beta) nu + sum(m2_beta) / 2


update_log_sig2_inv_vb_ <- function(lambda_vb, nu_vb) digamma(lambda_vb) - log(nu_vb)


###################
## tau's updates ##
###################

update_eta_vb_ <- function(n, eta, gam_vb) eta + n / 2 + colSums(gam_vb) / 2

update_eta_z_vb_ <- function(n, q, eta, gam_vb) eta + n / 2 + colSums(gam_vb) / 2 + q / 2

update_kappa_vb_ <- function(Y, X, kappa, mat_x_m1, m1_beta, m2_beta, sig2_inv_vb) {

  n <- nrow(Y)

  kappa + (colSums(Y^2) - 2 * colSums(Y * mat_x_m1)  +
             (n - 1 + sig2_inv_vb) * colSums(m2_beta) +
             colSums(mat_x_m1^2) - (n - 1) * colSums(m1_beta^2))/ 2

}

update_kappa_z_vb_ <- function(Y, X, Z, kappa, mu_alpha_vb, m1_beta, m2_alpha,
                               m2_beta, mat_x_m1, mat_z_mu, sig2_inv_vb,
                               zeta2_inv_vb, intercept = FALSE) {
  n <- nrow(Y)

  kappa_vb <- kappa + (colSums(Y^2) - 2 * colSums(Y * (mat_x_m1 + mat_z_mu))  +
                         (n - 1 + sig2_inv_vb) * colSums(m2_beta) +
                         colSums(mat_x_m1^2) - (n - 1) * colSums(m1_beta^2) +
                         (n - 1) * colSums(m2_alpha) +
                         crossprod(m2_alpha, zeta2_inv_vb) +
                         colSums(mat_z_mu^2) - (n - 1) * colSums(mu_alpha_vb^2) +
                         2 * colSums(mat_x_m1 * mat_z_mu))/ 2

  if (intercept)
    kappa_vb <- kappa_vb + (m2_alpha[1, ] - (mu_alpha_vb[1, ])^2) / 2

  kappa_vb
}

update_log_tau_vb_ <- function(eta_vb, kappa_vb) digamma(eta_vb) - log(kappa_vb)


#################
## W's updates ##
#################

update_W_info_ <- function(gam_vb, mat_v_mu) {

  gam_vb * (inv_mills_ratio_(1, mat_v_mu) - inv_mills_ratio_(0, mat_v_mu)) +
    mat_v_mu + inv_mills_ratio_(0, mat_v_mu)

}

update_W_probit_ <- function(Y, mat_z_mu, mat_x_m1) mat_z_mu + mat_x_m1 + inv_mills_ratio_(Y, mat_z_mu + mat_x_m1)


####################
## zeta's updates ##
####################

update_phi_z_vb_ <- function(phi, d) phi + d / 2

update_xi_z_vb_ <- function(xi, tau_vb, m2_alpha) xi + m2_alpha %*% tau_vb / 2

update_xi_bin_vb_ <- function(xi, m2_alpha) xi + rowSums(m2_alpha) / 2

update_log_zeta2_inv_vb_ <- function(phi_vb, xi_vb) digamma(phi_vb) - log(xi_vb)
