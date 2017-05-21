# This file is part of the `locus` R package:
#     https://github.com/hruffieux/locus
#
# Internal functions gathering the ELBO terms common to core algorithms.
#

############################################
## E log p(alpha | rest) - E log q(alpha) ##
############################################

e_alpha_ <- function(m2_alpha, log_tau_vb, log_zeta2_inv_vb, sig2_alpha_vb,
                     tau_vb, zeta2_inv_vb) {

  1 / 2 * sum(sweep( sweep( sweep( sweep(m2_alpha, 2, tau_vb, `*`), 1,
                                   -zeta2_inv_vb, `*`), 2, log_tau_vb, `+`), 1,
                     log_zeta2_inv_vb , `+`) + log(sig2_alpha_vb) + 1)


}


e_alpha_logit_ <- function(m2_alpha, log_zeta2_inv_vb, sig2_alpha_vb, zeta2_inv_vb) {

  1 / 2 * sum( sweep(-sweep(m2_alpha, 1, zeta2_inv_vb, `*`), 1,
                     log_zeta2_inv_vb, `+`) + log(sig2_alpha_vb) + 1)

}


e_alpha_probit_ <- function(m2_alpha, log_zeta2_inv_vb, sig2_alpha_vb, zeta2_inv_vb) {

  1 / 2 * sum(sweep(-sweep(m2_alpha, 1, zeta2_inv_vb, `*`), 1,
                    log_zeta2_inv_vb + log(sig2_alpha_vb), `+`) + 1)

}


########################################################
## E log p(beta, gamma | rest) - E log q(beta, gamma) ##
########################################################


e_beta_gamma_ <- function(gam_vb, log_om_vb, log_1_min_om_vb, log_sig2_inv_vb,
                          log_tau_vb, m2_beta, sig2_beta_vb, sig2_inv_vb, tau_vb) {

  eps <- .Machine$double.eps # to control the argument of the log when gamma is very small
  sum(log_sig2_inv_vb * gam_vb / 2 + sweep(gam_vb, 2, log_tau_vb, `*`) / 2 -
        sweep(m2_beta, 2, tau_vb, `*`) * sig2_inv_vb / 2 +
        sweep(gam_vb, 1, log_om_vb, `*`) +
        sweep(1 - gam_vb, 1, log_1_min_om_vb, `*`) +
        1 / 2 * sweep(gam_vb, 2, log(sig2_beta_vb) + 1, `*`) -
        gam_vb * log(gam_vb + eps) - (1 - gam_vb) * log(1 - gam_vb + eps))

}


e_g_beta_gamma_ <- function(gam_vb, g_sizes, log_om_vb, log_1_min_om_vb,
                            log_sig2_inv_vb, log_tau_vb, list_m1_btb,
                            list_sig2_beta_star, sig2_inv_vb, tau_vb, vec_log_det) {

  eps <- .Machine$double.eps # to control the argument of the log when gamma is very small

  G <- length(list_m1_btb)

  sum(unlist(lapply(1:G, function(g) {

    sum(g_sizes[g] / 2 * gam_vb[g, ] * (log_sig2_inv_vb + log_tau_vb) -
          list_m1_btb[[g]] * tau_vb * sig2_inv_vb / 2 +
          gam_vb[g, ] * log_om_vb[g] + (1 - gam_vb[g, ]) * log_1_min_om_vb[g] +
          1 / 2 * gam_vb[g, ] * (vec_log_det[g] - g_sizes[g] * (log(tau_vb) - 1)) -
          gam_vb[g, ] * log(gam_vb[g, ] + eps) - (1 - gam_vb[g, ]) * log(1 - gam_vb[g, ] + eps))
  })))

}


e_beta_gamma_bin_ <- function(gam_vb, log_om_vb, log_1_min_om_vb, log_sig2_inv_vb,
                                m2_beta, sig2_beta_vb, sig2_inv_vb) {

  eps <- .Machine$double.eps # to control the argument of the log when gamma is very small
  sum(gam_vb * log_sig2_inv_vb / 2 - m2_beta * sig2_inv_vb / 2 +
           sweep(gam_vb, 1, log_om_vb, `*`) +
           sweep(1 - gam_vb, 1, log_1_min_om_vb, `*`) +
           gam_vb * (log(sig2_beta_vb) + 1) / 2 -
           gam_vb * log(gam_vb + eps) - (1 - gam_vb) * log(1 - gam_vb + eps))

}


e_beta_gamma_info_ <- function(V, gam_vb, log_sig2_inv_vb, log_tau_vb, mat_v_mu,
                               m2_beta, sig2_beta_vb, sig2_c0_vb, sig2_c_vb,
                               sig2_inv_vb, tau_vb) {

  eps <- .Machine$double.eps^0.75 # to control the argument of the log when gamma is very small

  sum(log_sig2_inv_vb * gam_vb / 2 +
        sweep(gam_vb, 2, log_tau_vb, `*`) / 2 -
        sweep(m2_beta, 2, tau_vb, `*`) * sig2_inv_vb / 2 +
        gam_vb * pnorm(mat_v_mu, log.p = TRUE) +
        sweep((1 - gam_vb) * pnorm(mat_v_mu, lower.tail = FALSE, log.p = TRUE), 1, sig2_c_vb * rowSums(V^2) / 2, `-`) -
        sig2_c0_vb / 2 + 1 / 2 * sweep(gam_vb, 2, log(sig2_beta_vb) + 1, `*`) -
        gam_vb * log(gam_vb + eps) - (1 - gam_vb) * log(1 - gam_vb + eps))

}


e_beta_gamma_info_bin_ <- function(V, gam_vb, log_sig2_inv_vb, mat_v_mu, m2_beta,
                                   sig2_beta_vb, sig2_c0_vb, sig2_c_vb, sig2_inv_vb) {

  eps <- .Machine$double.eps^0.75 # to control the argument of the log when gamma is very small

  sum(gam_vb * log_sig2_inv_vb / 2 - m2_beta * sig2_inv_vb / 2 +
        gam_vb * pnorm(mat_v_mu, log.p = TRUE) +
        sweep((1 - gam_vb) * pnorm(mat_v_mu, lower.tail = FALSE, log.p = TRUE), 1, sig2_c_vb * rowSums(V^2) / 2, `-`) -
        sig2_c0_vb / 2 + gam_vb * (log(sig2_beta_vb) + 1) / 2 -
        gam_vb * log(gam_vb + eps) - (1 - gam_vb) * log(1 - gam_vb + eps))

}


e_beta_gamma_struct_ <- function(gam_vb, log_sig2_inv_vb, log_tau_vb, mu_theta_vb,
                                 m2_beta, sig2_beta_vb, list_sig2_theta_vb,
                                 sig2_inv_vb, tau_vb) {

  eps <- .Machine$double.eps^0.75 # to control the argument of the log when gamma is very small

  diag_sig2_theta_vb <- unlist(lapply(list_sig2_theta_vb, diag))

  sum(log_sig2_inv_vb * gam_vb / 2 +
        sweep(gam_vb, 2, log_tau_vb, `*`) / 2 -
        sweep(m2_beta, 2, tau_vb, `*`) * sig2_inv_vb / 2 +
        sweep(gam_vb, 1, pnorm(mu_theta_vb, log.p = TRUE), `*`) +
        sweep(sweep((1 - gam_vb), 1, pnorm(mu_theta_vb, lower.tail = FALSE, log.p = TRUE), `*`), 1, diag_sig2_theta_vb / 2, `-`) +
        1 / 2 * sweep(gam_vb, 2, log(sig2_beta_vb) + 1, `*`) -
        gam_vb * log(gam_vb + eps) - (1 - gam_vb) * log(1 - gam_vb + eps))

}


######################################
## E log p(c0 | rest) - E log q(c0) ##
######################################

e_c0_ <- function(m0, mu_c0_vb, s02, sig2_c0_vb) {

  sum(log(sig2_c0_vb) + 1 - log(s02) -
        (mu_c0_vb^2 + sig2_c0_vb - 2*mu_c0_vb * m0 + m0^2) / s02) / 2

}


######################################
## E log p(c | rest) - E log q(c) ##
######################################

e_c_ <- function(mu_c_vb, s2, sig2_c_vb) {

  sum(log(sig2_c_vb) + 1 - log(s2) - (mu_c_vb^2 + sig2_c_vb) / s2) / 2

}


############################################
## E log p(omega | rest) - E log q(omega) ##
############################################

e_omega_ <- function(a, a_vb, b, b_vb, log_om_vb, log_1_min_om_vb) {

  sum((a - a_vb) * log_om_vb + (b - b_vb) * log_1_min_om_vb - lbeta(a, b) + lbeta(a_vb, b_vb))

}


##################################################
## E log p(sig2_inv | rest) - E log q(sig2_inv) ##
##################################################

e_sig2_inv_ <- function(lambda, lambda_vb, log_sig2_inv_vb, nu, nu_vb, sig2_inv_vb) {

  (lambda - lambda_vb) * log_sig2_inv_vb - (nu - nu_vb) * sig2_inv_vb +
    lambda * log(nu) - lambda_vb * log(nu_vb) - lgamma(lambda) + lgamma(lambda_vb)

}


########################################
## E log p(tau | rest) - E log q(tau) ##
########################################

e_tau_ <- function(eta, eta_vb, kappa, kappa_vb, log_tau_vb, tau_vb) {

  sum((eta - eta_vb) * log_tau_vb - (kappa - kappa_vb) * tau_vb +
        eta * log(kappa) - eta_vb * log(kappa_vb) - lgamma(eta) + lgamma(eta_vb))

}


############################################
## E log p(theta | rest) - E log q(theta) ##
############################################

# S0_inv is assumed to be block-diagonal
e_theta_ <- function(m0, mu_theta_vb, list_S0_inv, list_sig2_theta_vb, vec_fac_st) {

  bl_ids <- unique(vec_fac_st)
  n_bl <- length(list_S0_inv)

  arg <- unlist(lapply(1:n_bl, function(bl) {

    mu_theta_vb_bl <- mu_theta_vb[vec_fac_st == bl_ids[bl]]
    m0_bl <- m0[vec_fac_st == bl_ids[bl]]
    S0_inv_bl <- list_S0_inv[[bl]]
    sig2_theta_vb_bl <- list_sig2_theta_vb[[bl]]

    (log(det(S0_inv_bl)) + log(det(sig2_theta_vb_bl)) -
       crossprod((mu_theta_vb_bl - m0_bl),
                 S0_inv_bl %*% (mu_theta_vb_bl - m0_bl)) -
       sum(S0_inv_bl * sig2_theta_vb_bl) + ncol(S0_inv_bl)) / 2 # trace of a product
  }))

  sum(arg)

}


#######################
## E log p(y | rest) ##
#######################

e_y_ <- function(n, kappa, kappa_vb, log_tau_vb, m2_beta, sig2_inv_vb, tau_vb,
                 m2_alpha = NULL, zeta2_inv_vb = NULL) {

  arg <- -n / 2 * log(2 * pi) + n / 2 * log_tau_vb - tau_vb *
    (kappa_vb - colSums(m2_beta) * sig2_inv_vb / 2 - kappa)

  if (!is.null(m2_alpha))
    arg <- arg + tau_vb * crossprod(m2_alpha, zeta2_inv_vb) / 2

  sum(arg)

}

e_g_y_ <- function(n, kappa, kappa_vb, list_m1_btb, log_tau_vb, sig2_inv_vb, tau_vb) {

  sum(-n / 2 * log(2 * pi) + n / 2 * log_tau_vb - tau_vb *
    (kappa_vb - Reduce(`+`, list_m1_btb) * sig2_inv_vb / 2 - kappa))

}


e_y_logit_ <- function(X, Y, Z, chi_vb, m1_beta, m2_alpha, m2_beta, mat_x_m1,
                       mat_z_mu, mu_alpha_vb, psi_vb) {

  sum(log_sigmoid_(chi_vb) + Y * (mat_x_m1 + mat_z_mu)  -  chi_vb / 2 -
        psi_vb * (X^2 %*% m2_beta + mat_x_m1^2 - X^2 %*% m1_beta^2 +
                    Z^2 %*% m2_alpha + mat_z_mu^2 - Z^2 %*% mu_alpha_vb^2 +
                    2 * mat_x_m1 * mat_z_mu - chi_vb^2))

}


e_y_probit_ <- function(X, Y, Z, m1_beta, m2_beta, mat_x_m1, mat_z_mu,
                        sig2_alpha_vb, sweep = TRUE) {

  U <- mat_x_m1 + mat_z_mu

  arg <- Y * pnorm(U, log.p = TRUE) +
    (1 - Y) * pnorm(U, lower.tail = FALSE, log.p = TRUE) -
    X^2 %*% (m2_beta - m1_beta^2) / 2

  if (sweep)
    arg <- sweep(arg, 1, Z^2 %*% sig2_alpha_vb / 2, `-`)
  else
    arg <- arg - Z^2 %*% sig2_alpha_vb / 2

  sum(arg)

}


####################################################
## E log p(zeta2_inv | rest) - E log q(zeta2_inv) ##
####################################################

e_zeta2_inv_ <- function(log_zeta2_inv_vb, phi, phi_vb, xi, xi_vb, zeta2_inv_vb) {

  sum((phi - phi_vb) * log_zeta2_inv_vb - (xi - xi_vb) * zeta2_inv_vb +
        phi * log(xi) - phi_vb * log(xi_vb) - lgamma(phi) + lgamma(phi_vb))

}
