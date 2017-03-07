## --- Y --- ##

update_nu_bin_vb_ <- function(nu, m2_beta) nu + sum(m2_beta) / 2

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

## --------- ##


## --- X --- ##

update_m1_beta_ <- function(gam_vb, mu_beta_vb) gam_vb * mu_beta_vb

update_m2_beta_ <- function(gam_vb, mu_beta_vb, sig2_beta_vb) sweep(mu_beta_vb ^ 2, 2, sig2_beta_vb, `+`) * gam_vb

update_sig2_beta_vb_ <- function(n, sig2_inv_vb, tau_vb) 1 / ((n - 1 + sig2_inv_vb) * tau_vb)

update_mat_x_m1_ <- function(X, m1_beta) X %*% m1_beta


##
update_lambda_vb_ <- function(lambda, sum_gam) lambda + sum_gam / 2

update_nu_vb_ <- function(nu, m2_beta, tau_vb) as.numeric(nu + crossprod(tau_vb, colSums(m2_beta)) / 2)

update_log_sig2_inv_vb_ <- function(lambda_vb, nu_vb) digamma(lambda_vb) - log(nu_vb)
##

## --------- ##


## --- V --- ##

update_mu_c0_vb_ <- function(W, mat_v_mu, m0, s02, sig2_c0_vb) sig2_c0_vb * (rowSums(W - mat_v_mu) + m0 / s02)

update_sig2_c0_vb_ <- function(d, s02) 1 / (d + (1/s02))

update_sig2_c_vb_ <- function(p, s2) 1 / (p - 1 + (1/s2))

update_mat_v_mu_ <- function(V, mu_c0_vb, mu_c_vb) sweep(V %*% mu_c_vb, 1, mu_c0_vb, `+`)

## --------- ##


## --- W --- ##

update_W_info_ <- function(gam_vb, mat_v_mu) {

  gam_vb * (inv_mills_ratio_(1, mat_v_mu) - inv_mills_ratio_(0, mat_v_mu)) +
    mat_v_mu + inv_mills_ratio_(0, mat_v_mu)

}

update_W_probit_ <- function(Y, mat_z_mu, mat_x_m1) mat_z_mu + mat_x_m1 + inv_mills_ratio_(Y, mat_z_mu + mat_x_m1)


## --------- ##



## --- Z --- ##

update_m2_alpha_ <- function(mu_alpha_vb, sig2_alpha_vb) sig2_alpha_vb + mu_alpha_vb ^ 2

update_mat_z_mu_ <- function(Z, mu_alpha_vb) Z %*% mu_alpha_vb

update_phi_z_vb_ <- function(phi, d) phi + d / 2

update_xi_z_vb_ <- function(xi, tau_vb, m2_alpha) xi + m2_alpha %*% tau_vb / 2

update_xi_bin_vb_ <- function(xi, m2_alpha) xi + rowSums(m2_alpha) / 2

## --------- ##


## --- auxiliary parameter logistic --- ##

update_psi_logit_vb_ <- function(chi_vb) {

  sig <- function(chi) {
    1 / (1 + exp(-chi))
  }

  (sig(chi_vb) - 1 / 2) / (2 * chi_vb)

}

## --------- ##




# update_mu_beta_ <- function(X, Y, mat_x_m1, sig2_beta_vb, tau_vb, j, k = NULL) {
#
#   #if (is.null(k)) {
#
#   sig2_beta_vb * (tau_vb * crossprod(Y - mat_x_m1, X[, j]))
#
#   # } else {  # not needed if batch = T
#   #
#   #   sig2_beta_vb[k] * tau_vb[k] * crossprod(Y[, k] - mat_x_m1[, k], X[, j])
#   #
#   # }
#
# }

# update_info_gam_vb_ <- function(log_sig2_inv_vb, log_tau_vb, mat_v_mu,
#                                 mu_beta_vb, sig2_beta_vb, j, k = NULL) {
#
#   # if (is.null(k)) {
#
#   log_part_gam_vb <- pnorm(mat_v_mu[j, ], log.p = TRUE) + log(sig2_beta_vb) / 2 +
#     mu_beta_vb[j, ] ^ 2 / (2 * sig2_beta_vb)
#
#   log_part2_gam_vb <- pnorm(mat_v_mu[j, ], lower.tail = FALSE, log.p = TRUE) -
#     log_tau_vb / 2 - log_sig2_inv_vb / 2
#
#   exp(log_part_gam_vb - log_sum_exp_mat_(list(log_part_gam_vb, log_part2_gam_vb)))
#
#   # } else { # not needed if batch = T
#   #
#   #   log_part_gam_vb <-  pnorm(mat_v_mu[j, k], log.p = TRUE) +
#   #     log(sig2_beta_vb[k]) / 2 + mu_beta_vb[j, k] ^ 2 / (2 * sig2_beta_vb[k])
#   #
#   #   log_part2_gam_vb <- pnorm(mat_v_mu[j, k], lower.tail = FALSE, log.p = TRUE) -
#   #     log_tau_vb[k] / 2 - log_sig2_inv_vb / 2
#   #
#   #   exp(log_part_gam_vb - log_sum_exp_(c(log_part_gam_vb, log_part2_gam_vb)))
#   # }
#
# }


# update_mu_c_vb_ <- function(V, W, mat_v_mu, sig2_c_vb, l, k = NULL) {
#
#   # if (is.null(k)) {
#
#   sig2_c_vb * crossprod(W - mat_v_mu, V[, l])
#
#   # } else {  # not needed if batch = T
#   #
#   #   sig2_c_vb * crossprod(W[, k] - mat_v_mu[, k], V[, l])
#   #
#   # }
# }

# update_m1_beta_ <- function(gam_vb, mu_beta_vb, j = NULL, k = NULL) {
#
#   nj <- is.null(j)
#   nk <- is.null(k)
#
#   if (nj & nk) {
#     gam_vb * mu_beta_vb
#   } else if (nk) {
#     gam_vb[j, ] * mu_beta_vb[j, ]
#   } else {
#     gam_vb[j, k] * mu_beta_vb[j, k] # not needed if batch = T
#   }
#
# }
