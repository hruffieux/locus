#' Gather model hyperparameters provided by the user.
#'
#' This function must be used to provide hyperparameter values for the model
#' used in \code{\link{locus}}.
#'
#' The \code{\link{locus}} function can also be used with default
#' hyperparameter choices (without using \code{\link{set_hyper}}) by
#' setting its argument \code{list_hyper} to \code{NULL}.
#'
#' @param d Number of responses.
#' @param p Number of candidate predictors.
#' @param lambda Vector of length 1 for \code{family = "gaussian"} and of length
#'   1 or d for \code{family = "binomial-logit"} and
#'   \code{family = "binomial-probit"} providing the values of hyperparameter
#'   \eqn{\lambda} for the prior distribution of \eqn{\sigma^{-2}}. If of length
#'   1 for \code{family = "binomial-logit"} or \code{family = "binomial-probit"},
#'   the provided value is repeated p times. \eqn{\sigma^2} represents the
#'   typical size of nonzero effects.
#' @param nu Vector of length 1 for \code{family = "gaussian"} and of length
#'   1 or d for \code{family = "binomial-logit"} and
#'   \code{family = "binomial-probit"} providing the values of hyperparameter
#'   \eqn{\nu} for the prior distribution of \eqn{\sigma^{-2}}. If of length 1
#'   for \code{family = "binomial-logit"} or \code{family = "binomial-probit"},
#'   the provided value is repeated p times. \eqn{\sigma^2} represents the
#'   typical size of nonzero effects.
#' @param a Vector of length 1 or p providing the values of hyperparameter
#'   \eqn{a} for the prior distributions for the proportion of responses
#'   associated with each candidate predictor, \eqn{\omega} (vector of size p).
#'   If of length 1, the provided value is repeated p times.
#' @param b Vector of length 1 or p providing the values of hyperparameter
#'   \eqn{b} for the prior distributions for the proportion of responses
#'   associated with each candidate predictor, \eqn{\omega} (vector of size p).
#'   If of length 1, the provided value is repeated p times.
#' @param eta Vector of length 1 or d for \code{family = "gaussian"},
#'   providing the values of hyperparameter \eqn{\eta} for the prior
#'   distributions of the response residual precisions, \eqn{\tau}
#'   (vector of size d). If of length 1, the provided value is repeated d times.
#'   Must be \code{NULL} for \code{family = "binomial-logit"} and
#'   \code{family = "binomial-probit"}.
#' @param kappa Vector of length 1 or d for \code{family = "gaussian"},
#'   providing the values of hyperparameter \eqn{\kappa} for the prior
#'   distributions of the response residual precisions, \eqn{\tau}
#'   (vector of size d). If of length 1, the provided value is repeated d times.
#'   Must be \code{NULL} for \code{family = "binomial-logit"} and
#'   \code{family = "binomial-probit"}.
#' @param family Response type. Must be either "\code{gaussian}" for linear
#'   regression, "\code{binomial-logit}" for logistic regression or
#'   "\code{binomial-probit}" for probit regression.
#' @param q Number of covariates. Default is \code{NULL}, for \code{Z}
#'   \code{NULL}.
#' @param phi For \code{family = "gaussian"}, vector of length 1 or q providing
#'   the values of hyperparameter \eqn{\phi} for the prior distributions for the
#'   sizes of the nonzero covariate effects, \eqn{\zeta}. If of length 1, the
#'   provided value is repeated q times. For \code{family = "binomial-logit"} or
#'   \code{family = "binomial-probit"}, matrix of dimension q x d as the values
#'   are also specific to the responses. Default is \code{NULL}, for \code{Z}
#'   \code{NULL}.
#' @param xi For \code{family = "gaussian"}, vector of length 1 or q providing
#'   the values of hyperparameter \eqn{\xi} for the prior distributions for the
#'   sizes of the nonzero covariate effects, \eqn{\zeta}. If of length 1, the
#'   provided value is repeated q times. For \code{family = "binomial-logit"} or
#'   \code{family = "binomial-probit"}, matrix of dimension q x d as the values
#'   are also specific to the responses. Default is \code{NULL}, for \code{Z}
#'   \code{NULL}.
#'
#' @return An object of class "\code{hyper}" preparing user hyperparameter in a
#'   form that can be passed to the \code{\link{locus}} function.
#'
#' @examples
#' user_seed <- 123
#' n <- 200; p <- 250; p0 <- 50; d <- 25; d0 <- 15
#' list_X <- generate_snps(n = n, p = p, user_seed = user_seed)
#' list_Y <- generate_phenos(n = n, d = d, var_err = 1, user_seed = user_seed)
#'
#' # Gaussian outcomes
#' dat_g <- generate_dependence(list_snps = list_X, list_phenos = list_Y,
#'                              ind_d0 = sample(1:d, d0),
#'                              ind_p0 = sample(1:p, p0), vec_prob_sh = 0.1,
#'                              family = "gaussian", max_tot_pve = 0.9,
#'                              user_seed = user_seed)
#'
#' # a and b chosen so that each candidate predictor has a prior probability to
#' # be included in the model of 1/4.
#' list_hyper_g <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
#'                           eta = 1, kappa = apply(dat_g$phenos, 2, var),
#'                           family = "gaussian")
#'
#' vb_g <- locus(Y = dat_g$phenos, X = dat_g$snps, p0_av = p0,
#'               family = "gaussian", list_hyper = list_hyper_g,
#'               user_seed = user_seed)
#'
#' # Gaussian outcomes with covariates
#' q <- 4
#' Z <- matrix(rnorm(n * q), nrow = n)
#'
#' phi <- xi <- rep(1, q)
#'
#' list_hyper_g_z <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
#'                             eta = 1, kappa = apply(dat_g$phenos, 2, var),
#'                             family = "gaussian", q = q, phi = phi, xi = xi)
#'
#' vb_g_z <- locus(Y = dat_g$phenos, X = dat_g$snps, p0_av = p0, Z = Z,
#'                 family = "gaussian", list_hyper = list_hyper_g_z,
#'                 user_seed = user_seed)
#'
#' # Binary outcomes
#' dat_b <- generate_dependence(list_snps = list_X, list_phenos = list_Y,
#'                              ind_d0 = sample(1:d, d0),
#'                              ind_p0 = sample(1:p, p0), vec_prob_sh = 0.1,
#'                              family = "binomial", max_tot_pve = 0.9,
#'                              user_seed = user_seed)
#'
#' # a and b chosen so that each candidate predictor has a prior probability to
#' # be included in the model of 1/4.
#' list_hyper_logit <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
#'                               eta = NULL, kappa = NULL,
#'                               family = "binomial-logit")
#'
#' p0_av <- floor(4*p/5) # overestimating the prior number of active covariates
#'                       # often leads to better inference
#'
#' vb_logit <- locus(Y = dat_b$phenos, X = dat_b$snps, p0_av = p0_av,
#'                   family = "binomial-logit", list_hyper = list_hyper_logit,
#'                   user_seed = user_seed)
#'
#' list_hyper_probit <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
#'                                eta = NULL, kappa = NULL,
#'                                family = "binomial-probit")
#'
#' vb_probit <- locus(Y = dat_b$phenos, X = dat_b$snps, p0_av = p0_av,
#'                   family = "binomial-probit", list_hyper = list_hyper_probit,
#'                   user_seed = user_seed)
#'
#' # Binary outcomes with covariates
#' phi <- xi <- matrix(1, nrow = q, ncol = d)
#'
#' list_hyper_logit_z <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
#'                                 eta = NULL, kappa = NULL,
#'                                 family = "binomial-logit", q = q, phi = phi,
#'                                 xi = xi)
#'
#' vb_logit_z <- locus(Y = dat_b$phenos, X = dat_b$snps, p0_av = p0_av, Z = Z,
#'                     family = "binomial-logit",
#'                     list_hyper = list_hyper_logit_z, user_seed = user_seed)
#'
#' list_hyper_probit_z <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
#'                                  eta = NULL, kappa = NULL,
#'                                  family = "binomial-probit", q = q, phi = phi,
#'                                  xi = xi)
#'
#' vb_probit_z <- locus(Y = dat_b$phenos, X = dat_b$snps, p0_av = p0_av, Z = Z,
#'                      family = "binomial-probit",
#'                      list_hyper = list_hyper_probit_z, user_seed = user_seed)
#'
#' @seealso  \code{\link{set_init}}, \code{\link{locus}}
#'
#' @export
#'
set_hyper <- function(d, p, lambda, nu, a, b, eta, kappa, family = "gaussian",
                      q = NULL, phi = NULL, xi = NULL) {

  check_structure_(d, "vector", "numeric", 1)
  check_natural_(d)

  check_structure_(p, "vector", "numeric", 1)
  check_natural_(p)

  check_structure_(q, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(q)) check_natural_(q)

  stopifnot(family %in% c("gaussian", "binomial-logit", "binomial-probit"))

  check_structure_(a, "vector", "double", c(1, p))
  check_positive_(a)
  if (length(a) == 1) a <- rep(a, p)

  check_structure_(b, "vector", "double", c(1, p))
  check_positive_(b)
  if (length(b) == 1) b <- rep(b, p)

  if (family == "gaussian") {

    check_structure_(lambda, "vector", "double", 1)
    check_positive_(lambda)

    check_structure_(nu, "vector", "double", 1)
    check_positive_(nu)

    check_structure_(eta, "vector", "double", c(1, d))
    check_positive_(eta)
    if (length(eta) == 1) eta <- rep(eta, d)

    check_structure_(kappa, "vector", "double", c(1, d))
    check_positive_(kappa)
    if (length(kappa) == 1) kappa <- rep(kappa, d)

  } else {

    check_structure_(lambda, "vector", "double", c(1, d))
    check_positive_(lambda)
    if (length(lambda) == 1) lambda <- rep(lambda, d)

    check_structure_(nu, "vector", "double", c(1, d))
    check_positive_(nu)
    if (length(nu) == 1) nu <- rep(nu, d)

    if (!is.null(eta) | !is.null(kappa))
      stop("Both eta and kappa must be NULL for logistic regression.")
  }

  if (!is.null(q)) {

    if (family == "gaussian") {

      check_structure_(phi, "vector", "double", c(1, q))
      check_positive_(phi)
      if (length(phi) == 1) phi <- rep(phi, q)

      check_structure_(xi, "vector", "double", c(1, q))
      check_positive_(xi)
      if (length(xi) == 1) xi <- rep(xi, q)

    } else {

      check_structure_(phi, "matrix", "double", c(q, d))
      check_positive_(phi)

      check_structure_(xi, "matrix", "double", c(q, d))
      check_positive_(xi)

    }


  } else if (!is.null(phi) | !is.null(xi)) {
    stop("Provided q = NULL, not consitent with phi or xi being non-null.")
  }

  d_hyper <- d
  p_hyper <- p
  q_hyper <- q

  family_hyper <- family

  list_hyper <- create_named_list_(d_hyper, p_hyper, q_hyper, family_hyper, eta,
                                   kappa, lambda, nu, a, b, phi, xi)

  class(list_hyper) <- "hyper"

  list_hyper

}


auto_set_hyper_ <- function(Y, p, p_star, q, family) {

  d <- ncol(Y)

  lambda <- 1e-2
  nu <- 1

  if (family == "gaussian") {
    # hyperparameter set using the data Y
    eta <- 1 / median(apply(Y, 2, var)) #median to be consistent when doing permutations
    if (!is.finite(eta)) eta <- 1e3
    eta <- rep(eta, d)
    kappa <- rep(1, d)
  } else {
    eta <- kappa <- NULL

    lambda <- rep(lambda, d)
    nu <- rep(nu, d)
  }

  # if p_star is of length 1, p_star is the prior average number of active
  # predictors else (p_star is of length p), p_star / p is the vector containg
  # the prior probabilities that each predictor is active and the sum of its
  # entries is the corresponding prior average number of active predictors
  if (length(p_star) == 1) p0 <- p_star
  else p0 <- sum(p_star / p)

  a <- rep(1, p)
  b <- d * (p - p_star) / p_star
  if (length(b) == 1) b <- rep(b, p)

  # hyperparameters of beta distributions
  check_positive_(a)
  check_positive_(b)

  if (!is.null(q)) {

    if (family == "gaussian") {

      phi <- xi <- rep(1, q)

    } else {

      phi <- xi <- matrix(1, nrow = q, ncol = d)

    }

  } else {

    phi <- NULL
    xi <- NULL

  }

  d_hyper <- d
  p_hyper <- p
  q_hyper <- q

  family_hyper <- family

  list_hyper <- create_named_list_(d_hyper, p_hyper, q_hyper, family_hyper, eta,
                                   kappa, lambda, nu, a, b, phi, xi)

  class(list_hyper) <- "out_hyper"

  list_hyper

}

#' Gather initial variational parameters provided by the user.
#'
#' This function must be used to provide initial values for the variational
#' parameters used in \code{\link{locus}}.
#'
#' The \code{\link{locus}} function can also be used with default initial
#' parameter choices (without using \code{\link{set_init}}) by setting
#' its argument \code{list_init} to \code{NULL}.
#'
#' @param d Number of responses.
#' @param p Number of candidate predictors.
#' @param gam_vb Matrix of size p x d with initial values for the variational
#'   parameter yielding posterior probabilities of inclusion.
#' @param mu_beta_vb Matrix of size p x d with initial values for the
#'   variationa parameter yielding regression coefficient estimates for
#'   predictor-response pairs included in the model.
#' @param sig2_beta_vb Vector of size d, for \code{family = "gaussian"} and
#'   \code{family = "binomial-probit"}, or matrix of size p x d, for
#'   \code{family = "binomial-logit"}, with initial values for the variational
#'   parameter yielding estimates of effect variances for predictor-response
#'   pairs included in the model. For \code{family = "gaussian"} and
#'   \code{family = "binomial-probit"}, these values are the same for all
#'   predictors (as a result of the predictor variables being standardized
#'   before the variational algorithm).
#' @param tau_vb  Vector of size d with initial values, for
#'   \code{family = "gaussian"}, for the variational parameter yielding
#'   estimates for the response residual precisions. Must be \code{NULL} for
#'   \code{family = "binomial-logit"} and \code{family = "binomial-probit"}.
#' @param family Response type. Must be either "\code{gaussian}" for linear
#'   regression, "\code{binomial-logit}" for logistic regression or
#'   \code{family = "binomial-probit"} for probit regression.
#' @param n Number of observations. Used only when
#'   \code{family = "binomial-logit"}.
#' @param q Number of covariates. Default is \code{NULL}, for \code{Z}
#'   \code{NULL}.
#' @param mu_alpha_vb Matrix of size p x q with initial values for the
#'   variational parameter yielding regression coefficient estimates for
#'   covariate-response pairs. Default is \code{NULL}, for \code{Z}
#'   \code{NULL}.
#' @param sig2_alpha_vb Matrix of size p x q with initial values for the
#'   variational parameter yielding estimates of effect variances for
#'   covariate-response pairs. Default is \code{NULL}, for \code{Z}
#'   \code{NULL}.
#'
#' @return An object of class "\code{init}" preparing user initial values for
#'   the variational parameters in a form that can be passed to the
#'   \code{\link{locus}} function.
#'
#' @examples
#' user_seed <- 123; set.seed(user_seed)
#' n <- 200; p <- 250; p0 <- 50; d <- 25; d0 <- 15
#' list_X <- generate_snps(n = n, p = p)
#' list_Y <- generate_phenos(n = n, d = d, var_err = 1)
#'
#' # Gaussian outcomes
#' dat_g <- generate_dependence(list_snps = list_X, list_phenos = list_Y,
#'                              ind_d0 = sample(1:d, d0),
#'                              ind_p0 = sample(1:p, p0),
#'                              vec_prob_sh = 0.1, family = "gaussian",
#'                              max_tot_pve = 0.9)
#'
#' # gam_vb chosen so that each candidate predictor has a prior probability to
#' # be included in the model of 1/4.
#' gam_vb <- matrix(rbeta(p * d, shape1 = 1, shape2 = 4*d-1), nrow = p)
#' mu_beta_vb <- matrix(rnorm(p * d), nrow = p)
#' tau_vb <- 1 / apply(dat_g$phenos, 2, var)
#' sig2_beta_vb <- 1 / rgamma(d, shape = 2, rate = 1 / tau_vb)
#'
#' list_init_g <- set_init(d, p, gam_vb, mu_beta_vb, sig2_beta_vb, tau_vb,
#'                         family = "gaussian")
#'
#' vb_g <- locus(Y = dat_g$phenos, X = dat_g$snps, p0_av = p0,
#'               family = "gaussian", list_init = list_init_g)
#'
#' # Gaussian outcomes with covariates
#' q <- 4
#' Z <- matrix(rnorm(n * q), nrow = n)
#'
#' mu_alpha_vb <- matrix(rnorm(q * d), nrow = q)
#' sig2_alpha_vb <- 1 / matrix(rgamma(q * d, shape = 2, rate = 1), nrow = q)
#'
#' list_init_g_z <- set_init(d, p, gam_vb, mu_beta_vb, sig2_beta_vb, tau_vb,
#'                           family = "gaussian", q = q,
#'                           mu_alpha_vb = mu_alpha_vb,
#'                           sig2_alpha_vb = sig2_alpha_vb)
#'
#' vb_g_z <- locus(Y = dat_g$phenos, X = dat_g$snps, p0_av = p0, Z = Z,
#'                 family = "gaussian", list_init = list_init_g_z)
#'
#' # Binary outcomes
#' dat_b <- generate_dependence(list_snps = list_X, list_phenos = list_Y,
#'                              ind_d0 = sample(1:d, d0),
#'                              ind_p0 = sample(1:p, p0),
#'                              vec_prob_sh = 0.1, family = "binomial",
#'                              max_tot_pve = 0.9)
#'
#' p0_av <- floor(4*p/5) # overestimating the prior number of active covariates
#'                       # often leads to better inference
#'
#' # gam_vb chosen so that each candidate predictor has a prior probability to
#' # be included in the model of 1/4.
#'
#' sig2_beta_vb_logit <- 1 / t(replicate(p, rgamma(d, shape = 2, rate = 1)))
#'
#' list_init_logit <- set_init(d, p, gam_vb, mu_beta_vb, sig2_beta_vb_logit,
#'                             tau_vb = NULL, family = "binomial-logit", n = n)
#'
#' vb_logit <- locus(Y = dat_b$phenos, X = dat_b$snps, p0_av = p0_av,
#'                   family = "binomial-logit", list_init = list_init_logit)
#'
#'
#' list_init_probit <- set_init(d, p, gam_vb, mu_beta_vb, sig2_beta_vb,
#'                              tau_vb = NULL, family = "binomial-probit")
#'
#' vb_probit <- locus(Y = dat_b$phenos, X = dat_b$snps, p0_av = p0_av,
#'                    family = "binomial-probit", list_init = list_init_probit)
#'
#' # Binary outcomes with covariates
#' list_init_logit_z <- set_init(d, p, gam_vb, mu_beta_vb, sig2_beta_vb_logit,
#'                               tau_vb = NULL, family = "binomial-logit",
#'                               n = n, q = q, mu_alpha_vb = mu_alpha_vb,
#'                               sig2_alpha_vb = sig2_alpha_vb)
#'
#' vb_logit_z <- locus(Y = dat_b$phenos, X = dat_b$snps, p0_av = p0_av, Z = Z,
#'                 family = "binomial-logit", list_init = list_init_logit_z)
#'
#'
#' list_init_probit_z <- set_init(d, p, gam_vb, mu_beta_vb, sig2_beta_vb,
#'                                tau_vb = NULL, family = "binomial-probit",
#'                                q = q, mu_alpha_vb = mu_alpha_vb,
#'                                sig2_alpha_vb = sig2_alpha_vb)
#'
#' vb_probit_z <- locus(Y = dat_b$phenos, X = dat_b$snps, p0_av = p0_av, Z = Z,
#'                 family = "binomial-probit", list_init = list_init_probit_z)
#'
#' @seealso  \code{\link{set_hyper}}, \code{\link{locus}}
#'
#' @export
#'
set_init <- function(d, p, gam_vb, mu_beta_vb, sig2_beta_vb, tau_vb,
                     family = "gaussian", n = NULL, q = NULL,
                     mu_alpha_vb = NULL, sig2_alpha_vb = NULL) {

  check_structure_(d, "vector", "numeric", 1)
  check_natural_(d)

  check_structure_(p, "vector", "numeric", 1)
  check_natural_(p)

  check_structure_(n, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(n)) check_natural_(n)

  check_structure_(q, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(q)) check_natural_(q)

  stopifnot(family %in% c("gaussian", "binomial-logit", "binomial-probit"))

  check_structure_(gam_vb, "matrix", "double", c(p, d))
  check_zero_one_(gam_vb)

  check_structure_(mu_beta_vb, "matrix", "double", c(p, d))


  if (family == "gaussian") {

    check_structure_(sig2_beta_vb, "vector", "double", d)

    check_structure_(tau_vb, "vector", "double", d)
    check_positive_(tau_vb)

    chi_vb <- NULL

  } else if (family == "binomial-logit"){

    check_structure_(sig2_beta_vb, "matrix", "double", c(p, d))

    if (!is.null(tau_vb))
      stop("tau_vb must be NULL for logistic regression.")

    if (is.null(n))
      stop("The number of observations, n, must be provided for logistic regression.")
    chi_vb <- matrix(1 / 2, nrow = n, ncol = d)

  } else {

    check_structure_(sig2_beta_vb, "vector", "double", d)

    if (!is.null(tau_vb))
      stop("tau_vb must be NULL for logistic regression.")

    chi_vb <- NULL

  }

  check_positive_(sig2_beta_vb)

  if (!is.null(q)) {

    check_structure_(mu_alpha_vb, "matrix", "double", c(q, d))

    check_structure_(sig2_alpha_vb, "matrix", "double", c(q, d))
    check_positive_(sig2_alpha_vb)

  } else if (!is.null(mu_alpha_vb) | !is.null(sig2_alpha_vb)) {

    stop(paste("Provided q = NULL, not consistent with mu_alpha_vb or ",
               "sig2_alpha_vb being non-null.", sep = ""))

  }

  d_init <- d
  p_init <- p
  q_init <- q
  n_init <- n

  family_init <- family

  list_init <- create_named_list_(d_init, n_init, p_init, q_init, family_init,
                                  gam_vb, mu_beta_vb, sig2_beta_vb, tau_vb,
                                  chi_vb, mu_alpha_vb, sig2_alpha_vb)

  class(list_init) <- "init"

  list_init
}


auto_set_init_ <- function(Y, p, p_star, q, user_seed, family) {

  d <- ncol(Y)
  n <- nrow(Y)

  if (!is.null(user_seed)) set.seed(user_seed)

  if (length(p_star) == 1) {
    shape1_gam <- 1
    p0 <- p_star
  } else {
    shape1_gam <- rep(1, p)
    p0 <- sum(p_star / p)
  }

  shape2_gam <- d * (p - p_star) / p_star

  gam_vb <- matrix(rbeta(p * d, shape1 = shape1_gam, shape2 = shape2_gam),
                   nrow = p)
  mu_beta_vb <- matrix(rnorm(p * d), nrow = p)


  sig2_inv_vb <- 1e-2

  if (family == "gaussian") {

    tau_vb <- 1 / median(apply(Y, 2, var))
    if (!is.finite(tau_vb)) tau_vb <- 1e3
    tau_vb <- rep(tau_vb, d)

    sig2_beta_vb <- 1 / rgamma(d, shape = 2, rate = 1 / (sig2_inv_vb * tau_vb))

    chi_vb <- NULL

  } else if (family == "binomial-logit") {

    sig2_beta_vb <- 1 / t(replicate(p, rgamma(d, shape = 2, rate = 1 / sig2_inv_vb)))

    chi_vb <- matrix(1 / 2, nrow = n, ncol = d)

    tau_vb <- NULL

  } else {

    sig2_beta_vb <- 1 / rgamma(d, shape = 2, rate = 1 / sig2_inv_vb)

    tau_vb <- chi_vb <- NULL

  }



  if (!is.null(q)) {

    mu_alpha_vb <- matrix(rnorm(q * d), nrow = q)

    if (family == "gaussian") {

      zeta2_inv_vb <- rgamma(q, shape = 1, rate = 1)

      sig2_alpha_vb <- 1 / sapply(tau_vb,
                                  function(tau_vb_t) {
                                    rgamma(q, shape = 2,
                                           rate = 1 / (zeta2_inv_vb * tau_vb_t))
                                  } )

    } else{

      zeta2_inv_vb <- matrix(rgamma(q * d, shape = 1, rate = 1), nrow = q)
      sig2_alpha_vb <- 1 / apply(zeta2_inv_vb, 2, function(zeta2_inv_vb_t) rgamma(q, shape = 2, rate = 1 / zeta2_inv_vb_t))

    }


  } else {

    mu_alpha_vb <- NULL
    sig2_alpha_vb <- NULL

  }


  d_init <- d
  p_init <- p
  q_init <- q
  n_init <- n


  family_init <- family

  list_init <- create_named_list_(d_init, n_init, p_init, q_init, family_init,
                                  gam_vb, mu_beta_vb, sig2_beta_vb, tau_vb,
                                  chi_vb, mu_alpha_vb, sig2_alpha_vb)

  class(list_init) <- "out_init"

  list_init
}
