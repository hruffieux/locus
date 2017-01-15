#' Gather model hyperparameters provided by the user.
#'
#' @export
feed_hyperparam <- function(d, p, eta, kappa, lambda, nu, a, b,
                            q = NULL, phi = NULL, xi = NULL) {

  check_structure_(eta, "vector", "double", c(1, d))
  check_positive_(eta)
  if (length(eta) == 1) eta <- rep(eta, d)

  check_structure_(kappa, "vector", "double", c(1,d))
  check_positive_(kappa)
  if (length(kappa) == 1) kappa <- rep(kappa, d)

  check_structure_(lambda, "vector", "double", 1)
  check_positive_(lambda)

  check_structure_(nu, "vector", "double", 1)
  check_positive_(nu)

  check_structure_(a, "vector", "double", c(1, p))
  check_positive_(a)
  if (length(a) == 1) a <- rep(a, p)

  check_structure_(b, "vector", "double", c(1, p))
  check_positive_(b)
  if (length(b) == 1) b <- rep(b, p)

  if (!is.null(q)) {

    check_structure_(phi, "vector", "double", c(1, q))
    check_positive_(phi)
    if (length(phi) == 1) phi <- rep(phi, q)

    check_structure_(xi, "vector", "double", c(1, q))
    check_positive_(xi)
    if (length(xi) == 1) xi <- rep(xi, q)

  } else if (!is.null(phi) | !is.null(xi)) {
    stop("Provided q = NULL, not consitent with phi or xi being non-null.")
  }

  d_hyper <- d
  p_hyper <- p
  q_hyper <- q

  list_hyper <- create_named_list_(d_hyper, p_hyper, q_hyper,
                                   eta, kappa, lambda, nu, a, b, phi, xi)

  class(list_hyper) <- "hyper"

  list_hyper

}


auto_set_hyperparam_ <- function(Y, p, p_star, q = NULL) {

  d <- ncol(Y)

  # hyperparameter set using the data Y
  eta <- 1 / median(apply(Y, 2, var)) #median to be consistent when doing permutations
  if (!is.finite(eta)) eta <- 1e3
  eta <- rep(eta, d)
  kappa <- rep(1, d)

  # if p_star is of length 1, p_star is the guessed number of active predictors
  # else (p_star is of length p), p_star / p is the vector of guessed
  # probabilities that each predictor is active and the sum of its entries is
  # the corresponding guessed number of active predictors
  if (length(p_star) == 1) p0 <- p_star
  else p0 <- sum(p_star / p)

  lambda <- 1e-2
  nu <- 1

  a <- rep(1, p)
  b <- d * (p - p_star) / p_star
  if (length(b) == 1) b <- rep(b, p)

  # hyperparameters of beta distributions
  check_positive_(a)
  check_positive_(b)

  if (!is.null(q)) {
    phi <- rep(1, q)
    xi <- rep(1, q)
  } else {
    phi <- NULL
    xi <- NULL
  }

  d_hyper <- d
  p_hyper <- p
  q_hyper <- q

  list_hyper <- create_named_list_(d_hyper, p_hyper, q_hyper,
                                   eta, kappa, lambda, nu, a, b, phi, xi)

  class(list_hyper) <- "out_hyper"

  list_hyper

}

#' Gather initial variational parameters provided by the user.
#'
#' @export
feed_init_param <- function(d, p, gam_vb, mu_beta_vb, sig2_beta_vb,
                            tau_vb, q = NULL, mu_alpha_vb = NULL,
                            sig2_alpha_vb = NULL, zeta2_inv_vb = NULL) {

  check_structure_(gam_vb, "matrix", "double", c(p, d))
  check_zero_one_(gam_vb)

  check_structure_(mu_beta_vb, "matrix", "double", c(p, d))

  # check_structure_(sig2_beta_vb, "matrix", "double", c(p, d))
  # check_positive_(sig2_beta_vb)

  check_structure_(sig2_beta_vb, "vector", "double", d)

  check_positive_(sig2_beta_vb)

  check_structure_(tau_vb, "vector", "double", d)
  check_positive_(tau_vb)


  if (!is.null(q)) {

    check_structure_(mu_alpha_vb, "matrix", "double", c(q, d))

    check_structure_(sig2_alpha_vb, "matrix", "double", c(q, d))
    check_positive_(sig2_alpha_vb)

    check_structure_(zeta2_inv_vb, "vector", "double", q)
    check_positive_(zeta2_inv_vb)

  } else if (!is.null(mu_alpha_vb) | !is.null(sig2_alpha_vb) | !is.null(zeta2_inv_vb) ) {
    stop(paste("Provided q = NULL, not consistent with mu_alpha_vb, sig2_alpha_vb",
               "or zeta2_inv_vb being non-null.", sep = ""))
  }

  d_init <- d
  p_init <- p
  q_init <- q

  list_init <- create_named_list_(d_init, p_init, q_init, gam_vb, mu_beta_vb,
                                  sig2_beta_vb, tau_vb, mu_alpha_vb,
                                  sig2_alpha_vb, zeta2_inv_vb)

  class(list_init) <- "init"

  list_init
}


auto_init_param_ <- function(Y, p, p_star, user_seed, q = NULL) {

  d <- ncol(Y)

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


  tau_vb <- 1 / median(apply(Y, 2, var))
  if (!is.finite(tau_vb)) tau_vb <- 1e3
  tau_vb <- rep(tau_vb, d)

  sig2_inv_vb <- 1e-2
  sig2_beta_vb <- 1 / rgamma(d, shape = 2, rate = 1 / (sig2_inv_vb * tau_vb))


  if (!is.null(q)) {
    mu_alpha_vb <- matrix(rnorm(q * d), nrow = q)
    zeta2_inv_vb <- rgamma(q, shape = 1, rate = 1)
    sig2_alpha_vb <- 1 / sapply(tau_vb,
                                function(tau_vb_t) {
                                  rgamma(q, shape = 2,
                                         rate = 1 / (zeta2_inv_vb * tau_vb_t))
                                } )
  } else {
    mu_alpha_vb <- NULL
    sig2_alpha_vb <- NULL
    zeta2_inv_vb <- NULL
  }


  d_init <- d
  p_init <- p
  q_init <- q

  list_init <- create_named_list_(d_init, p_init, q_init, gam_vb, mu_beta_vb,
                                  sig2_beta_vb, tau_vb, mu_alpha_vb,
                                  sig2_alpha_vb, zeta2_inv_vb)

  class(list_init) <- "out_init"

  list_init
}
