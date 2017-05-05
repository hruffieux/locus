# This file is part of the `locus` R package:
#     https://github.com/hruffieux/locus
#

#' Gather model hyperparameters provided by the user.
#'
#' This function must be used to provide hyperparameter values for the model
#' used in \code{\link{locus}}.
#'
#' The \code{\link{locus}} function can also be used with default hyperparameter
#' choices (without using \code{\link{set_hyper}}) by setting its argument
#' \code{list_hyper} to \code{NULL}.
#'
#' @param d Number of responses.
#' @param p Number of candidate predictors.
#' @param lambda Vector of length 1 providing the values of hyperparameter
#'   \eqn{\lambda} for the prior distribution of \eqn{\sigma^{-2}}. \eqn{\sigma}
#'   represents the typical size of nonzero effects.
#' @param nu Vector of length 1 providing the values of hyperparameter \eqn{\nu}
#'   for the prior distribution of \eqn{\sigma^{-2}}. \eqn{\sigma} represents
#'   the typical size of nonzero effects.
#' @param a Vector of length 1 or p providing the values of hyperparameter
#'   \eqn{a} for the prior distributions for the proportion of responses
#'   associated with each candidate predictor, \eqn{\omega} (vector of length p).
#'   If of length 1, the provided value is repeated p times.
#' @param b Vector of length 1 or p providing the values of hyperparameter
#'   \eqn{b} for the prior distributions for the proportion of responses
#'   associated with each candidate predictor, \eqn{\omega} (vector of length p).
#'   If of length 1, the provided value is repeated p times.
#' @param eta Vector of length 1 or d for \code{link = "identity"}, and of
#'   length 1 or d_cont = d - length(ind_bin) (the number of continuous response
#'   variables) for \code{link = "mix"}. Provides the values of
#'   hyperparameter \eqn{\eta} for the prior distributions of the continuous
#'   response residual precisions, \eqn{\tau}. If of length 1, the provided
#'   value is repeated d, resp. d_cont, times. Must be \code{NULL} for
#'   \code{link = "logit"} and \code{link = "probit"}.
#' @param kappa Vector of length 1 or d for \code{link = "identity"}, and of
#'   length 1 or d_cont = d - length(ind_bin) (the number of continuous response
#'   variables) for \code{link = "mix"}. Provides the values of hyperparameter
#'   \eqn{\kappa} for the prior distributions of the response residual
#'   precisions, \eqn{\tau}. If of length 1, the provided value is repeated d,
#'   resp. d_cont, times. Must be \code{NULL} for \code{link = "logit"} and
#'   \code{link = "probit"}.
#' @param link Response link. Must be "\code{identity}" for linear regression,
#'   "\code{logit}" for logistic regression, "\code{probit}"
#'   for probit regression, or "\code{mix}" for a mix of identity and probit
#'   link functions (in this case, the indices of the binary responses must be
#'   gathered in argument \code{ind_bin}, see below).
#' @param ind_bin If \code{link = "mix"}, vector of indices corresponding to
#'   the binary variables in \code{Y}. Must be \code{NULL} if
#'   \code{link != "mix"}.
#' @param q Number of covariates. Default is \code{NULL}, for \code{Z}
#'   \code{NULL}.
#' @param phi Vector of length 1 or q providing the values of hyperparameter
#'   \eqn{\phi} for the prior distributions for the sizes of the nonzero
#'   covariate effects, \eqn{\zeta}. If of length 1, the provided value is
#'   repeated q times. Default is \code{NULL}, for \code{Z} \code{NULL}.
#' @param xi Vector of length 1 or q providing the values of hyperparameter
#'   \eqn{\xi} for the prior distributions for the sizes of the nonzero
#'   covariate effects, \eqn{\zeta}. If of length 1, the provided value is
#'   repeated q times. Default is \code{NULL}, for \code{Z} \code{NULL}.
#' @param r Number of variables representing external information on the
#'   candidate predictors. Default is \code{NULL}, for \code{V} \code{NULL}.
#' @param m0 Vector of length 1 or p providing the values of hyperparameter
#'   \eqn{m0} for the prior distribution of \eqn{c0} linked to the proportion of
#'   responses associated with each candidate predictor. If of length 1, the
#'   provided value is repeated p times. Default is \code{NULL}, for \code{V}
#'   \code{NULL}.
#'
#' @return An object of class "\code{hyper}" preparing user hyperparameter in a
#'   form that can be passed to the \code{\link{locus}} function.
#'
#' @examples
#' seed <- 123; set.seed(seed)
#'
#' ###################
#' ## Simulate data ##
#' ###################
#'
#' ## Examples using small problem sizes:
#' ##
#' n <- 200; p <- 200; p0 <- 20; d <- 20; d0 <- 15; q <- 2; r <- 3
#'
#' ## Candidate predictors (subject to selection)
#' ##
#' # Here we simulate common genetic variants (but any type of candidate
#' # predictors can be supplied).
#' # 0 = homozygous, major allele, 1 = heterozygous, 2 = homozygous, minor allele
#'
#' X_act <- matrix(rbinom(n * p0, size = 2, p = 0.25), nrow = n)
#' X_inact <- matrix(rbinom(n * (p - p0), size = 2, p = 0.25), nrow = n)
#'
#' shuff_x_ind <- sample(p)
#' X <- cbind(X_act, X_inact)[, shuff_x_ind]
#'
#' bool_x_act <- shuff_x_ind <= p0
#'
#' pat_act <- beta <- matrix(0, nrow = p0, ncol = d0)
#' pat_act[sample(p0*d0, floor(p0*d0/5))] <- 1
#' beta[as.logical(pat_act)] <-  rnorm(sum(pat_act))
#'
#' ## Covariates (not subject to selection)
#' ##
#' Z <- matrix(rnorm(n * q), nrow = n)
#'
#' alpha <-  matrix(rnorm(q * d), nrow = q)
#'
#' ## Gaussian responses
#' ##
#' Y_act <- matrix(rnorm(n * d0, mean = X_act %*% beta, sd = 0.5), nrow = n)
#' Y_inact <- matrix(rnorm(n * (d - d0), sd = 0.5), nrow = n)
#' shuff_y_ind <- sample(d)
#' Y <- cbind(Y_act, Y_inact)[, shuff_y_ind] + Z %*% alpha
#'
#' ## Binary responses
#' ##
#' Y_bin <- ifelse(Y > 0, 1, 0)
#' ## Informative annotation variables
#' ##
#' V <- matrix(rnorm(p * r), nrow = p)
#' V[bool_x_act, ] <- rnorm(p0 * r, mean = 2)
#'
#' ########################
#' ## Infer associations ##
#' ########################
#'
#' ## Continuous responses
#' ##
#'
#' # No covariate
#' #
#' # a and b chosen so that the prior mean number of responses associated with
#' # each candidate predictor is 1/4.
#' list_hyper_g <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
#'                           eta = 1, kappa = apply(Y, 2, var),
#'                           link = "identity")
#'
#' # We take p0_av = p0 (known here); this choice may result in variable
#' # selections that are (too) conservative in some cases. In practice, it is
#' # advised to set p0_av as a slightly overestimated guess of p0, or perform
#' # cross-validation using function `set_cv'.
#'
#' vb_g <- locus(Y = Y, X = X, p0_av = p0, link = "identity",
#'               list_hyper = list_hyper_g, user_seed = seed)
#'
#' # With covariates
#' #
#' list_hyper_g_z <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
#'                             eta = 1, kappa = apply(Y, 2, var),
#'                             link = "identity", q = q, phi = 1, xi = 1)
#'
#' vb_g_z <- locus(Y = Y, X = X, p0_av = p0, Z = Z, link = "identity",
#'                 list_hyper = list_hyper_g_z, user_seed = seed)
#'
#'
#' # With external annotation variables
#' #
#' list_hyper_g_v <- set_hyper(d, p, lambda = 1, nu = 1, a = NULL, b = NULL,
#'                             eta = 1, kappa = apply(Y, 2, var),
#'                             link = "identity", r = r, m0 = 0)
#'
#' vb_g_v <- locus(Y = Y, X = X, p0_av = p0,  V = V, link = "identity",
#'                 list_hyper = list_hyper_g_v, user_seed = seed)
#'
#' ## Binary responses
#' ##
#' list_hyper_logit <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
#'                               eta = NULL, kappa = NULL, link = "logit",
#'                               q = q, phi = 1, xi = 1)
#'
#' vb_logit <- locus(Y = Y_bin, X = X, p0_av = p0, Z = Z, link = "logit",
#'                   list_hyper = list_hyper_logit, user_seed = seed)
#'
#' list_hyper_probit <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
#'                                eta = NULL, kappa = NULL, link = "probit",
#'                                q = q, phi = 1, xi = 1)
#'
#' vb_probit <- locus(Y = Y_bin, X = X, p0_av = p0, Z = Z, link = "probit",
#'                    list_hyper = list_hyper_probit, user_seed = seed)
#'
#'
#' ## Mix of continuous and binary responses
#' ##
#' Y_mix <- cbind(Y, Y_bin)
#' ind_bin <- (d+1):(2*d)
#'
#' list_hyper_mix <- set_hyper(2*d, p, lambda = 1, nu = 1, a = 1, b = 8*d-1,
#'                             eta = 1, kappa = apply(Y, 2, var), link = "mix",
#'                             ind_bin = ind_bin, q = q, phi = 1, xi = 1)
#'
#' vb_mix <- locus(Y = Y_mix, X = X, p0_av = p0, Z = Z, link = "mix",
#'                 ind_bin = ind_bin, list_hyper = list_hyper_mix,
#'                 user_seed = seed)
#'
#' @seealso  \code{\link{set_init}}, \code{\link{locus}}
#'
#' @export
#'
set_hyper <- function(d, p, lambda, nu, a, b, eta, kappa, link = "identity",
                      ind_bin = NULL, q = NULL, phi = NULL, xi = NULL,
                      r = NULL, m0 = NULL, G = NULL) {

  check_structure_(d, "vector", "numeric", 1)
  check_natural_(d)

  check_structure_(p, "vector", "numeric", 1)
  check_natural_(p)

  check_structure_(q, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(q)) check_natural_(q)

  check_structure_(r, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(r)) check_natural_(r)

  check_structure_(G, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(G)) check_natural_(G)

  stopifnot(link %in% c("identity", "logit", "probit", "mix"))

  if(!is.null(G) && (link != "identity" | !is.null(q) | !is.null(r)))
    stop("Group selection(G non-NULL) implemented only for identity link, Z = NULL and V = NULL. Exit.")

  ind_bin <- prepare_ind_bin_(d, ind_bin, link)

  if (is.null(r)) {

    if (is.null(G)) {

      check_structure_(a, "vector", "double", c(1, p))
      if (length(a) == 1) a <- rep(a, p)

      check_structure_(b, "vector", "double", c(1, p))
      if (length(b) == 1) b <- rep(b, p)

    } else {

      check_structure_(a, "vector", "double", c(1, G))
      if (length(a) == 1) a <- rep(a, G)

      check_structure_(b, "vector", "double", c(1, G))
      if (length(b) == 1) b <- rep(b, G)

    }

    check_positive_(a)
    check_positive_(b)


    s02 <- s2 <- NULL

    if (!is.null(m0))
      stop("Provided r = NULL, not consitent with m0 being non-null.")

  } else {

    check_structure_(m0, "vector", "double", c(1, p))
    if (length(m0) == 1) m0 <- rep(m0, p)

    # prior info
    s02 <- 0.1 # prior variance for the intercept, bernoulli-probit
    s2 <- 1e-2 # prior variance for external info coefficients
    # (effects likely to be concentrated around zero)

    if (!is.null(a) | !is.null(b))
      stop("Provided r != NULL, not consitent with a and b being non-null.")

  }

  check_structure_(lambda, "vector", "double", 1)
  check_positive_(lambda)

  check_structure_(nu, "vector", "double", 1)
  check_positive_(nu)

  if (link %in% c("identity", "mix")) {

    d_cont <- d - length(ind_bin) # length(NULL) = 0 for link = "identity"

    check_structure_(eta, "vector", "double", c(1, d_cont))
    check_positive_(eta)
    if (length(eta) == 1) eta <- rep(eta, d_cont)

    check_structure_(kappa, "vector", "double", c(1, d_cont))
    check_positive_(kappa)
    if (length(kappa) == 1) kappa <- rep(kappa, d_cont)

  } else {

    if (!is.null(eta) | !is.null(kappa))
      stop("Both eta and kappa must be NULL for logistic and probit regression.")

  }

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
  r_hyper <- r
  G_hyper <- G

  ind_bin_hyper <- ind_bin

  link_hyper <- link

  list_hyper <- create_named_list_(d_hyper, G_hyper, p_hyper, q_hyper, r_hyper,
                                   link_hyper, ind_bin_hyper, eta, kappa,
                                   lambda, nu, a, b, phi, xi, m0, s02, s2)

  class(list_hyper) <- "hyper"

  list_hyper

}


# Internal function setting default model hyperparameters when not provided by
# the user.
#
auto_set_hyper_ <- function(Y, G, p, p_star, q, r, link, ind_bin) {

  d <- ncol(Y)

  lambda <- 1e-2
  nu <- 1

  if (link %in% c("identity", "mix")) {

    if (link == "mix") Y <- Y[, -ind_bin, drop = FALSE]

    d_cont <- d - length(ind_bin) # length(NULL) = 0 for link = "identity"

    # hyperparameter set using the data Y
    eta <- 1 / median(apply(Y, 2, var)) #median to be consistent when doing permutations
    if (!is.finite(eta)) eta <- 1e3
    eta <- rep(eta, d_cont)
    kappa <- rep(1, d_cont)

  } else {

    eta <- kappa <- NULL

  }

  if (is.null(r)) {

    b <- d * (p - p_star) / p_star

    if(is.null(G)) {
      a <- rep(1, p)
      if (length(b) == 1) b <- rep(b, p)
    } else {
      a <- rep(1, G)
      if (length(b) == 1) b <- rep(b, G)
    }

    # hyperparameters of beta distributions
    check_positive_(a)
    check_positive_(b)

    m0 <- s02 <- s2 <- NULL

  } else {

    # hyperparameters external info model
    s02 <- 0.1 # prior variance for the intercept, bernoulli-probit
    s2 <- 1e-2 # prior variance for external info coefficients (effects likely to be concentrated around zero)

    up <- 1e6 # we solve the equation below numerically as no closed form exists.
    m0 <- - sapply(p_star, function(ps) uniroot(function(x)    # sparsity control, m0 = - m0_star
      pnorm(x / sqrt(1 + d*s02)) * (pnorm(x))^(d-1) - (p - ps)/p,
      interval = c(-up, up))$root)

    if (length(m0) == 1) m0 <- rep(m0, p)

    a <- b <- NULL

  }

  if (!is.null(q)) {

    phi <- xi <- rep(1, q)

  } else {

    phi <- xi <- NULL

  }

  d_hyper <- d
  G_hyper <- G
  p_hyper <- p
  q_hyper <- q
  r_hyper <- r

  ind_bin_hyper <- ind_bin

  link_hyper <- link

  list_hyper <- create_named_list_(d_hyper, G_hyper, p_hyper, q_hyper, r_hyper,
                                   link_hyper, ind_bin_hyper, eta, kappa,
                                   lambda, nu, a, b, phi, xi, m0, s02, s2)

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
#'   variational parameter yielding regression coefficient estimates for
#'   predictor-response pairs included in the model.
#' @param sig2_beta_vb Vector of length d, for \code{link = "identity"} and
#'   for \code{link = "mix"}, of length 1 for \code{link = "probit"}, and a
#'   matrix of size p x d, for \code{link = "logit"}, with initial values for
#'   the variational parameter yielding estimates of effect variances for
#'   predictor-response pairs included in the model. For
#'   \code{link = "identity"} and \code{link = "mix"}, these values are the same
#'   for all the predictors (as a result of the predictor variables being
#'   standardized before the variational algorithm). For \code{link = "probit"},
#'   they are the same for all the predictors and responses.
#' @param tau_vb Vector of length d, for \code{link = "identity"}, and of
#'   length d_cont = d - length(ind_bin) (number of continuous responses), for
#'   \code{link = "mix"}, with initial values for the variational parameter
#'   yielding estimates for the continuous response residual precisions. Must be
#'   \code{NULL} for \code{link = "logit"} and \code{link = "probit"}.
#' @param link Response link. Must be "\code{identity}" for linear regression,
#'   "\code{logit}" for logistic regression, "\code{probit}" for probit
#'   regression, or "\code{mix}" for a mix of identity and probit link functions
#'   (in this case, the indices of the binary responses must be gathered in
#'   argument \code{ind_bin}, see below).
#' @param ind_bin If \code{link = "mix"}, vector of indices corresponding to the
#'   binary variables in \code{Y}. Must be \code{NULL} if \code{link != "mix"}.
#' @param q Number of covariates. Default is \code{NULL}, for \code{Z}
#'   \code{NULL}.
#' @param mu_alpha_vb Matrix of size q x d with initial values for the
#'   variational parameter yielding regression coefficient estimates for
#'   covariate-response pairs. Default is \code{NULL}, for \code{Z} \code{NULL}.
#' @param sig2_alpha_vb Matrix of size q x d for \code{link = "identity"},
#'   for \code{link = "logit"} and for \code{link = "mix"} with initial values
#'   for the variational parameter yielding estimates of effect variances for
#'   covariate-response pairs. Vector of length q for \code{link = "probit"}.
#'   Default is \code{NULL}, for \code{Z} \code{NULL}.
#' @param r Number of variables representing external information on the
#'   candidate predictors. Default is \code{NULL}, for \code{V} \code{NULL}.
#' @param mu_c0_vb Vector of length p with initial values for the variational
#'   parameter linked to the proportion of responses associated with each
#'   candidate predictor. Default is \code{NULL}, for \code{V} \code{NULL}.
#' @param mu_c_vb Matrix of size r x d with initial values for the variational
#'   parameter yielding regression coefficient estimates for the influence of
#'   external information on the candidate predictors on their selection.
#'   Default is \code{NULL}, for \code{V} \code{NULL}.
#'
#'
#' @return An object of class "\code{init}" preparing user initial values for
#'   the variational parameters in a form that can be passed to the
#'   \code{\link{locus}} function.
#'
#' @examples
#' seed <- 123; set.seed(seed)
#'
#' ###################
#' ## Simulate data ##
#' ###################
#'
#' ## Examples using small problem sizes:
#' ##
#' n <- 200; p <- 200; p0 <- 20; d <- 20; d0 <- 15; q <- 2; r <- 3
#'
#' ## Candidate predictors (subject to selection)
#' ##
#' # Here we simulate common genetic variants (but any type of candidate
#' # predictors can be supplied).
#' # 0 = homozygous, major allele, 1 = heterozygous, 2 = homozygous, minor allele
#'
#' X_act <- matrix(rbinom(n * p0, size = 2, p = 0.25), nrow = n)
#' X_inact <- matrix(rbinom(n * (p - p0), size = 2, p = 0.25), nrow = n)
#'
#' shuff_x_ind <- sample(p)
#' X <- cbind(X_act, X_inact)[, shuff_x_ind]
#'
#' bool_x_act <- shuff_x_ind <= p0
#'
#' pat_act <- beta <- matrix(0, nrow = p0, ncol = d0)
#' pat_act[sample(p0*d0, floor(p0*d0/5))] <- 1
#' beta[as.logical(pat_act)] <-  rnorm(sum(pat_act))
#'
#' ## Covariates (not subject to selection)
#' ##
#' Z <- matrix(rnorm(n * q), nrow = n)
#'
#' alpha <-  matrix(rnorm(q * d), nrow = q)
#'
#' ## Gaussian responses
#' ##
#' Y_act <- matrix(rnorm(n * d0, mean = X_act %*% beta, sd = 0.5), nrow = n)
#' Y_inact <- matrix(rnorm(n * (d - d0), sd = 0.5), nrow = n)
#' shuff_y_ind <- sample(d)
#' Y <- cbind(Y_act, Y_inact)[, shuff_y_ind] + Z %*% alpha
#'
#' ## Binary responses
#' ##
#' Y_bin <- ifelse(Y > 0, 1, 0)
#' ## Informative annotation variables
#' ##
#' V <- matrix(rnorm(p * r), nrow = p)
#' V[bool_x_act, ] <- rnorm(p0 * r, mean = 2)
#'
#' ########################
#' ## Infer associations ##
#' ########################
#'
#' ## Continuous responses
#' ##
#'
#' # No covariate
#' #
#' # gam_vb chosen so that the prior mean number of responses associated with
#' # each candidate predictor is 1/4.
#' gam_vb <- matrix(rbeta(p * d, shape1 = 1, shape2 = 4*d-1), nrow = p)
#' mu_beta_vb <- matrix(rnorm(p * d), nrow = p)
#' tau_vb <- 1 / apply(Y, 2, var)
#' sig2_beta_vb <- 1 / rgamma(d, shape = 2, rate = 1 / tau_vb)
#'
#' list_init_g <- set_init(d, p, gam_vb, mu_beta_vb, sig2_beta_vb, tau_vb,
#'                         link = "identity")
#'
#' # We take p0_av = p0 (known here); this choice may result in variable
#' # selections that are (too) conservative in some cases. In practice, it is
#' # advised to set p0_av as a slightly overestimated guess of p0, or perform
#' # cross-validation using function `set_cv'.
#'
#' vb_g <- locus(Y = Y, X = X, p0_av = p0, link = "identity",
#'               list_init = list_init_g)
#'
#' # With covariates
#' #
#' mu_alpha_vb <- matrix(rnorm(q * d), nrow = q)
#' sig2_alpha_vb <- 1 / matrix(rgamma(q * d, shape = 2, rate = 1), nrow = q)
#'
#' list_init_g_z <- set_init(d, p, gam_vb, mu_beta_vb, sig2_beta_vb, tau_vb,
#'                           link = "identity", q = q,
#'                           mu_alpha_vb = mu_alpha_vb,
#'                           sig2_alpha_vb = sig2_alpha_vb)
#'
#' vb_g_z <- locus(Y = Y, X = X, p0_av = p0, Z = Z, link = "identity",
#'                 list_init = list_init_g_z)
#'
#' # With external annotation variables
#' #
#' mu_c0_vb <- rnorm(p, mean = -1)
#' mu_c_vb <- matrix(rnorm(r * d, mean = 0, sd = 0.01), nrow = r)
#'
#' list_init_g_v <- set_init(d, p, gam_vb, mu_beta_vb, sig2_beta_vb, tau_vb,
#'                           link = "identity", r = r, mu_c0_vb = mu_c0_vb,
#'                           mu_c_vb = mu_c_vb)
#'
#' vb_g_v <- locus(Y = Y, X = X, p0_av = p0,  V = V, link = "identity",
#'                 list_init = list_init_g_v)
#'
#' ## Binary responses
#' ##
#' # gam_vb chosen so that the prior mean number of responses associated with
#' # each candidate predictor is 1/4.
#' sig2_beta_vb_logit <- 1 / t(replicate(p, rgamma(d, shape = 2, rate = 1)))
#'
#' list_init_logit <- set_init(d, p, gam_vb, mu_beta_vb, sig2_beta_vb_logit,
#'                             tau_vb = NULL, link = "logit", q = q,
#'                             mu_alpha_vb = mu_alpha_vb,
#'                             sig2_alpha_vb = sig2_alpha_vb)
#'
#' vb_logit <- locus(Y = Y_bin, X = X, p0_av = p0, Z = Z, link = "logit",
#'                   list_init = list_init_logit)
#'
#' sig2_alpha_vb_probit <- sig2_alpha_vb[, 1]
#' sig2_beta_vb_probit <- sig2_beta_vb[1]
#' list_init_probit <- set_init(d, p, gam_vb, mu_beta_vb, sig2_beta_vb_probit,
#'                              tau_vb = NULL, link = "probit", q = q,
#'                              mu_alpha_vb = mu_alpha_vb,
#'                              sig2_alpha_vb = sig2_alpha_vb_probit)
#'
#' vb_probit <- locus(Y = Y_bin, X = X, p0_av = p0, Z = Z, link = "probit",
#'                    list_init = list_init_probit)
#'
#' ## Mix of continuous and binary responses
#' ##
#' Y_mix <- cbind(Y, Y_bin)
#' ind_bin <- (d+1):(2*d)
#'
#' # gam_vb chosen so that the prior mean number of responses associated with
#' # each candidate predictor is 1/4.
#' gam_vb_mix <- matrix(rbeta(p * 2*d, shape1 = 1, shape2 = 8*d-1), nrow = p)
#' mu_beta_vb_mix <- matrix(rnorm(p * 2*d), nrow = p)
#' sig2_beta_vb_mix <- 1 / c(rgamma(d, shape = 2, rate = 1 / tau_vb),
#'                           rgamma(d, shape = 2, rate = 1))
#' mu_alpha_vb_mix <- matrix(rnorm(q * 2*d), nrow = q)
#' sig2_alpha_vb_mix <- 1 / matrix(rgamma(q * 2*d, shape = 2, rate = 1), nrow = q)
#'
#' list_init_mix <- set_init(2*d, p, gam_vb_mix, mu_beta_vb_mix,
#'                           sig2_beta_vb_mix, tau_vb, link = "mix",
#'                           ind_bin = ind_bin, q = q,
#'                           mu_alpha_vb = mu_alpha_vb_mix,
#'                           sig2_alpha_vb = sig2_alpha_vb_mix)
#'
#' vb_mix <- locus(Y = Y_mix, X = X, p0_av = p0, Z = Z, link = "mix",
#'                 ind_bin = ind_bin, list_init = list_init_mix)
#'
#' @seealso  \code{\link{set_hyper}}, \code{\link{locus}}
#'
#' @export
#'
set_init <- function(d, p, gam_vb, mu_beta_vb, sig2_beta_vb, tau_vb,
                     link = "identity", ind_bin = NULL, q = NULL,
                     mu_alpha_vb = NULL, sig2_alpha_vb = NULL, r = NULL,
                     mu_c0_vb = NULL, mu_c_vb = NULL, sig2_inv_vb = NULL,
                     G = NULL) {

  check_structure_(d, "vector", "numeric", 1)
  check_natural_(d)

  check_structure_(p, "vector", "numeric", 1)
  check_natural_(p)

  check_structure_(q, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(q)) check_natural_(q)

  check_structure_(r, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(r)) check_natural_(r)

  check_structure_(G, "vector", "numeric", p, null_ok = TRUE)
  if (!is.null(G)) check_natural_(G)

  stopifnot(link %in% c("identity", "logit", "probit", "mix"))

  if(!is.null(G) && (link != "identity" | !is.null(q) | !is.null(r)))
    stop("Group selection (G non-NULL) implemented only for identity link, Z = NULL and V = NULL. Exit.")

  ind_bin <- prepare_ind_bin_(d, ind_bin, link)


  if(is.null(G)) {

    check_structure_(gam_vb, "matrix", "double", c(p, d))

  } else {

    check_structure_(gam_vb, "matrix", "double", c(G, d))

  }
  check_zero_one_(gam_vb)

  check_structure_(mu_beta_vb, "matrix", "double", c(p, d))

  if (link %in% c("identity", "mix")) {

    if(is.null(G)) {

      check_structure_(sig2_beta_vb, "vector", "double", d)
      if (!is.null(sig2_inv_vb))
        stop("sig2_inv_vb is not used when G is NULL.")

    } else {

      check_structure_(sig2_inv_vb, "vector", "double", d)
      if (!is.null(sig2_beta_vb))
        stop("sig2_beta_vb is not used when G is non-NULL.")

    }

    d_cont <- d - length(ind_bin) # length(NULL) = 0 for link = "identity"

    check_structure_(tau_vb, "vector", "double", d_cont)
    check_positive_(tau_vb)

    if (link == "mix") {
      tmp_tau_vb <- tau_vb
      tau_vb <- rep(1, d) # tau_vb is set to 1 for binary responses.
      tau_vb[-ind_bin] <- tmp_tau_vb
      rm(tmp_tau_vb)
    }

  } else if (link == "logit"){

    check_structure_(sig2_beta_vb, "matrix", "double", c(p, d))

    if (!is.null(sig2_inv_vb))
      stop("sig2_inv_vb is not used when G is NULL.")

    if (!is.null(tau_vb))
      stop("tau_vb must be NULL for logistic regression.")

  } else {

    check_structure_(sig2_beta_vb, "vector", "double", 1)

    if (!is.null(sig2_inv_vb))
      stop("sig2_inv_vb is not used when G is NULL.")

    if (!is.null(tau_vb))
      stop("tau_vb must be NULL for probit regression.")

  }

  check_positive_(sig2_beta_vb)

  if (!is.null(q)) {

    check_structure_(mu_alpha_vb, "matrix", "double", c(q, d))

    if (link == "probit") {

      check_structure_(sig2_alpha_vb, "vector", "double", q)

    } else {

      check_structure_(sig2_alpha_vb, "matrix", "double", c(q, d))

    }

    check_positive_(sig2_alpha_vb)

  } else if (!is.null(mu_alpha_vb) | !is.null(sig2_alpha_vb)) {

    stop(paste("Provided q = NULL, not consistent with mu_alpha_vb or ",
               "sig2_alpha_vb being non-null.", sep = ""))

  }


  if (!is.null(r)) {

    check_structure_(mu_c0_vb, "vector", "double", p)
    check_structure_(mu_c_vb, "matrix", "double", c(r, d))

  } else if (!is.null(mu_c0_vb) | !is.null(mu_c_vb)) {

    stop("Provided r = NULL, not consistent with mu_c0_vb or mu_c_vb being non-null.")

  }


  d_init <- d
  G_init <- G
  p_init <- p
  q_init <- q
  r_init <- r
  ind_bin_init <- ind_bin

  link_init <- link

  list_init <- create_named_list_(d_init, G_init, p_init, q_init, r_init, link_init,
                                  ind_bin_init, gam_vb, mu_beta_vb,
                                  sig2_beta_vb, sig2_inv_vb, tau_vb, mu_alpha_vb,
                                  sig2_alpha_vb, mu_c0_vb, mu_c_vb)

  class(list_init) <- "init"

  list_init
}


# Internal function setting default starting values when not provided by the user.
#
auto_set_init_ <- function(Y, G, p, p_star, q, r, user_seed, link, ind_bin) {

  d <- ncol(Y)

  if (!is.null(user_seed)) set.seed(user_seed)

  shape1_gam <- 1
  shape2_gam <- d * (p - p_star) / p_star


  if (is.null(G)) {

    if (length(p_star) > 1)
      shape1_gam <- rep(shape1_gam, p)

    gam_vb <- matrix(rbeta(p * d, shape1 = shape1_gam, shape2 = shape2_gam),
                     nrow = p)

  } else {

    if (length(p_star) > 1)
      shape1_gam <- rep(shape1_gam, G)

    gam_vb <- matrix(rbeta(G * d, shape1 = shape1_gam, shape2 = shape2_gam),
                     nrow = G)

  }

  mu_beta_vb <- matrix(rnorm(p * d), nrow = p)

  sig2_inv_vb <- 1e-2

  if (link %in% c("identity", "mix")) {

    if (link == "mix") Y <- Y[, -ind_bin, drop = FALSE]

    d_cont <- d - length(ind_bin) # length(NULL) = 0 for link = "identity"

    tau_vb <- 1 / median(apply(Y, 2, var))
    if (!is.finite(tau_vb)) tau_vb <- 1e3
    tau_vb <- rep(tau_vb, d_cont)

    if (link == "mix") {

      tmp_tau_vb <- tau_vb
      tau_vb <- rep(1, d) # tau_vb is set to 1 for binary responses.
      tau_vb[-ind_bin] <- tmp_tau_vb
      rm(tmp_tau_vb)

    }

    if (is.null(G)) {
      sig2_beta_vb <- 1 / rgamma(d, shape = 2, rate = 1 / (sig2_inv_vb * tau_vb))
    } else {
      sig2_beta_vb <- NULL
    }

  } else if (link == "logit") {

    sig2_beta_vb <- 1 / t(replicate(p, rgamma(d, shape = 2, rate = 1 / sig2_inv_vb)))

    tau_vb <- NULL

  } else {

    sig2_beta_vb <- 1 / rgamma(1, shape = 2, rate = 1 / sig2_inv_vb)

    tau_vb <- NULL

  }



  if (!is.null(q)) {

    mu_alpha_vb <- matrix(rnorm(q * d), nrow = q)

    if (link %in% c("identity", "mix")) {

      zeta2_inv_vb <- rgamma(q, shape = 1, rate = 1)

      sig2_alpha_vb <- 1 / sapply(tau_vb,
                                  function(tau_vb_t) {
                                    rgamma(q, shape = 2,
                                           rate = 1 / (zeta2_inv_vb * tau_vb_t))
                                  } )

    } else if (link == "logit"){

      zeta2_inv_vb <- matrix(rgamma(q * d, shape = 1, rate = 1), nrow = q)
      sig2_alpha_vb <- 1 / apply(zeta2_inv_vb, 2, function(zeta2_inv_vb_t) rgamma(q, shape = 2, rate = 1 / zeta2_inv_vb_t))

    } else {

      zeta2_inv_vb <- rgamma(q, shape = 1, rate = 1)
      sig2_alpha_vb <- 1 / rgamma(q, shape = 2, rate = 1 / zeta2_inv_vb)

    }


  } else {

    mu_alpha_vb <- NULL
    sig2_alpha_vb <- NULL

  }

  if (!is.null(r)) {

    mu_c0_vb <- rnorm(p)
    mu_c_vb <- matrix(rnorm(r * d), nrow = r)

  } else {

    mu_c0_vb <- mu_c_vb <- NULL

  }


  d_init <- d
  G_init <- G
  p_init <- p
  q_init <- q
  r_init <- r
  ind_bin_init <- ind_bin

  link_init <- link

  list_init <- create_named_list_(d_init, G_init, p_init, q_init, r_init,
                                  link_init, ind_bin_init, gam_vb, mu_beta_vb,
                                  sig2_inv_vb, sig2_beta_vb, tau_vb, mu_alpha_vb,
                                  sig2_alpha_vb, mu_c0_vb, mu_c_vb)

  class(list_init) <- "out_init"

  list_init
}
