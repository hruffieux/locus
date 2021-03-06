% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_hyper_init.R
\name{set_hyper}
\alias{set_hyper}
\title{Gather model hyperparameters provided by the user.}
\usage{
set_hyper(
  d,
  p,
  lambda,
  nu,
  a,
  b,
  eta,
  kappa,
  link = "identity",
  ind_bin = NULL,
  q = NULL,
  phi = NULL,
  xi = NULL,
  m0 = NULL,
  s02 = NULL,
  G = NULL,
  struct = FALSE
)
}
\arguments{
\item{d}{Number of responses.}

\item{p}{Number of candidate predictors.}

\item{lambda}{Vector of length 1 providing the values of hyperparameter
\eqn{\lambda} for the prior distribution of \eqn{\sigma^{-2}}. \eqn{\sigma}
represents the typical size of nonzero effects.}

\item{nu}{Vector of length 1 providing the values of hyperparameter \eqn{\nu}
for the prior distribution of \eqn{\sigma^{-2}}. \eqn{\sigma} represents
the typical size of nonzero effects.}

\item{a}{Vector of length 1 or p providing the values of hyperparameter
\eqn{a} for the prior distributions for the proportion of responses
associated with each candidate predictor, \eqn{\omega} (vector of length p).
If of length 1, the provided value is repeated p times.}

\item{b}{Vector of length 1 or p providing the values of hyperparameter
\eqn{b} for the prior distributions for the proportion of responses
associated with each candidate predictor, \eqn{\omega} (vector of length p).
If of length 1, the provided value is repeated p times.}

\item{eta}{Vector of length 1 or d for \code{link = "identity"}, and of
length 1 or d_cont = d - length(ind_bin) (the number of continuous response
variables) for \code{link = "mix"}. Provides the values of
hyperparameter \eqn{\eta} for the prior distributions of the continuous
response residual precisions, \eqn{\tau}. If of length 1, the provided
value is repeated d, resp. d_cont, times. Must be \code{NULL} for
\code{link = "logit"} and \code{link = "probit"}.}

\item{kappa}{Vector of length 1 or d for \code{link = "identity"}, and of
length 1 or d_cont = d - length(ind_bin) (the number of continuous response
variables) for \code{link = "mix"}. Provides the values of hyperparameter
\eqn{\kappa} for the prior distributions of the response residual
precisions, \eqn{\tau}. If of length 1, the provided value is repeated d,
resp. d_cont, times. Must be \code{NULL} for \code{link = "logit"} and
\code{link = "probit"}.}

\item{link}{Response link. Must be "\code{identity}" for linear regression,
"\code{logit}" for logistic regression, "\code{probit}"
for probit regression, or "\code{mix}" for a mix of identity and probit
link functions (in this case, the indices of the binary responses must be
gathered in argument \code{ind_bin}, see below).}

\item{ind_bin}{If \code{link = "mix"}, vector of indices corresponding to
the binary variables in \code{Y}. Must be \code{NULL} if
\code{link != "mix"}.}

\item{q}{Number of covariates. Default is \code{NULL}, for \code{Z}
\code{NULL}.}

\item{phi}{Vector of length 1 or q providing the values of hyperparameter
\eqn{\phi} for the prior distributions for the sizes of the nonzero
covariate effects, \eqn{\zeta}. If of length 1, the provided value is
repeated q times. Default is \code{NULL}, for \code{Z} \code{NULL}.}

\item{xi}{Vector of length 1 or q providing the values of hyperparameter
\eqn{\xi} for the prior distributions for the sizes of the nonzero
covariate effects, \eqn{\zeta}. If of length 1, the provided value is
repeated q times. Default is \code{NULL}, for \code{Z} \code{NULL}.}

\item{m0}{Vector of length 1 or p. Hyperparameter when \code{list_struct} 
non-\code{NULL}. Default is \code{NULL}.}

\item{s02}{Variance hyperparameter when \code{list_struct} is 
non-\code{NULL}. Default is \code{NULL}.}

\item{G}{Number of candidate predictor groups when using the group selection
model from the \code{\link{locus}} function. Default is \code{NULL},
for no group selection.}

\item{struct}{Boolean indicating the use of structured sparse priors
set through the \code{\link{set_struct}} function. Default is \code{FALSE},
for no structured selection.}
}
\value{
An object of class "\code{hyper}" preparing user hyperparameter in a
  form that can be passed to the \code{\link{locus}} function.
}
\description{
This function must be used to provide hyperparameter values for the model
used in \code{\link{locus}}.
}
\details{
The \code{\link{locus}} function can also be used with default hyperparameter
choices (without using \code{\link{set_hyper}}) by setting its argument
\code{list_hyper} to \code{NULL}.
}
\examples{
seed <- 123; set.seed(seed)

###################
## Simulate data ##
###################

## Examples using small problem sizes:
##
n <- 200; p <- 200; p0 <- 20; d <- 20; d0 <- 15; q <- 2

## Candidate predictors (subject to selection)
##
# Here we simulate common genetic variants (but any type of candidate
# predictors can be supplied).
# 0 = homozygous, major allele, 1 = heterozygous, 2 = homozygous, minor allele

X_act <- matrix(rbinom(n * p0, size = 2, p = 0.25), nrow = n)
X_inact <- matrix(rbinom(n * (p - p0), size = 2, p = 0.25), nrow = n)

shuff_x_ind <- sample(p)
X <- cbind(X_act, X_inact)[, shuff_x_ind]

bool_x_act <- shuff_x_ind <= p0

pat_act <- beta <- matrix(0, nrow = p0, ncol = d0)
pat_act[sample(p0*d0, floor(p0*d0/5))] <- 1
beta[as.logical(pat_act)] <-  rnorm(sum(pat_act))

## Covariates (not subject to selection)
##
Z <- matrix(rnorm(n * q), nrow = n)

alpha <-  matrix(rnorm(q * d), nrow = q)

## Gaussian responses
##
Y_act <- matrix(rnorm(n * d0, mean = X_act \%*\% beta, sd = 0.5), nrow = n)
Y_inact <- matrix(rnorm(n * (d - d0), sd = 0.5), nrow = n)
shuff_y_ind <- sample(d)
Y <- cbind(Y_act, Y_inact)[, shuff_y_ind] + Z \%*\% alpha

## Binary responses
##
Y_bin <- ifelse(Y > 0, 1, 0)

########################
## Infer associations ##
########################

## Continuous responses
##

# No covariate
#
# a and b chosen so that the prior mean number of responses associated with
# each candidate predictor is 1/4.
list_hyper_g <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
                          eta = 1, kappa = apply(Y, 2, var),
                          link = "identity")

# We take p0_av = p0 (known here); this choice may result in variable
# selections that are (too) conservative in some cases. In practice, it is
# advised to set p0_av as a slightly overestimated guess of p0, or perform
# cross-validation using function `set_cv'.

vb_g <- locus(Y = Y, X = X, p0_av = p0, link = "identity",
              list_hyper = list_hyper_g, user_seed = seed)

# With covariates
#
list_hyper_g_z <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
                            eta = 1, kappa = apply(Y, 2, var),
                            link = "identity", q = q, phi = 1, xi = 1)

vb_g_z <- locus(Y = Y, X = X, p0_av = p0, Z = Z, link = "identity",
                list_hyper = list_hyper_g_z, user_seed = seed)

## Binary responses
##
list_hyper_logit <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
                              eta = NULL, kappa = NULL, link = "logit",
                              q = q, phi = 1, xi = 1)

vb_logit <- locus(Y = Y_bin, X = X, p0_av = p0, Z = Z, link = "logit",
                  list_hyper = list_hyper_logit, user_seed = seed)

list_hyper_probit <- set_hyper(d, p, lambda = 1, nu = 1, a = 1, b = 4*d-1,
                               eta = NULL, kappa = NULL, link = "probit",
                               q = q, phi = 1, xi = 1)

vb_probit <- locus(Y = Y_bin, X = X, p0_av = p0, Z = Z, link = "probit",
                   list_hyper = list_hyper_probit, user_seed = seed)


## Mix of continuous and binary responses
##
Y_mix <- cbind(Y, Y_bin)
ind_bin <- (d+1):(2*d)

list_hyper_mix <- set_hyper(2*d, p, lambda = 1, nu = 1, a = 1, b = 8*d-1,
                            eta = 1, kappa = apply(Y, 2, var), link = "mix",
                            ind_bin = ind_bin, q = q, phi = 1, xi = 1)

vb_mix <- locus(Y = Y_mix, X = X, p0_av = p0, Z = Z, link = "mix",
                ind_bin = ind_bin, list_hyper = list_hyper_mix,
                user_seed = seed)

}
\seealso{
\code{\link{set_init}}, \code{\link{locus}}
}
