rm(list = ls())

set.seed(123)

############################
## simulate basic dataset ##
############################

n <- 100; p <- 75; d <- 20; q <- 3; p0 <- 10

# covariates (not subject to selection)
Z <- matrix(rnorm(n * q), nrow = n)

alpha <-  matrix(rnorm(q * d), nrow = q)

# candidate predictors (subject to selection)
X_act <- matrix(rbinom(n * p0, size = 2, p = 0.2), nrow = n)
X_inact <- matrix(rbinom(n * (p - p0), size = 2, p = 0.2), nrow = n)
X <- cbind(X_act, X_inact)[, sample(p)]

beta <-  matrix(rnorm(p0 * d), nrow = p0)

# Gaussian outcomes
Y <- matrix(rnorm(n * d, mean = Z %*% alpha + X_act %*% beta, sd = 1), nrow = n)

# Binary outcomes
Y_bin <- ifelse(Y > 0, 1, 0)

# remove constant variables (needed for checking dimension consistency)
X <- scale(X); Z <- scale(Z)
rm_cst <- function(mat_sc) mat_sc[, !is.nan(colSums(mat_sc))]
rm_coll <- function(mat_sc) mat_sc[, !duplicated(mat_sc, MARGIN = 2)]

X <- rm_cst(X); Z <- rm_cst(Z)
X <- rm_coll(X); Z <- rm_coll(Z)

p <- ncol(X); q <- ncol(Z)

####################
## locus settings ##
####################

# hyperparameter (prior number of active predictors)
p0_av <- p0


#####################
## locus inference ##
#####################

# Gaussian outcomes, no covariates
vb <- locus(Y = Y, X = X, p0_av = p0_av, family = "gaussian")

# Gaussian outcomes, with covariates
vb_z <- locus(Y = Y, X = X, p0_av = p0_av, Z = Z, family = "gaussian")

# Binomial outcomes, no covariates
vb_bin <- locus(Y = Y_bin, X = X, p0_av = p0_av, family = "binomial")

# Binomial outcomes, with covariates
vb_bin_z <- locus(Y = Y_bin, X = X, p0_av = p0_av, Z = Z, family = "binomial")

