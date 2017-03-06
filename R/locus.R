#' Fit sparse multivariate regression models using variational inference.
#'
#' Variational approximation procedure fitting sparse multivariate regression
#' models for combined selection of predictors and associated responses in
#' high-dimensional set-ups. Dependence across responses linked to the same
#' predictors is captured through the model hierarchical structure.
#' The responses can be purely continuous, purely binary (logit or probit link
#' fits), or a mix of continuous and binary variables.
#'
#' The continuous response variables in \code{Y} (if any) will be centered
#' before application of the variational algorithm, and the candidate predictors
#' and covariates resp. in \code{X} and \code{Z} will be standardized. An
#' intercept will be added if \code{link} is \code{"logit"}, \code{"probit"} or
#' \code{"mix"} (do not supply it in \code{X} or \code{Z}).
#'
#' @param Y Response data matrix of dimension n x d, where n is the number of
#'   observations and d is the number of response variables.
#' @param X Input matrix of dimension n x p, where p is the number of candidate
#'   predictors. \code{X} cannot contain NAs. No intercept must be supplied.
#' @param p0_av Prior average number of predictors expected to be included in
#'   the model. Must be \code{NULL} if \code{list_init} and \code{list_hyper}
#'   are both non-\code{NULL} or if \code{list_cv} is non-\code{NULL}. Can also
#'   be a vector of length p with entry s corresponding to the prior probability
#'   that candidate predictor s is associated with at least one response.
#' @param Z Covariate matrix of dimension n x q, where q is the number of
#'   covariates. Variables in \code{Z} are not subject to selection. \code{NULL}
#'   if no covariate. Factor covariates must be supplied after transformation to
#'   dummy coding. No intercept must be supplied.
#' @param V Covariates matrix of dimension p x r, where r is the number of
#'   variables representing external information on the candidate predictors
#'   which may make their selection more or less likely. \code{NULL} if no such
#'   information.
#' @param link Response link. Must be "\code{identity}" for linear regression,
#'   "\code{logit}" for logistic regression, "\code{probit}" for probit
#'   regression, or "\code{mix}" for a mix of identity and probit link functions
#'   (in this case, the indices of the binary responses must be gathered in
#'   argument \code{ind_bin}, see below).
#' @param ind_bin If \code{link = "mix"}, vector of indices corresponding to
#'   the binary variables in \code{Y}. Must be \code{NULL} if
#'   \code{link != "mix"}.
#' @param list_hyper An object of class "\code{hyper}" containing the model
#'   hyperparameters. Must be filled using the \code{\link{set_hyper}}
#'   function or must be \code{NULL} for default hyperparameters.
#' @param list_init An object of class "\code{init}" containing the initial
#'   variational parameters. Must be filled using the \code{\link{set_init}}
#'   function or be \code{NULL} for a default initialization.
#' @param list_cv An object of class "\code{cv}" containing settings for
#'   choosing the prior average number of predictors expected to be included in
#'   the model, \code{p0_av}, by cross-validation. Must be filled using the
#'   \code{\link{set_cv}} function or must be \code{NULL} for no
#'   cross-validation. If non-\code{NULL}, \code{p0_av}, \code{list_init} and
#'   \code{list_hyper} must all be \code{NULL}. Cross-validation only available
#'   for \code{link = "identity"}.
#' @param list_blocks An object of class "\code{blocks}" containing settings for
#'   parallel inference on a partitioned predictor space. Must be filled using
#'   the \code{\link{set_blocks}} function or must be \code{NULL} for no
#'   partitioning.
#' @param user_seed Seed set for reproducible default choices of hyperparameters
#'   (if \code{list_hyper} is \code{NULL}) and initial variational parameters
#'   (if \code{list_init} is \code{NULL}). Also used at the cross-validation
#'   stage (if \code{list_cv} is non-\code{NULL}). Default is \code{NULL}, no
#'   seed set.
#' @param tol Tolerance for the stopping criterion.
#' @param maxit Maximum number of iterations allowed.
#' @param batch If \code{TRUE} a fast batch updating scheme is used
#'   (recommended).
#' @param save_hyper If \code{TRUE}, the hyperparameters used for the model are
#'   saved as output.
#' @param save_init If \code{TRUE}, the initial variational parameters used for
#'   the inference are saved as output.
#' @param verbose If \code{TRUE}, messages are displayed during execution.
#'
#' @return An object of class "\code{vb}" containing the following variational
#'   estimates and settings:
#'  \item{lb_opt}{Optimized variational lower bound for the marginal
#'                log-likelihood.}
#'  \item{gam_vb}{Posterior inclusion probability matrix of dimension p x d.
#'                Entry (s, t) corresponds to the posterior probability of
#'                association between candidate predictor s and response t.}
#'  \item{mu_alpha_vb}{Matrix of dimension q x d whose entries are the posterior
#'                     mean regression coefficients for the covariates provided
#'                     in \code{Z} (if \code{link = "logit"},
#'                     \code{link = "logit"} or
#'                     \code{link = "mix"} also for the intercept).
#'                     \code{NULL} if \code{Z} is \code{NULL}.}
#'  \item{om_vb}{Vector of length p containing the posterior mean of omega.
#'               Entry s controls the proportion of responses associated with
#'               candidate predictor s.}
#'  \item{rmvd_cst_x, rmvd_cst_z}{Vectors containing the indices of constant
#'                                variables in \code{X} (resp. \code{Z}) removed
#'                                prior to the analysis.}
#'  \item{rmvd_coll_x, rmvd_coll_z}{Vectors containing the indices of variables
#'                                  in \code{X} (resp. \code{Z}) removed prior
#'                                  to the analysis because collinear to other
#'                                  variables. The entry name indicates the
#'                                  corresponding variable kept in the analysis
#'                                  (i.e., that causing the collinearity for the
#'                                  entry in question).}
#'  \item{list_hyper, list_init}{If \code{save_hyper}, resp. \code{save_init},
#'                               \code{TRUE}, hyperparameters, resp. initial
#'                               variational parameters, used for inference are
#'                               saved as output.}
#' @examples
#' user_seed <- 123; set.seed(user_seed)
#' n <- 200; p <- 300; p0 <- 50; d <- 40; d0 <- 30
#' list_X <- generate_snps(n = n, p = p)
#' list_Y <- generate_phenos(n = n, d = d, var_err = 1)
#'
#' # Continuous outcomes
#' #
#' dat_g <- generate_dependence(list_snps = list_X, list_phenos = list_Y,
#'                              ind_d0 = sample(1:d, d0),
#'                              ind_p0 = sample(1:p, p0),
#'                              vec_prob_sh = 0.1, family = "gaussian",
#'                              max_tot_pve = 0.9)
#'
#' # we take p0_av = p0 (known here); this choice may result in variable
#' # selections that are (too) conservative in some cases. In practice, often
#' # p0_av as a slightly overestimated guess of p0.
#' vb_g <- locus(Y = dat_g$phenos, X = dat_g$snps, p0_av = p0,
#'               link = "identity", user_seed = user_seed)
#'
#' # Continuous outcomes with covariates
#' #
#' q <- 4
#' Z <- matrix(rnorm(n * q), nrow = n)
#' vb_g_z <- locus(Y = dat_g$phenos, X = dat_g$snps, p0_av = p0,  Z = Z,
#'                 link = "identity", user_seed = user_seed)
#'
#' # Continuous outcomes with external annotation
#' #
#' r <- 4
#' V <- matrix(rnorm(p * r), nrow = p)
#' bool_p0 <- rowSums(dat_g$pat) > 0
#' V[bool_p0, ] <- rnorm(sum(bool_p0) * r, mean = 2) # informative annotations
#' vb_g_v <- locus(Y = dat_g$phenos, X = dat_g$snps, p0_av = p0,  V = V,
#'                 link = "identity", user_seed = user_seed)
#'
#' # Binary outcomes
#' #
#' dat_b <- generate_dependence(list_snps = list_X, list_phenos = list_Y,
#'                              ind_d0 = sample(1:d, d0),
#'                              ind_p0 = sample(1:p, p0),
#'                              vec_prob_sh = 0.1, family = "binomial",
#'                              max_tot_pve = 0.9)
#'
#' vb_logit <- locus(Y = dat_b$phenos, X = dat_b$snps, p0_av = p0,
#'                   link = "logit", user_seed = user_seed)
#'
#' vb_probit <- locus(Y = dat_b$phenos, X = dat_b$snps, p0_av = p0,
#'                    link = "probit", user_seed = user_seed)
#'
#' # Binary outcomes with covariates
#' #
#' vb_logit_z <- locus(Y = dat_b$phenos, X = dat_b$snps, p0_av = p0,  Z = Z,
#'                     link = "logit", user_seed = user_seed)
#'
#' vb_probit_z <- locus(Y = dat_b$phenos, X = dat_b$snps, p0_av = p0,  Z = Z,
#'                      link = "probit", user_seed = user_seed)
#'
#' # Mix of continuous and binary outcomes
#' #
#' Y_mix <- cbind(dat_g$phenos, dat_b$phenos)
#' ind_bin <- (d+1):(2*d)
#' p0_mix <- sum(rowSums(cbind(dat_g$pat, dat_b$pat)) > 0)
#'
#' vb_mix <- locus(Y = Y_mix, X = dat_b$snps, p0_av = p0, link = "mix",
#'                 ind_bin = ind_bin, user_seed = user_seed)
#'
#' # Mix of continuous and binary outcomes with covariates
#' #
#' vb_mix_z <- locus(Y = Y_mix, X = dat_b$snps, p0_av = p0,  Z = Z,
#'                   link = "mix", ind_bin = ind_bin, user_seed = user_seed)
#'
#' @seealso \code{\link{set_hyper}}, \code{\link{set_init}},
#'   \code{\link{set_cv}}, \code{\link{set_blocks}}
#'
#' @export
#'
locus <- function(Y, X, p0_av, Z = NULL, V = NULL, link = "identity",
                  ind_bin = NULL, list_hyper = NULL, list_init = NULL,
                  list_cv = NULL, list_blocks = NULL, user_seed = NULL,
                  tol = 1e-4, maxit = 1000, batch = TRUE, save_hyper = FALSE,
                  save_init = FALSE, verbose = TRUE) { ##

  if (verbose) cat("== Preparing the data ... \n")
  dat <- prepare_data_(Y, X, Z, V, link, ind_bin, user_seed, tol, maxit, batch,
                       verbose)

  bool_rmvd_x <- dat$bool_rmvd_x
  bool_rmvd_z <- dat$bool_rmvd_z
  bool_rmvd_v <- dat$bool_rmvd_v

  X <- dat$X
  Y <- dat$Y
  Z <- dat$Z
  V <- dat$V

  n <- nrow(X)
  p <- ncol(X)
  d <- ncol(Y)
  r <- ncol(V)

  names_x <- colnames(X)
  names_y <- colnames(Y)

  if (!is.null(Z)) {
    q <- ncol(Z)
    names_z <- colnames(Z)
  } else {
    q <- NULL
    names_z <- NULL
  }

  if (!is.null(V)) {
    r <- ncol(V)
    names_v <- colnames(V)
  } else {
    r <- NULL
    names_v <- NULL
  }

  if (verbose) cat("... done. == \n\n")



  if (!is.null(list_cv) & is.null(list_blocks)) { ## TODO: allow cross-validation when list_blocks is used.

      if (verbose) {
        cat("=============================== \n")
        cat("===== Cross-validation... ===== \n")
        cat("=============================== \n")
      }
      list_cv <- prepare_cv_(list_cv, n, p, r, bool_rmvd_x, p0_av, link,
                             list_hyper, list_init, verbose)

      p_star <- cross_validate_(Y, X, Z, list_cv, user_seed, verbose)

  } else {

    if (!is.null(list_blocks)) {

      list_blocks <- prepare_blocks_(list_blocks, r, bool_rmvd_x, list_cv)

      n_bl <- list_blocks$n_bl
      n_cpus <- list_blocks$n_cpus
      vec_fac_bl <- list_blocks$vec_fac_bl

    }

    if (is.null(list_hyper) | is.null(list_init)) {

      p_star <- convert_p0_av_(p0_av, p, verbose)

      # remove the entries corresponding to the removed constant covariates in X (if any)
      if (length(p_star) > 1) p_star <- p_star[!bool_rmvd_x]

    } else {

      if (!is.null(p0_av))
        warning(paste("Provided argument p0_av not used, as both list_hyper ",
                      "and list_init were provided.", sep = ""))

      p_star <- NULL

    }

  }

  if (verbose) cat("== Preparing the hyperparameters ... \n\n")
  list_hyper <- prepare_list_hyper_(list_hyper, Y, p, p_star, q, r, link, ind_bin,
                                    bool_rmvd_x, bool_rmvd_z, bool_rmvd_v,
                                    names_x, names_y, names_z, verbose)
  if (verbose) cat("... done. == \n\n")

  if (verbose) cat("== Preparing the parameter initialization ... \n\n")
  list_init <- prepare_list_init_(list_init, Y, p, p_star, q, r, link, ind_bin,
                                  bool_rmvd_x, bool_rmvd_z, bool_rmvd_v,
                                  user_seed, verbose)
  if (verbose) cat("... done. == \n\n")


  if (link != "identity") { # adds an intercept for logistic/probit regression

    if (is.null(q)) {

      Z <- matrix(1, nrow = n, ncol = 1)

      # uninformative prior
      list_hyper$phi <- list_hyper$xi <- 1e-3

      list_init$mu_alpha_vb <- matrix(0, nrow = 1, ncol = d)

      if (link == "probit") {

        list_init$sig2_alpha_vb <- 1

      } else {

        list_init$sig2_alpha_vb <- matrix(1, nrow = 1, ncol = d)

      }

    } else{

      Z <- cbind(rep(1, n), Z)

      # uninformative prior
      list_hyper$phi <- c(1e-3, list_hyper$phi)
      list_hyper$xi <- c(1e-3, list_hyper$xi)

      list_init$mu_alpha_vb <- rbind(rep(0, d), list_init$mu_alpha_vb)

      if (link == "probit") {

        list_init$sig2_alpha_vb <- c(1, list_init$sig2_alpha_vb)

      } else {

        list_init$sig2_alpha_vb <- rbind(rep(1, d), list_init$sig2_alpha_vb)

      }


    }
    colnames(Z)[1] <- "Intercept"
  }


  if (verbose){
    cat(paste("============================================================== \n",
              "== Variational inference for sparse multivariate regression == \n",
              "============================================================== \n\n",
              sep = ""))
  }


  if (is.null(list_blocks)) {

    if (link == "identity") {

      if (is.null(q)) {

        if (is.null(r)) {
          vb <- locus_core_(Y, X, list_hyper, list_init$gam_vb,
                            list_init$mu_beta_vb, list_init$sig2_beta_vb,
                            list_init$tau_vb, tol, maxit, batch, verbose)
        } else {
          vb <- locus_info_core_(Y, X, V, list_hyper, list_init$gam_vb,
                                list_init$mu_beta_vb, list_init$mu_c0_vb,
                                list_init$mu_c_vb, list_init$sig2_beta_vb,
                                list_init$tau_vb, tol, maxit, batch, verbose)
        }

      } else {
        vb <- locus_z_core_(Y, X, Z, list_hyper, list_init$gam_vb,
                            list_init$mu_alpha_vb, list_init$mu_beta_vb,
                            list_init$sig2_alpha_vb, list_init$sig2_beta_vb,
                            list_init$tau_vb, tol, maxit, batch, verbose)
      }

    } else if (link == "logit"){

      vb <- locus_logit_core_(Y, X, Z, list_hyper, list_init$chi_vb,
                              list_init$gam_vb, list_init$mu_alpha_vb,
                              list_init$mu_beta_vb, list_init$sig2_alpha_vb,
                              list_init$sig2_beta_vb, tol, maxit, batch, verbose)

    } else if (link == "probit"){

      vb <- locus_probit_core_(Y, X, Z, list_hyper, list_init$gam_vb,
                               list_init$mu_alpha_vb, list_init$mu_beta_vb,
                               list_init$sig2_alpha_vb, list_init$sig2_beta_vb,
                               tol, maxit, batch, verbose)

    } else {

      vb <- locus_mix_core_(Y, X, Z, ind_bin, list_hyper, list_init$gam_vb,
                            list_init$mu_alpha_vb, list_init$mu_beta_vb,
                            list_init$sig2_alpha_vb, list_init$sig2_beta_vb,
                            list_init$tau_vb, tol, maxit, batch, verbose)
    }

  } else {


    list_pos_bl <- split(1:p, vec_fac_bl)

    split_bl_hyper <- lapply(list_pos_bl, function(pos_bl) {
      list_hyper$p_hyper <- length(n_bl)
      list_hyper$a <- list_hyper$a[pos_bl]
      list_hyper$b <- list_hyper$b[pos_bl]
      list_hyper
    })

    split_bl_init <- lapply(list_pos_bl, function(pos_bl) {
      list_init$p_init <- length(pos_bl)
      list_init$gam_vb <- list_init$gam_vb[pos_bl,, drop = FALSE]
      list_init$mu_beta_vb <- list_init$mu_beta_vb[pos_bl,, drop = FALSE]
      if (link == "logit")
        list_init$sig2_beta_vb <- list_init$sig2_beta_vb[pos_bl,, drop = FALSE]

      list_init
    })

    locus_bl_ <- function(k) {

      X_bl <- X[, list_pos_bl[[k]], drop = FALSE]

      list_hyper_bl <- split_bl_hyper[[k]]
      list_init_bl <- split_bl_init[[k]]

      if (link == "identity") {
        if (is.null(q))
          vb_bl <- locus_core_(Y, X_bl, list_hyper_bl,
                               list_init_bl$gam_vb, list_init_bl$mu_beta_vb,
                               list_init_bl$sig2_beta_vb, list_init_bl$tau_vb,
                               tol, maxit, batch, verbose = FALSE)
        else
          vb_bl <- locus_z_core_(Y, X_bl, Z, list_hyper_bl, list_init_bl$gam_vb,
                                 list_init_bl$mu_alpha_vb,list_init_bl$mu_beta_vb,
                                 list_init_bl$sig2_alpha_vb,
                                 list_init_bl$sig2_beta_vb, list_init_bl$tau_vb,
                                 tol, maxit, batch, verbose = FALSE)

      } else if (link == "logit") {

        vb_bl <- locus_logit_core_(Y, X_bl, Z, list_hyper_bl,
                                   list_init_bl$chi_vb, list_init_bl$gam_vb,
                                   list_init_bl$mu_alpha_vb, list_init_bl$mu_beta_vb,
                                   list_init_bl$sig2_alpha_vb,
                                   list_init_bl$sig2_beta_vb, tol, maxit, batch,
                                   verbose = FALSE)

      } else  if (link == "probit") {

        vb_bl <- locus_probit_core_(Y, X_bl, Z, list_hyper_bl,
                                    list_init_bl$gam_vb, list_init_bl$mu_alpha_vb,
                                    list_init_bl$mu_beta_vb,
                                    list_init_bl$sig2_alpha_vb,
                                    list_init_bl$sig2_beta_vb, tol, maxit, batch,
                                    verbose = FALSE)

      } else {

        vb_bl <- locus_mix_core_(Y, X_bl, Z, ind_bin, list_hyper_bl,
                                 list_init_bl$gam_vb, list_init_bl$mu_alpha_vb,
                                 list_init_bl$mu_beta_vb,
                                 list_init_bl$sig2_alpha_vb,
                                 list_init_bl$sig2_beta_vb, list_init_bl$tau_vb,
                                 tol, maxit, batch, verbose = FALSE)
      }

      vb_bl

    }

    list_vb <- parallel::mclapply(1:n_bl, function(k) locus_bl_(k), mc.cores = n_cpus)

    names_vec <- c("lb_opt", "om_vb", "p_star")
    names_mat <- "gam_vb"

    vb <- c(lapply(names_vec, function(key) do.call(c, lapply(list_vb, `[[`, key))),
            lapply(names_mat, function(key) do.call(rbind, lapply(list_vb, `[[`, key))))

    names(vb) <- c(names_vec, names_mat)

  }

  vb$p_star <- p_star

  vb$rmvd_cst_x <- dat$rmvd_cst_x
  vb$rmvd_coll_x <- dat$rmvd_coll_x
  if (!is.null(Z)) {
    vb$rmvd_cst_z <- dat$rmvd_cst_z
    vb$rmvd_coll_z <- dat$rmvd_coll_z
  }
  if (!is.null(V)) {
    vb$rmvd_cst_v <- dat$rmvd_cst_v
    vb$rmvd_coll_v <- dat$rmvd_coll_v
  }

  if (save_hyper) vb$list_hyper <- list_hyper
  if (save_init) vb$list_init <- list_init

  class(vb) <- "locus"


  if (verbose) cat("... done. == \n\n")

  vb

}
