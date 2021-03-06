# This file is part of the `locus` R package:
#     https://github.com/hruffieux/locus
#

#' Fit sparse multivariate regression models using variational inference.
#'
#' Variational approximation procedure fitting sparse multivariate regression
#' models for combined selection of predictors and associated responses in
#' high-dimensional set-ups. Dependence across responses linked to the same
#' predictors is modelled through the model hierarchical structure.
#' The responses can be purely continuous, purely binary (logit or probit link
#' fits), or a mix of continuous and binary variables.
#'
#'
#' The optimization uses efficient block coordinate ascent schemes, for
#' which convergence is ensured as the objective (elbo) is multiconcave
#' for the selected blocks, i.e., it is concave in each block of parameters
#' whose updates are made simultaneously, see Wu et al. (reference Section
#' below).
#'
#' The continuous response variables in \code{Y} (if any) will be centered
#' before application of the variational algorithm, and the candidate predictors
#' and covariates resp. in \code{X} and \code{Z} will be standardized. An
#' intercept will be added if \code{link} is \code{"logit"}, \code{"probit"} or
#' \code{"mix"} (do not supply it in \code{X} or \code{Z}).
#'
#'
#' @param Y Response data matrix of dimension n x d, where n is the number of
#'   samples and d is the number of response variables.
#' @param X Input matrix of dimension n x p, where p is the number of candidate
#'   predictors. \code{X} cannot contain NAs. No intercept must be supplied.
#' @param p0_av Prior average number of predictors (or groups of predictors if 
#'   \code{list_groups} is non-\code{NULL}) expected to be included in the 
#'   model.  Can also be a vector of length p (resp. of length the number of 
#'   groups) with entry s corresponding to the prior probability that candidate 
#'   predictor s (resp. group s) is associated with at least one response. Must 
#'   be \code{NULL} if \code{list_init} and \code{list_hyper} are both 
#'   non-\code{NULL} or if \code{list_cv} is non-\code{NULL}.
#' @param Z Covariate matrix of dimension n x q, where q is the number of
#'   covariates. Variables in \code{Z} are not subject to selection. \code{NULL}
#'   if no covariate. Factor covariates must be supplied after transformation to
#'   dummy coding. No intercept must be supplied.
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
#' @param list_groups An object of class "\code{groups}" containing settings for
#'   group selection of candidate predictors. Must be filled using the
#'   \code{\link{set_groups}} function or must be \code{NULL} for group
#'   selection.
#' @param list_struct An object of class "\code{struct}" containing settings for
#'   structure sparsity priors. Must be filled using the
#'   \code{\link{set_struct}} function or must be \code{NULL} for structured
#'   selection.
#' @param user_seed Seed set for reproducible default choices of hyperparameters
#'   (if \code{list_hyper} is \code{NULL}) and initial variational parameters
#'   (if \code{list_init} is \code{NULL}). Also used at the cross-validation
#'   stage (if \code{list_cv} is non-\code{NULL}). Default is \code{NULL}, no
#'   seed set.
#' @param tol Tolerance for the stopping criterion.
#' @param maxit Maximum number of iterations allowed.
#' @param anneal Parameters for annealing scheme. Must be a vector whose first
#'   entry is sets the type of ladder: 1 = geometric spacing, 2 = harmonic
#'   spacing or 3 = linear spacing, the second entry is the initial temperature,
#'   and the third entry is the ladder size. If \code{NULL} (default), no
#'   annealing is performed.
#' @param save_hyper If \code{TRUE}, the hyperparameters used for the model are
#'   returned.
#' @param save_init If \code{TRUE}, the initial variational parameters used for
#'   the inference are returned (note that the size of the resulting objects is
#'   likely to be large). Default is \code{FALSE}.
#' @param full_output If \code{TRUE}, the inferred variational parameters for 
#'   all parameters are returned.
#' @param verbose If \code{TRUE}, messages are displayed during execution.
#' @param checkpoint_path Path where to save temporary checkpoint outputs. 
#'   Default is \code{NULL}, for no checkpointing.
#'
#' @return An object of class "\code{vb}" containing the following variational
#'   estimates and settings:
#'  \item{gam_vb}{Posterior inclusion probability matrix of dimension p x d.
#'                Entry (s, t) corresponds to the posterior probability of
#'                association between candidate predictor s and response t.}
#'  \item{alpha_vb}{Matrix of dimension q x d whose entries are the posterior
#'                     mean regression coefficients for the covariates provided
#'                     in \code{Z} (if \code{link = "logit"},
#'                     \code{link = "logit"} or
#'                     \code{link = "mix"} also for the intercept).
#'                     \code{NULL} if \code{Z} is \code{NULL}.}
#'  \item{om_vb}{Vector of length p containing the posterior mean of omega.
#'               Entry s controls the proportion of responses associated with
#'               candidate predictor s.}
#'  \item{converged}{A boolean indicating whether the algorithm has converged
#'                   before reaching \code{maxit} iterations.}
#'  \item{it}{Final number of iterations.}
#'  \item{lb_opt}{Optimized variational lower bound for the marginal
#'                log-likelihood (ELBO).}
#'  \item{diff_lb}{Difference in ELBO between the last and penultimate
#'                 iterations. This may be a useful diagnostic information when
#'                 convergence has not been reached before \code{maxit}.}
#'  \item{p_star}{Vector of length 1 or p defining the applied sparsity control.}
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
#'  \item{group_labels}{If \code{list_groups} is non-\code{NULL}, labels of the
#'                      groups to which the candidate predictor belong (these
#'                      labels are gathered after removal of constant and
#'                      collinear predictors, whose indices are stored in
#'                      \code{rmvd_cst_x} and \code{rmvd_coll_x}).}
#'  \item{...}{Other specific outputs are possible depending on the model used.}
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
#' n <- 200; p <- 250; p0 <- 25; d <- 30; d0 <- 25; q <- 3
#'
#' ## Candidate predictors (subject to selection)
#' ##
#' # Here we simulate common genetic variants (but any type of candidate
#' # predictors can be supplied).
#' # 0 = homozygous, major allele, 1 = heterozygous, 2 = homozygous, minor allele
#' #
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
#'
#' ########################
#' ## Infer associations ##
#' ########################
#'
#' ## Continuous responses
#' ##
#' # We take p0_av = p0 (known here); this choice may, in some cases, result in
#' # (too) conservative variable selections. In practice, it is advised to set
#' # p0_av as a slightly overestimated guess of p0, or perform cross-validation
#' # using function `set_cv'.
#'
#' # No covariate
#' #
#' vb_g <- locus(Y = Y, X = X, p0_av = p0, link = "identity", user_seed = seed)
#'
#' # With covariates
#' #
#' vb_g_z <- locus(Y = Y, X = X, p0_av = p0,  Z = Z, link = "identity",
#'                 user_seed = seed)
#'
#' ## Binary responses
#' ##
#' vb_logit <- locus(Y = Y_bin, X = X, p0_av = p0, Z = Z, link = "logit",
#'                   user_seed = seed)
#'
#' vb_probit <- locus(Y = Y_bin, X = X, p0_av = p0, Z = Z, link = "probit",
#'                    user_seed = seed)
#'
#' ## Mix of continuous and binary responses
#' ##
#' Y_mix <- cbind(Y, Y_bin)
#' ind_bin <- (d+1):(2*d)
#'
#' vb_mix <- locus(Y = Y_mix, X = X, p0_av = p0, Z = Z, link = "mix",
#'                 ind_bin = ind_bin, user_seed = seed)
#'
#' @references
#' H. Ruffieux, A. C. Davison, J. Hager, I. Irincheeva. Efficient inference for
#'   genetic association studies with multiple outcomes. Biostatistics, 2017.
#'
#' Y. Xu, and W. Yin. A block coordinate descent method for
#'   regularized multiconvex optimization with applications to nonnegative
#'   tensor factorization and completion. SIAM Journal on imaging sciences, 6,
#'   pp.1758-1789, 2013.
#'
#' @seealso \code{\link{set_hyper}}, \code{\link{set_init}},
#'   \code{\link{set_cv}}, \code{\link{set_blocks}}, \code{\link{set_groups}}
#'   and \code{\link{set_struct}}.
#'
#' @export
#'
locus <- function(Y, X, p0_av, Z = NULL, link = "identity",
                  ind_bin = NULL, list_hyper = NULL, list_init = NULL, 
                  list_cv = NULL, list_blocks = NULL, list_groups = NULL, 
                  list_struct = NULL, user_seed = NULL, tol = 0.1, maxit = 1000, 
                  anneal = NULL, save_hyper = FALSE, save_init = FALSE, 
                  full_output = FALSE, verbose = TRUE, checkpoint_path = NULL) {
  
  check_structure_(verbose, "vector", "logical", 1)
  
  if (verbose) cat("== Preparing the data ... \n")
  
  check_annealing_(anneal, link, Z, list_groups, list_struct)
  
  dat <- prepare_data_(Y, X, Z, link, ind_bin, user_seed, tol, maxit, verbose, 
                       checkpoint_path)
  
  bool_rmvd_x <- dat$bool_rmvd_x
  bool_rmvd_z <- dat$bool_rmvd_z
  
  X <- dat$X
  Y <- dat$Y
  Z <- dat$Z
  
  n <- nrow(X)
  p <- ncol(X)
  d <- ncol(Y)
  
  names_x <- colnames(X)
  names_y <- colnames(Y)
  
  if (!is.null(Z)) {
    q <- ncol(Z)
    names_z <- colnames(Z)
  } else {
    q <- NULL
    names_z <- NULL
  }
  
  if (verbose) cat("... done. == \n\n")
  
  if (!is.null(list_cv) & is.null(list_blocks) & is.null(list_groups) & is.null(list_struct)) { ## TODO: allow cross-validation when list_blocks is used.
    
    if (verbose) {
      cat("=============================== \n")
      cat("===== Cross-validation... ===== \n")
      cat("=============================== \n")
    }
    list_cv <- prepare_cv_(list_cv, n, p, bool_rmvd_x, p0_av, link, list_hyper, 
                           list_init, verbose)
    
    p_star <- cross_validate_(Y, X, Z, link, ind_bin, list_cv, user_seed, verbose)
    
    vec_fac_gr <- vec_fac_st <- NULL
    
  } else {
    
    if (!is.null(list_blocks)) {
      
      list_blocks <- prepare_blocks_(list_blocks, bool_rmvd_x, list_cv, 
                                     list_groups, list_struct)
      
      n_bl <- list_blocks$n_bl
      n_cpus <- list_blocks$n_cpus
      vec_fac_bl <- list_blocks$vec_fac_bl
      
    }
    
    
    if (!is.null(list_groups)) {
      
      list_groups <- prepare_groups_(list_groups, X, q, bool_rmvd_x, link, list_cv)
      
      X <- list_groups$X
      vec_fac_gr <- list_groups$vec_fac_gr
      
    } else {
      
      vec_fac_gr <- NULL
      
    }
    
    
    if (!is.null(list_struct)) {
      
      list_struct <- prepare_struct_(list_struct, n, q, bool_rmvd_x, link, 
                                     list_cv, list_groups)
      
      vec_fac_st <- list_struct$vec_fac_st
      
    } else {
      
      vec_fac_st <- NULL
      
    }
    
    
    if (is.null(list_hyper) | is.null(list_init)) {
      
      if (is.null(list_groups)) p_tot <- p
      else p_tot <- length(unique(vec_fac_gr))
      
      p_star <- convert_p0_av_(p0_av, p_tot, list_blocks, verbose)
      
      # remove the entries corresponding to the removed constant covariates in X (if any)
      if (length(p_star) > 1) {
        if (is.null(list_groups)) p_star <- p_star[!bool_rmvd_x]
        else p_star <- p_star[unique(vec_fac_gr)]
      }
      
      
    } else {
      
      if (!is.null(p0_av))
        warning(paste0("Provided argument p0_av not used, as both list_hyper ",
                       "and list_init were provided."))
      
      p_star <- NULL
      
    }
    
  }
  
  if (verbose) cat("== Preparing the hyperparameters ... \n\n")
  
  list_hyper <- prepare_list_hyper_(list_hyper, Y, p, p_star, q, link, ind_bin,
                                    vec_fac_gr, vec_fac_st, bool_rmvd_x, 
                                    bool_rmvd_z, names_x, names_y, names_z, 
                                    verbose)
  
  if (verbose) cat("... done. == \n\n")
  
  if (verbose) cat("== Preparing the parameter initialization ... \n\n")
  
  list_init <- prepare_list_init_(list_init, Y, p, p_star, q, link, ind_bin, 
                                  vec_fac_gr, bool_rmvd_x, bool_rmvd_z,
                                  user_seed, verbose)
  
  if (verbose) cat("... done. == \n\n")
  
  nq <- is.null(q)
  
  if (link != "identity") { # adds an intercept for logistic/probit regression
    
    if (nq) {
      
      Z <- matrix(1, nrow = n, ncol = 1)
      
      # uninformative prior
      list_hyper$phi <- list_hyper$xi <- 1e-3
      
      list_init$alpha_vb <- matrix(0, nrow = 1, ncol = d)
      
      if (link == "probit") {
        
        list_init$sig2_alpha_vb <- 1
        
      } else {
        
        list_init$sig2_alpha_vb <- matrix(1, nrow = 1, ncol = d)
        
      }
      
    } else {
      
      Z <- cbind(rep(1, n), Z)
      
      # uninformative prior
      list_hyper$phi <- c(1e-3, list_hyper$phi)
      list_hyper$xi <- c(1e-3, list_hyper$xi)
      
      list_init$alpha_vb <- rbind(rep(0, d), list_init$alpha_vb)
      
      if (link == "probit") {
        
        list_init$sig2_alpha_vb <- c(1, list_init$sig2_alpha_vb)
        
      } else {
        
        list_init$sig2_alpha_vb <- rbind(rep(1, d), list_init$sig2_alpha_vb)
        
      }
      
      
    }
    colnames(Z)[1] <- "Intercept"
  }
  
  
  if (verbose){
    cat(paste0("============================================================== \n",
               "== Variational inference for sparse multivariate regression == \n",
               "============================================================== \n\n"))
  }
  
  
  if (is.null(list_blocks)) {
    
    if (link == "identity") {
      
      ng  <- is.null(list_groups)
      ns <- is.null(list_struct)
      
      if (ng & ns){
        
        if (nq) {
          
          vb <- locus_core_(Y, X, list_hyper, list_init$gam_vb,
                            list_init$mu_beta_vb, list_init$sig2_beta_vb,
                            list_init$tau_vb, tol, maxit, anneal, verbose, 
                            checkpoint_path = checkpoint_path, 
                            full_output = full_output)
          
        } else { # q non-null
          
          vb <- locus_z_core_(Y, X, Z, list_hyper, list_init$gam_vb,
                              list_init$alpha_vb, list_init$mu_beta_vb,
                              list_init$sig2_alpha_vb, list_init$sig2_beta_vb,
                              list_init$tau_vb, tol, maxit, anneal, verbose, 
                              full_output = full_output)
          
        } 
        
      } else if (!ng){
        
        # X is a list (transformed in prepare_data)
        # mu_beta_vb is a list (transformed in prepare_init)
        vb <- locus_group_core_(Y, X, list_hyper, list_init$gam_vb,
                                list_init$mu_beta_vb, list_init$sig2_inv_vb,
                                list_init$tau_vb, tol, maxit, verbose, 
                                full_output = full_output)
        
      } else { # list_struct non-null, and only predictor propensity control.
        
        vb <- locus_struct_core_(Y, X, list_hyper, list_init$gam_vb,
                                 list_init$mu_beta_vb, list_init$sig2_beta_vb,
                                 list_init$tau_vb, list_struct, tol, maxit, 
                                 verbose, full_output = full_output)
      }
      
    } else if (link == "logit"){
      
      
      vb <- locus_logit_core_(Y, X, Z, list_hyper, list_init$chi_vb,
                              list_init$gam_vb, list_init$alpha_vb,
                              list_init$mu_beta_vb, list_init$sig2_alpha_vb,
                              list_init$sig2_beta_vb, tol, maxit, verbose, 
                              full_output = full_output)
      
      
    } else if (link == "probit"){
      
      vb <- locus_probit_core_(Y, X, Z, list_hyper, list_init$gam_vb,
                               list_init$alpha_vb, list_init$mu_beta_vb,
                               list_init$sig2_alpha_vb, list_init$sig2_beta_vb,
                               tol, maxit, verbose, full_output = full_output)
      
    } else {
      
      vb <- locus_mix_core_(Y, X, Z, ind_bin, list_hyper, list_init$gam_vb,
                            list_init$alpha_vb, list_init$mu_beta_vb,
                            list_init$sig2_alpha_vb, list_init$sig2_beta_vb,
                            list_init$tau_vb, tol, maxit, verbose, 
                            full_output = full_output)
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
        
        if (nq) {
          
          vb_bl <- locus_core_(Y, X_bl, list_hyper_bl,
                               list_init_bl$gam_vb, list_init_bl$mu_beta_vb,
                               list_init_bl$sig2_beta_vb, list_init_bl$tau_vb,
                               tol, maxit, anneal, verbose = FALSE, 
                               full_output = full_output)
          
        } else {
          
          vb_bl <- locus_z_core_(Y, X_bl, Z, list_hyper_bl, list_init_bl$gam_vb,
                                 list_init_bl$alpha_vb,list_init_bl$mu_beta_vb,
                                 list_init_bl$sig2_alpha_vb,
                                 list_init_bl$sig2_beta_vb, list_init_bl$tau_vb,
                                 tol, maxit, anneal, verbose = FALSE, 
                                 full_output = full_output)
          
        }
        
      } else if (link == "logit") {
        
        vb_bl <- locus_logit_core_(Y, X_bl, Z, list_hyper_bl,
                                   list_init_bl$chi_vb, list_init_bl$gam_vb,
                                   list_init_bl$alpha_vb, list_init_bl$mu_beta_vb,
                                   list_init_bl$sig2_alpha_vb,
                                   list_init_bl$sig2_beta_vb, tol, maxit,
                                   verbose = FALSE, full_output = full_output)
        
        
      } else  if (link == "probit") {
        
        vb_bl <- locus_probit_core_(Y, X_bl, Z, list_hyper_bl,
                                    list_init_bl$gam_vb, list_init_bl$alpha_vb,
                                    list_init_bl$mu_beta_vb,
                                    list_init_bl$sig2_alpha_vb,
                                    list_init_bl$sig2_beta_vb, tol, maxit,
                                    verbose = FALSE, full_output = full_output)
        
      } else {
        
        vb_bl <- locus_mix_core_(Y, X_bl, Z, ind_bin, list_hyper_bl,
                                 list_init_bl$gam_vb, list_init_bl$alpha_vb,
                                 list_init_bl$mu_beta_vb,
                                 list_init_bl$sig2_alpha_vb,
                                 list_init_bl$sig2_beta_vb, list_init_bl$tau_vb,
                                 tol, maxit, verbose = FALSE, 
                                 full_output = full_output)
        
        if (verbose) {
          if (vb_bl$converged) {
            cat(paste0("The algorithm reached convergence on one block after ",
                       format(vb_bl$it), " iterations. \n", "Optimal marginal ",
                       "log-likelihood variational lower bound (ELBO) = ",
                       format(vb_bl$lb_opt), ". \n\n"))
          } else {
            cat(paste0("The algorithm reached the maximal number of iterations ",
                       "on one block before converging.\n", "Difference in ELBO ",
                       "between last and penultimate iterations: ",
                       format(vb_bl$diff_lb), ".\n\n"))
          }
        }
        
      } 
      
      vb_bl
    }
    
    list_vb <- parallel::mclapply(1:n_bl, function(k) locus_bl_(k), mc.cores = n_cpus)
    
    names_vec <- c("converged", "it", "lb_opt")
    names_vec <- c(names_vec, "om_vb")
    
    names_mat <- "gam_vb"
    
    vb <- c(lapply(names_vec, function(key) {
      vec <- do.call(c, lapply(list_vb, `[[`, key))
      names(vec) <- paste0("bl_", 1:n_bl)
      vec}),
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
  
  if (!is.null(list_groups)) {
    vb$group_labels <- vec_fac_gr # after removal of constant or collinear covariates
  }
  
  if (save_hyper) vb$list_hyper <- list_hyper
  if (save_init) vb$list_init <- list_init
  
  class(vb) <- "locus"
  
  
  if (verbose) cat("... done. == \n\n")
  
  vb
  
}
