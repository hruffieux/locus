# This file is part of the `locus` R package:
#     https://github.com/hruffieux/locus
#
# Internal core function to call the variational algorithm for dual propensity
# control. Sparse regression with identity link, no fixed covariates.
# See help of `locus` function for details.
#
locus_dual_horseshoe_core_ <- function(Y, X, list_hyper, gam_vb, mu_beta_vb, sig2_beta_vb,
                                       tau_vb, list_struct, tol, maxit, anneal, verbose,
                                       batch = "y", full_output = FALSE, debug = TRUE) {

  # Y must have been centered, and X standardized.
  
  d <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  
  with(list_hyper, { # list_init not used with the with() function to avoid
                     # copy-on-write for large objects
    
    # Preparing annealing if any
    #
    anneal_scale <- TRUE # if TRUE, scale parameters s02 and b_vb also annealed.
    
    if (is.null(anneal)) {
      annealing <- FALSE
      c <- c_s <- 1 # c_s for scale parameters
    } else {
      annealing <- TRUE
      ladder <- get_annealing_ladder_(anneal, verbose)
      c <- ladder[1]
      c_s <- ifelse(anneal_scale, c, 1)
    }
    
    eps <- .Machine$double.eps^0.5
    
    # Covariate-specific parameters: objects derived from s02, list_struct (possible block-wise in parallel)
    #
    if (is.null(list_struct)) { # /! here list_struct is used to define block-specific s02, not to inject structure!
      n_bl <- 1
      bl_lgths <- p
    } else {
      vec_fac_bl <- list_struct$vec_fac_st
      bl_ids <- list_struct$bl_ids <- unique(vec_fac_bl)
      n_bl <- list_struct$n_bl <- length(bl_ids) 
      bl_lgths <- list_struct$bl_lgths <- table(vec_fac_bl)
    }
    
    # Variance initialization
    #
    S0_inv_vb <- rgamma(n_bl, shape = (bl_lgths + 1) / 2, rate = 1) # initial guess
    
    
    # Some hyperparameters
    #
    A2_inv <- 1 #  hyperparameter # TODO: see how to fix, sensitivity analysis
    
    # Choose m0 so that, `a priori' (i.e. before optimization), E_p_gam is as specified by the user:
    n0_star <- - n0[1]
    m0_star <- n0_star * (sqrt(1 + (1/S0_inv_vb[1] / d) / (1 + t02)) - 1) # assumes b equiv 1. see hyperparameter_setting document. m0 not needed anymore in set_hyper.
    
    m0 <- - rep(m0_star, p)
    
    # Parameter initialization here for the top level 
    #
    mu_theta_vb <- rnorm(p, mean = m0, sd = abs(m0) / 5) # m0 ########################### see how to set m0 
    sig2_theta_vb <- 1 / (d + rgamma(p, shape = S0_inv_vb[1] * d, rate = 1)) # initial guess assuming b_vb = 1

    mu_rho_vb <- rnorm(d, mean = n0, sd = abs(n0) / 5) # n0
    
    
    # Response-specific parameters: objects derived from t02
    #
    T0_inv <- 1 / t02
    sig2_rho_vb <- update_sig2_c0_vb_(p, t02, c = c) # stands for a diagonal matrix of size d with this value on the (constant) diagonal
    
    vec_sum_log_det_rho <- - d * (log(t02) + log(p + T0_inv))
    
    
    # Stored/precomputed objects
    #
    m1_beta <- update_m1_beta_(gam_vb, mu_beta_vb)
    m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE)
    
    mat_x_m1 <- update_mat_x_m1_(X, m1_beta)
    
    # Fixed VB parameter
    #
    lambda_a_inv_vb <- 1 # no change with annealing 
    
    converged <- FALSE
    lb_new <- -Inf
    it <- 0
    
    
    while ((!converged) & (it < maxit)) {
      
      lb_old <- lb_new
      it <- it + 1
      
      if (verbose & (it == 1 | it %% 5 == 0))
        cat(paste("Iteration ", format(it), "... \n", sep = ""))
      
      # % #
      lambda_vb <- update_lambda_vb_(lambda, sum(gam_vb), c = c)
      nu_vb <- update_nu_vb_(nu, m2_beta, tau_vb, c = c)
      
      sig2_inv_vb <- lambda_vb / nu_vb
      # % #
      
      # % #
      eta_vb <- update_eta_vb_(n, eta, gam_vb, c = c)
      kappa_vb <- update_kappa_vb_(Y, kappa, mat_x_m1, m1_beta, m2_beta, sig2_inv_vb, c = c)
      
      tau_vb <- eta_vb / kappa_vb
      # % #
      
      sig2_beta_vb <- update_sig2_beta_vb_(n, sig2_inv_vb, tau_vb, c = c)
      
      log_tau_vb <- update_log_tau_vb_(eta_vb, kappa_vb)
      log_sig2_inv_vb <- update_log_sig2_inv_vb_(lambda_vb, nu_vb)
      
      
      # different possible batch-coordinate ascent schemes:
      #
      if (batch == "y") { # optimal scheme
        
        log_Phi_mu_theta_plus_rho <- sapply(mu_rho_vb, function(mu_rho_k) {
          pnorm(mu_theta_vb + mu_rho_k, log.p = TRUE)})
        
        log_1_min_Phi_mu_theta_plus_rho <- sapply(mu_rho_vb, function(mu_rho_k) {
          pnorm(mu_theta_vb + mu_rho_k, lower.tail = FALSE, log.p = TRUE)})
        
        # C++ Eigen call for expensive updates
        shuffled_ind <- as.numeric(sample(0:(p-1))) # Zero-based index in C++
        
        coreDualLoop(X, Y, gam_vb, log_Phi_mu_theta_plus_rho,
                     log_1_min_Phi_mu_theta_plus_rho, log_sig2_inv_vb,
                     log_tau_vb, m1_beta, mat_x_m1, mu_beta_vb,
                     sig2_beta_vb, tau_vb, shuffled_ind, c = c)
        
      } else if (batch == "0"){ # no batch, used only internally
                                # schemes "x" of "x-y" are not batch concave
                                # hence not implemented as they may diverge
        
        for (k in sample(1:d)) {
          
          for (j in sample(1:p)) {
            
            mat_x_m1[, k] <- mat_x_m1[, k] - X[, j] * m1_beta[j, k]
            
            mu_beta_vb[j, k] <- c * sig2_beta_vb[k] * tau_vb[k] * crossprod(Y[, k] - mat_x_m1[, k], X[, j])
            
            gam_vb[j, k] <- exp(-log_one_plus_exp_(c * (pnorm(mu_theta_vb[j] + mu_rho_vb[k], lower.tail = FALSE, log.p = TRUE) -
                                                          pnorm(mu_theta_vb[j] + mu_rho_vb[k], log.p = TRUE) -
                                                          log_tau_vb[k] / 2 - log_sig2_inv_vb / 2 -
                                                          mu_beta_vb[j, k] ^ 2 / (2 * sig2_beta_vb[k]) -
                                                          log(sig2_beta_vb[k]) / 2)))
            
            m1_beta[j, k] <- gam_vb[j, k] * mu_beta_vb[j, k]
            
            mat_x_m1[, k] <- mat_x_m1[, k] + X[, j] * m1_beta[j, k]
            
          }
        }
        
      } else {
        
        stop ("Batch scheme not defined. Exit.")
        
      }
      
      
      m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE)
      
      W <- update_W_info_(gam_vb, sweep(tcrossprod(mu_theta_vb, rep(1, d)), 2,
                                        mu_rho_vb, `+`), c = c) # we use info_ so that the second argument is a matrix
      
      # keep this order!
      #
      
      if (is.null(list_struct)) {
        
        G_vb <- c_s * S0_inv_vb * d * (mu_theta_vb^2 + sig2_theta_vb - 2 * mu_theta_vb * m0 + m0^2) / 2 # b_vb not annealed. possibly no closed form.
        nu_a_inv_vb <- c_s * (A2_inv + S0_inv_vb)
        
      } else {
        
        G_vb <- unlist(lapply(1:n_bl, function(bl) { c_s * S0_inv_vb[bl] * 
            (mu_theta_vb[vec_fac_bl == bl_ids[bl]]^2 + sig2_theta_vb[vec_fac_bl == bl_ids[bl]] - 
               2 * mu_theta_vb[vec_fac_bl == bl_ids[bl]] * m0[vec_fac_bl == bl_ids[bl]] + 
               m0[vec_fac_bl == bl_ids[bl]]^2) / 2 }))
        
        nu_a_inv_vb <- sapply(1:n_bl, function(bl) c_s * (A2_inv + S0_inv_vb[bl]))
        
      } 
      
      
      if (annealing & anneal_scale) {
        
        b_vb <- gsl::gamma_inc(- c_s + 2, G_vb) / (gsl::gamma_inc(- c_s + 1, G_vb) * G_vb) - 1
        
      } else {
        
        b_vb <- 1 / (sapply(G_vb, function(G_vb_s) Q_approx(G_vb_s)) * G_vb) - 1 # TODO implement a Q_approx for vectors
        
      }
      
      a_inv_vb <- lambda_a_inv_vb / nu_a_inv_vb
      

      if (is.null(list_struct)) {
        
        sig2_theta_vb <- update_sig2_c0_vb_(d, 1 / (S0_inv_vb * b_vb * d), c = c)
        mu_theta_vb <- update_mu_theta_vb_(W, m0, S0_inv_vb * b_vb * d, sig2_theta_vb,
                                           vec_fac_st = NULL, mu_rho_vb, is_mat = FALSE, c = c)
        
        lambda_s0_vb <- update_lambda_vb_(1 / 2, p, c = c_s)
        
        nu_s0_vb <- c_s * (a_inv_vb + 
                           sum(b_vb * d * (mu_theta_vb^2 + sig2_theta_vb - 2 * mu_theta_vb * m0 + m0^2)) / 2) 
        
      } else {
        
        sig2_theta_vb <- unlist(lapply(1:n_bl, function(bl) { 
          update_sig2_c0_vb_(d, 1 / (S0_inv_vb[bl] * b_vb[vec_fac_bl == bl_ids[bl]] * d), c = c) }))
        
        mu_theta_vb <- unlist(lapply(1:n_bl, function(bl) {
          update_mu_theta_vb_(W[vec_fac_bl == bl_ids[bl], , drop = FALSE], m0[vec_fac_bl == bl_ids[bl]], 
                              S0_inv_vb[bl] * b_vb[vec_fac_bl == bl_ids[bl]] * d, sig2_theta_vb[vec_fac_bl == bl_ids[bl]],
                              vec_fac_st = NULL, mu_rho_vb, is_mat = FALSE, c = c)
        }))
        
        lambda_s0_vb <- sapply(1:n_bl, function(bl) update_lambda_vb_(1 / 2, bl_lgths[bl], c = c_s))
        
        nu_s0_vb <- sapply(1:n_bl, function(bl) { c_s * (a_inv_vb[bl] + 
                           sum(b_vb[vec_fac_bl == bl_ids[bl]] * d * 
                                 (mu_theta_vb[vec_fac_bl == bl_ids[bl]]^2 + 
                                    sig2_theta_vb[vec_fac_bl == bl_ids[bl]] - 
                                    2 * mu_theta_vb[vec_fac_bl == bl_ids[bl]] * m0[vec_fac_bl == bl_ids[bl]] + 
                                    m0[vec_fac_bl == bl_ids[bl]]^2)) / 2)}) 
        
      }
      
      S0_inv_vb <- as.numeric(lambda_s0_vb / nu_s0_vb)
      
      mu_rho_vb <- update_mu_rho_vb_(W, mu_theta_vb, n0, sig2_rho_vb, T0_inv,
                                     is_mat = FALSE, c = c) 
      

      if (verbose & (it == 1 | it %% 5 == 0)) {
        
        if (is.null(list_struct)) {
          cat(paste0("Updated s02 / d: ", format(1 / S0_inv_vb / d, digits = 4), ".\n"))
          cat("Updated 1 / b: \n")
          print(summary(1 / b_vb))
          cat("\n")
        } else {
          cat("Updated block-specific s02 / d: \n")
          print(summary(1 / S0_inv_vb / d))
          cat("\n")
          cat("Updated 1 / b: \n")
          print(summary(1 / b_vb))
          cat("\n")
        }
        
      }
      
      if (annealing) {
        
        if (verbose & (it == 1 | it %% 5 == 0))
          cat(paste("Temperature = ", format(1 / c, digits = 4), "\n\n", sep = ""))
        
        sig2_rho_vb <- c * sig2_rho_vb
        
        c <- ifelse(it < length(ladder), ladder[it + 1], 1)
        c_s <- ifelse(anneal_scale, c, 1)
        
        sig2_rho_vb <- sig2_rho_vb / c
        
        if (isTRUE(all.equal(c, 1))) {
          
          annealing <- FALSE
          
          if (verbose)
            cat("** Exiting annealing mode. **\n\n")
        }
        
        
      } else {
        
        lb_new <- elbo_dual_horseshoe_(Y, a_inv_vb, A2_inv, b_vb, eta, eta_vb, G_vb, gam_vb, kappa, kappa_vb, lambda,
                                       lambda_vb, lambda_a_inv_vb, lambda_s0_vb, m0, n0, mu_rho_vb,
                                       mu_theta_vb, nu, nu_vb, nu_a_inv_vb, nu_s0_vb, sig2_beta_vb,
                                       S0_inv_vb, sig2_theta_vb, sig2_inv_vb, sig2_rho_vb,
                                       T0_inv, tau_vb, m1_beta, m2_beta, mat_x_m1,
                                       vec_sum_log_det_rho, list_struct)
        
        if (verbose & (it == 1 | it %% 5 == 0))
          cat(paste("ELBO = ", format(lb_new), "\n\n", sep = ""))
        
        
        if (debug && lb_new + eps < lb_old)
          stop("ELBO not increasing monotonically. Exit. ")
        
        converged <- (abs(lb_new - lb_old) < tol)
        
      }
      
    }
    
    if (verbose) {
      if (converged) {
        cat(paste("Convergence obtained after ", format(it), " iterations. \n",
                  "Optimal marginal log-likelihood variational lower bound ",
                  "(ELBO) = ", format(lb_new), ". \n\n", sep = ""))
      } else {
        warning("Maximal number of iterations reached before convergence. Exit.")
      }
    }
    
    lb_opt <- lb_new
    
    if (full_output) { # for internal use only
      
      create_named_list_(a_inv_vb, A2_inv, b_vb, eta, eta_vb, G_vb, gam_vb, kappa, kappa_vb, lambda,
                         lambda_vb, lambda_a_inv_vb, lambda_s0_vb, m0, n0, mu_rho_vb,
                         mu_theta_vb, nu, nu_vb, nu_a_inv_vb, nu_s0_vb, sig2_beta_vb,
                         S0_inv_vb, sig2_theta_vb, sig2_inv_vb, sig2_rho_vb,
                         T0_inv, tau_vb, m1_beta, m2_beta, mat_x_m1,
                         vec_sum_log_det_rho, list_struct)
      
    } else {
      
      names_x <- colnames(X)
      names_y <- colnames(Y)
      
      rownames(gam_vb) <- names_x
      colnames(gam_vb) <- names_y
      names(mu_theta_vb) <- names_x
      names(mu_rho_vb) <- names_y
      names(b_vb) <- names_x
      
      diff_lb <- abs(lb_opt - lb_old)
      
      create_named_list_(gam_vb, mu_theta_vb, mu_rho_vb, converged, it, lb_opt,
                         diff_lb, S0_inv_vb, b_vb)
      
    }
  })
  
}



# Internal function which implements the marginal log-likelihood variational
# lower bound (ELBO) corresponding to the `locus_struct_core` algorithm.
#
elbo_dual_horseshoe_ <- function(Y, a_inv_vb, A2_inv, b_vb, eta, eta_vb, G_vb, gam_vb, kappa, kappa_vb, lambda,
                                 lambda_vb, lambda_a_inv_vb, lambda_s0_vb, m0, n0, mu_rho_vb,
                                 mu_theta_vb, nu, nu_vb, nu_a_inv_vb, nu_s0_vb, sig2_beta_vb,
                                 S0_inv_vb, sig2_theta_vb, sig2_inv_vb, sig2_rho_vb,
                                 T0_inv, tau_vb, m1_beta, m2_beta, mat_x_m1,
                                 vec_sum_log_det_rho, list_struct) {
  
  n <- nrow(Y)
  d <- nrow(Y)
  p <- length(m0)
  
  # needed for monotonically increasing elbo.
  #
  eta_vb <- update_eta_vb_(n, eta, gam_vb)
  kappa_vb <- update_kappa_vb_(Y, kappa, mat_x_m1, m1_beta, m2_beta, sig2_inv_vb)
  
  lambda_vb <- update_lambda_vb_(lambda, sum(gam_vb))
  nu_vb <- update_nu_vb_(nu, m2_beta, tau_vb)
  
  log_tau_vb <- update_log_tau_vb_(eta_vb, kappa_vb)
  log_sig2_inv_vb <- update_log_sig2_inv_vb_(lambda_vb, nu_vb)
  
  log_S0_inv_vb <- update_log_sig2_inv_vb_(lambda_s0_vb, nu_s0_vb)
  log_a_inv_vb <- update_log_sig2_inv_vb_(lambda_a_inv_vb, nu_a_inv_vb)
  
  if (!is.null(list_struct)) {
    n_bl <- list_struct$n_bl
    bl_ids <- list_struct$bl_ids
    bl_lgths <- list_struct$bl_lgths
    vec_fac_bl <- list_struct$vec_fac_st
  }
  
  
  elbo_A <- e_y_(n, kappa, kappa_vb, log_tau_vb, m2_beta, sig2_inv_vb, tau_vb)
  
  
  if (is.null(list_struct)) {
    elbo_B <- e_beta_gamma_dual_(gam_vb, log_sig2_inv_vb, log_tau_vb,
                                 mu_rho_vb, mu_theta_vb, m2_beta,
                                 sig2_beta_vb, sig2_rho_vb,
                                 sig2_theta_vb, sig2_inv_vb, tau_vb)
    
    elbo_C <- e_theta_hs_(b_vb, G_vb, log_S0_inv_vb + log(d), m0, mu_theta_vb, S0_inv_vb * d, sig2_theta_vb)
    
  } else {
    elbo_B <- sum(sapply(1:n_bl, function(bl) {
      e_beta_gamma_dual_(gam_vb[vec_fac_bl == bl_ids[bl], , drop = FALSE], log_sig2_inv_vb, log_tau_vb,
                         mu_rho_vb, mu_theta_vb[vec_fac_bl == bl_ids[bl]], 
                         m2_beta[vec_fac_bl == bl_ids[bl], , drop = FALSE],
                         sig2_beta_vb, sig2_rho_vb,
                         sig2_theta_vb[vec_fac_bl == bl_ids[bl]], sig2_inv_vb, tau_vb)}))
    
    elbo_C <- sum(sapply(1:n_bl, function(bl) { 
      e_theta_hs_(b_vb[vec_fac_bl == bl_ids[bl]], G_vb[vec_fac_bl == bl_ids[bl]], log_S0_inv_vb[bl] + log(d), 
                  m0[vec_fac_bl == bl_ids[bl]], mu_theta_vb[vec_fac_bl == bl_ids[bl]], S0_inv_vb[bl] * d, 
                  sig2_theta_vb[vec_fac_bl == bl_ids[bl]])}))
  }
  
  elbo_D <- e_rho_(mu_rho_vb, n0, sig2_rho_vb, T0_inv, vec_sum_log_det_rho)
  
  elbo_E <- e_tau_(eta, eta_vb, kappa, kappa_vb, log_tau_vb, tau_vb)
  
  elbo_F <- sum(e_sig2_inv_hs_(a_inv_vb, lambda_s0_vb, log_a_inv_vb, log_S0_inv_vb, nu_s0_vb, S0_inv_vb)) # S0_inv_vb
  
  elbo_G <- sum(e_sig2_inv_(1 / 2, lambda_a_inv_vb, log_a_inv_vb, A2_inv, nu_a_inv_vb, a_inv_vb)) # a_inv_vb
  
  elbo_H <- e_sig2_inv_(lambda, lambda_vb, log_sig2_inv_vb, nu, nu_vb, sig2_inv_vb)
  
  elbo_A + elbo_B + elbo_C + elbo_D + elbo_E + elbo_F + elbo_G + elbo_H
  
}

