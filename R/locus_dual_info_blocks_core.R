# This file is part of the `locus` R package:
#     https://github.com/hruffieux/locus
#
# Internal core function to call the variational algorithm for dual propensity
# control with external information variables. Sparse regression with identity
# link, no fixed covariates. See help of `locus` function for details.
#
locus_dual_info_blocks_core_ <- function(Y, X, list_V, vec_fac_bl, list_hyper, gam_vb, mu_beta_vb,
                                         sig2_beta_vb, tau_vb, list_struct, tol, maxit,
                                         anneal, verbose, batch = "y", full_output = FALSE, debug = TRUE) {
  
  # Y centered, and X and V standardized.
  
  d <- ncol(Y)
  n <- nrow(Y)
  p <- ncol(X)
  
  with(list_hyper, { # list_init not used with the with() function to avoid
                     # copy-on-write for large objects
    
    # Preparing annealing if any
    #
    if (is.null(anneal)) {
      annealing <- FALSE
      c <- 1
    } else {
      annealing <- TRUE
      ladder <- get_annealing_ladder_(anneal, verbose)
      c <- ladder[1]
    }
    
    eps <- .Machine$double.eps^0.5
    
    bl_ids <- unique(vec_fac_bl)
    n_bl <- length(bl_ids)

    vec_p_bl <- table(vec_fac_bl)
    vec_r_bl <- sapply(list_V, function(V) ncol(V))# size of the V_bl (after removal of cst and coll annotations in each)
    
    # sanity checks
    #
    stopifnot(length(om_vb) == n_bl & length(s2) == n_bl)
    stopifnot(all(sapply(1:n_bl, function(bl) length(om_vb[[bl]]) == vec_r_bl[bl])))
    stopifnot(sum(sapply(list_V, function(V) nrow(V))) == p)
    
    # Parameter initialization here for the top level only
    #
    mu_theta_vb <- rnorm(p, mean = m0, sd = abs(m0) / 5) 
    mu_rho_vb <- rnorm(d, mean = n0, sd = abs(n0) / 5) 
    
    
    mu_c_vb <- lapply(vec_r_bl, function(r) rnorm(r, sd = 0.1)) 
    
    # om_vb is a list of length n_bl
    zeta_vb <- lapply(om_vb, function(om_bl) rbeta(length(om_bl), shape1 = om_bl + eps, shape2 = 1 - om_bl + eps))
    
    log_om_vb <- lapply(om_vb, function(om_bl) log(om_bl + eps))
    log_1_min_om_vb <- lapply(om_vb, function(om_bl) log(1 - om_bl + eps))
    
    
    # Covariate-specific parameters: objects derived from s02, list_struct (possible block-wise in parallel)
    #
    obj_theta_vb <- update_sig2_theta_vb_(d, p, list_struct, s02, X, c = c)
    
    S0_inv <- obj_theta_vb$S0_inv
    sig2_theta_vb <- obj_theta_vb$sig2_theta_vb
    vec_sum_log_det_theta <- obj_theta_vb$vec_sum_log_det_theta
    
    vec_fac_st <- obj_theta_vb$vec_fac_st
    
    
    # Response-specific parameters: objects derived from t02
    #
    T0_inv <- 1 / t02
    sig2_rho_vb <- update_sig2_c0_vb_(p, t02, c = c) # stands for a diagonal matrix of size d with this value on the (constant) diagonal
    vec_sum_log_det_rho <- - d * (log(t02) + log(p + T0_inv))
    
    
    # External information effects
    #
    sig2_c_vb <- sapply(1:n_bl, function(bl) update_sig2_c_vb_(vec_p_bl[bl], s2[bl], d, c = c))
    
    # Stored/precomputed objects
    #
    m1_beta <- update_m1_beta_(gam_vb, mu_beta_vb)
    m2_beta <- update_m2_beta_(gam_vb, mu_beta_vb, sig2_beta_vb, sweep = TRUE)
    m1_c <- lapply(1:n_bl, function(bl) update_m1_beta_(zeta_vb[[bl]], mu_c_vb[[bl]])) # matrix of size r x n_bl
    
    mat_x_m1 <- update_mat_x_m1_(X, m1_beta)
    mat_v_mu <- update_mat_v_mu_block_(list_V, mu_theta_vb, mu_rho_vb, m1_c, vec_fac_bl)
    
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
      
      if (batch == "y") { # optimal scheme
        
        log_Phi_mat_v_mu <- pnorm(mat_v_mu, log.p = TRUE)
        
        log_1_min_Phi_mat_v_mu <- pnorm(mat_v_mu, lower.tail = FALSE, log.p = TRUE)
        
        # C++ Eigen call for expensive updates
        shuffled_ind <- as.numeric(sample(0:(p-1))) # Zero-based index in C++
        
        coreDualLoop(X, Y, gam_vb, log_Phi_mat_v_mu,
                     log_1_min_Phi_mat_v_mu, log_sig2_inv_vb,
                     log_tau_vb, m1_beta, mat_x_m1, mu_beta_vb,
                     sig2_beta_vb, tau_vb, shuffled_ind, c = c)
        
      } else if (batch == "0"){ # no batch, used only internally
        # schemes "x" of "x-y" are not batch concave
        # hence not implemented as they may diverge
        
        for (k in sample(1:d)) {
          
          for (j in sample(1:p)) {
            
            mat_x_m1[, k] <- mat_x_m1[, k] - X[, j] * m1_beta[j, k]
            
            mu_beta_vb[j, k] <- c * sig2_beta_vb[k] * tau_vb[k] * crossprod(Y[, k] - mat_x_m1[, k], X[, j])
            
            gam_vb[j, k] <- exp(-log_one_plus_exp_(c * (pnorm(mat_v_mu[j, k], lower.tail = FALSE, log.p = TRUE) -
                                                          pnorm(mat_v_mu[j, k], log.p = TRUE) -
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
      
      W <- update_W_info_(gam_vb, mat_v_mu, c = c) # we use info_ so that the second argument is a matrix
      
      mat_v_mu <- sweep(mat_v_mu, 1, mu_theta_vb, `-`)
      mu_theta_vb <- update_mu_theta_vb_(W, m0, S0_inv, sig2_theta_vb,
                                         vec_fac_st, mat_v_mu, is_mat = TRUE, c = c)
      
      mat_v_mu <- sweep(sweep(mat_v_mu, 1, mu_theta_vb, `+`), 2, mu_rho_vb, `-`)
      
      mu_rho_vb <- update_mu_rho_vb_(W, mat_v_mu, n0, sig2_rho_vb, T0_inv, is_mat = TRUE, c = c)
      mat_v_mu <- sweep(mat_v_mu, 2, mu_rho_vb, `+`)
        
      for (bl in sample(1:n_bl)) {
      
        for (l in sample(1:vec_r_bl[bl])) {
          
          mat_v_mu[vec_fac_bl == bl_ids[bl], ] <- sweep(mat_v_mu[vec_fac_bl == bl_ids[bl], , drop = FALSE], 1,
                                                        list_V[[bl]][, l] * m1_c[[bl]][l], `-`)
          
          mu_c_vb[[bl]][l] <- c * sig2_c_vb[bl] *
            sum(crossprod(W[vec_fac_bl == bl_ids[bl], , drop = FALSE] - mat_v_mu[vec_fac_bl == bl_ids[bl], , drop = FALSE], list_V[[bl]][, l]))
          
          zeta_vb[[bl]][l] <- exp(-log_one_plus_exp_(c * (log_1_min_om_vb[[bl]][l] - log_om_vb[[bl]][l] +
                                                                 log(s2[bl]) / 2 - log(sig2_c_vb[bl]) / 2 -
                                                                 mu_c_vb[[bl]][l] ^ 2 / (2 * sig2_c_vb[bl]))))
          
          m1_c[[bl]][l] <- mu_c_vb[[bl]][l] * zeta_vb[[bl]][l]
          
          
          mat_v_mu[vec_fac_bl == bl_ids[bl], ] <- sweep(mat_v_mu[vec_fac_bl == bl_ids[bl], , drop = FALSE], 1,
                                                        list_V[[bl]][, l] * m1_c[[bl]][l], `+`)
          
        }
        
      }
      
      
      if (annealing) {
        
        if (verbose & (it == 1 | it %% 5 == 0))
          cat(paste("Temperature = ", format(1 / c, digits = 4), "\n\n", sep = ""))
        
        sig2_theta_vb <- c * sig2_theta_vb
        sig2_rho_vb <- c * sig2_rho_vb
        sig2_c_vb <- c * sig2_c_vb
        
        c <- ifelse(it < length(ladder), ladder[it + 1], 1)
        
        sig2_theta_vb <- sig2_theta_vb / c
        sig2_rho_vb <- sig2_rho_vb / c
        sig2_c_vb <- sig2_c_vb / c
        
        if (isTRUE(all.equal(c, 1))) {
          
          annealing <- FALSE
          
          if (verbose)
            cat("** Exiting annealing mode. **\n\n")
        }
        
      } else {
        
        lb_new <- elbo_dual_info_blocks_(Y, list_V, eta, eta_vb, gam_vb, kappa, kappa_vb, lambda,
                                  lambda_vb, log_1_min_om_vb, log_om_vb, m0, n0, mu_c_vb, 
                                  mu_rho_vb, mu_theta_vb, nu, nu_vb, 
                                  sig2_beta_vb, S0_inv, s2, sig2_c_vb, sig2_theta_vb,
                                  sig2_inv_vb, sig2_rho_vb, T0_inv, tau_vb, zeta_vb, m1_beta,
                                  m2_beta, mat_x_m1, mat_v_mu, vec_fac_st, vec_sum_log_det_rho,
                                  vec_sum_log_det_theta)

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
      
      create_named_list_(Y, list_V, eta, eta_vb, gam_vb, kappa, kappa_vb, lambda,
                         lambda_vb, m0, n0, mu_c_vb, mu_rho_vb, mu_theta_vb, nu, nu_vb, om_vb,
                         sig2_beta_vb, S0_inv, s2, sig2_c_vb, sig2_theta_vb,
                         sig2_inv_vb, sig2_rho_vb, T0_inv, tau_vb, zeta_vb, m1_beta,
                         m2_beta, mat_x_m1, mat_v_mu, vec_fac_st, vec_sum_log_det_rho,
                         vec_sum_log_det_theta)
      
    } else {
      
      names_x <- colnames(X)
      names_y <- colnames(Y)
      names_v <- lapply(list_V, function(V) colnames(V))
      
      rownames(gam_vb) <- names_x
      colnames(gam_vb) <- names_y
      
      names(mu_theta_vb) <- names_x
      names(mu_rho_vb) <- names_y
      
      mu_c_vb <- lapply(1:n_bl, function(bl) {
        names(mu_c_vb[[bl]]) <- colnames(list_V[[bl]])
        mu_c_vb[[bl]]})
      
      om_vb <- lapply(1:n_bl, function(bl) {
        names(om_vb[[bl]]) <- colnames(list_V[[bl]])
        om_vb[[bl]]})
      
      zeta_vb <- lapply(1:n_bl, function(bl) {
        names(zeta_vb[[bl]]) <- colnames(list_V[[bl]])
        zeta_vb[[bl]]})
        
      names(zeta_vb) <- names(mu_c_vb) <- names(om_vb) <- paste0("bl_", 1:n_bl)
      
      diff_lb <- abs(lb_opt - lb_old)
      
      create_named_list_(mu_c_vb, om_vb, gam_vb, mu_theta_vb, mu_rho_vb, zeta_vb, converged, it,
                         lb_opt, diff_lb)
      
    }
  })
  
}



# Internal function which implements the marginal log-likelihood variational
# lower bound (ELBO) corresponding to the `locus_struct_core` algorithm.
#
elbo_dual_info_blocks_ <- function(Y, list_V, eta, eta_vb, gam_vb, kappa, kappa_vb, lambda,
                                   lambda_vb, log_1_min_om_vb, log_om_vb, m0, n0, mu_c_vb, 
                                   mu_rho_vb, mu_theta_vb, nu, nu_vb, 
                                   sig2_beta_vb, S0_inv, s2, sig2_c_vb, sig2_theta_vb,
                                   sig2_inv_vb, sig2_rho_vb, T0_inv, tau_vb, zeta_vb, m1_beta,
                                   m2_beta, mat_x_m1, mat_v_mu, vec_fac_st, vec_sum_log_det_rho,
                                   vec_sum_log_det_theta) {
  
  n <- nrow(Y)
  
  # needed for monotonically increasing elbo.
  #
  eta_vb <- update_eta_vb_(n, eta, gam_vb)
  kappa_vb <- update_kappa_vb_(Y, kappa, mat_x_m1, m1_beta, m2_beta, sig2_inv_vb)
  
  lambda_vb <- update_lambda_vb_(lambda, sum(gam_vb))
  nu_vb <- update_nu_vb_(nu, m2_beta, tau_vb)
  
  log_tau_vb <- update_log_tau_vb_(eta_vb, kappa_vb)
  log_sig2_inv_vb <- update_log_sig2_inv_vb_(lambda_vb, nu_vb)
  
  elbo_A <- e_y_(n, kappa, kappa_vb, log_tau_vb, m2_beta, sig2_inv_vb, tau_vb)
  
  elbo_B <- e_beta_gamma_dual_info_(list_V, gam_vb, log_sig2_inv_vb, log_tau_vb,
                                    mat_v_mu, mu_c_vb, m2_beta,
                                    sig2_beta_vb, sig2_c_vb, sig2_rho_vb,
                                    sig2_theta_vb, sig2_inv_vb, tau_vb, zeta_vb, bool_blocks = TRUE)
  
  elbo_C <- e_theta_(m0, mu_theta_vb, S0_inv, sig2_theta_vb, vec_fac_st,
                     vec_sum_log_det_theta)
  
  elbo_D <- e_rho_(mu_rho_vb, n0, sig2_rho_vb, T0_inv, vec_sum_log_det_rho)
  
  elbo_E <- sum(sapply(1:n_bl, function(bl) e_c_zeta_(log_om_vb[[bl]], log_1_min_om_vb[[bl]], 
                                                      mu_c_vb[[bl]], s2[bl], sig2_c_vb[bl], zeta_vb[[bl]])))
  
  elbo_F <- e_tau_(eta, eta_vb, kappa, kappa_vb, log_tau_vb, tau_vb)
  
  elbo_G <- e_sig2_inv_(lambda, lambda_vb, log_sig2_inv_vb, nu, nu_vb, sig2_inv_vb)
  
  elbo_A + elbo_B + elbo_C + elbo_D + elbo_E + elbo_F + elbo_G 
  
}

