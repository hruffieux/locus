locus_dual_vbem_core_ <- function(Y, X, list_hyper, gam_vb, mu_beta_vb,
                                       sig2_beta_vb, tau_vb, list_struct, tol, maxit,
                                       anneal, verbose) {

  maxit_em <- ceiling(maxit / 5)
  
  converged_em <- FALSE
  it_em <- 0
  tol_em <- 10 * tol # TODO: rather monitor the elbo across VB runs
  s02_min <- 1e-6
  list_hyper$om_vb <- rep(1/20, r) # prior proportion of hotspots
  
  vb <- create_named_list_(gam_vb, mu_beta_vb, sig2_beta_vb, tau_vb)
  
  
  while ((!converged_em) & list_hyper$s02 > s02_min & (it_em < maxit_em)) {
    
    it_em <- it_em + 1
    
    s02_old <- list_hyper$s02
    om_old <- list_hyper$om_vb
    
    if (verbose)
      cat("---------- VB updates ----------\n")
    
    vb <- locus_dual_pleio_core_(Y, X, list_hyper, vb$gam_vb, vb$mu_beta_vb,
                           vb$sig2_beta_vb, vb$tau_vb, list_struct, tol,
                           maxit, anneal, verbose = FALSE, full_output = TRUE, bool_eb = TRUE)
    
    if (verbose)
      cat("--- EM hyperparameter updates ---\n")
    
    list_hyper$om_vb <- vb$zeta_vb
    list_hyper$s02 <- sum(vb$zeta_vb * (vb$sig2_theta_vb + vb$mu_theta_vb^2)) / sum(vb$zeta_vb) ## check!
    
    if (verbose) {
      cat(paste0("EM iteration ", it_em, ". \n New value for hyperparameter s02 : ", format(list_hyper$s02, digits = 4), 
                 ". \n New values for hyperparameter omega : \n"))
      print(summary(list_hyper$om_vb))
      cat("\n\n")
    }
    
    converged_em <- (abs(list_hyper$s02 - s02_old) / s02_old < tol_em) && (max(abs(list_hyper$om_vb - om_old) / om_old) < tol_em)
    
  }
  
  if (verbose) {
    if (converged_em) {
      cat(paste0("Convergence of the EM hyperparameter optimization run obtained after ", format(it_em), " EM iterations. \n\n"))
    } else if (list_hyper$s02 <= s02_min) {
      cat(paste0("EM hyperparameter optimization run obtained after ", format(list_hyper$s02), " getting below 1e-5. \n\n"))
    } else {
      warning("Maximal number of EM iterations reached before convergence. Exit EM run. \n\n")
    }
    
    cat("======= Final VB run =======\n") 
    cat(paste0("Empirical-Bayes hyperparameters, s02 : ", format(list_hyper$s02, digits = 4), ", omega :\n"))
    print(summary(list_hyper$om_vb))
    cat("\n\n")
  }
  
  
  out <- locus_dual_pleio_core_(Y, X,  list_hyper, vb$gam_vb, vb$mu_beta_vb,
                               vb$sig2_beta_vb, vb$tau_vb, list_struct, tol,
                               maxit, anneal, verbose, bool_eb = TRUE)
  
  out$s02 <- list_hyper$s02
  out
  
}

