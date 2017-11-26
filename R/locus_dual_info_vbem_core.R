locus_dual_info_vbem_core_ <- function(Y, X, V, list_hyper, gam_vb, mu_beta_vb,
                                       sig2_beta_vb, tau_vb, list_struct, tol, maxit,
                                       anneal, verbose) {
  
  r <- ncol(V)
  
  maxit_em <- ceiling(maxit / 5)
  
  converged_em <- FALSE
  it_em <- 0
  tol_em <- tol
  
  vb <- create_named_list_(gam_vb, mu_beta_vb, sig2_beta_vb, tau_vb)
  
  while ((!converged_em) & (it_em < maxit_em)) {
    
    it_em <- it_em + 1
    s2_old <- list_hyper$s2
    
  
    if (it_em == 1) list_hyper$om_vb <- rep(1/5, r) # prior proportion of active annotations # 1 / 5, of length 1 doesn't work well of length 1 and not r as before!
    
    om_old <- list_hyper$om_vb
    
    if (verbose)
      cat("---------- VB updates ----------\n")
    
    vb <- locus_dual_info_core_(Y, X, V, list_hyper, vb$gam_vb, vb$mu_beta_vb,
                                vb$sig2_beta_vb, vb$tau_vb, list_struct, tol,
                                maxit, anneal, verbose = FALSE, full_output = TRUE, bool_eb = TRUE)
    
    if (verbose)
      cat("--- EM hyperparameter updates ---\n")
    
    list_hyper$om_vb <- vb$zeta_vb
    list_hyper$s2 <- sum(vb$zeta_vb * (vb$sig2_c_vb + vb$mu_c_vb^2)) / sum(vb$zeta_vb)
  
    if (verbose) {
      cat(paste0("EM iteration ", it_em, ". \n New value for hyperparameter s2 : ", format(list_hyper$s2, digits = 4), 
                 ". \n New values for hyperparameter omega : \n"))
      print(summary(list_hyper$om_vb))
      cat("\n\n")
    }
     
    converged_em <- (abs(list_hyper$s2 - s2_old) / s2_old < tol_em) && (max(abs(list_hyper$om_vb - om_old) / om_old) < tol_em)
    
  }
  
  if (verbose) {
    if (converged_em) {
      cat(paste0("Convergence of the EM hyperparameter optimization run obtained after ", format(it_em), " EM iterations. \n"))
    } else {
      warning("Maximal number of EM iterations reached before convergence. Exit EM run. \n")
    }
    
    cat("======= Final VB run =======\n") 
    cat(paste0("Empirical-Bayes hyperparameters, s2 : ", format(list_hyper$s2, digits = 4), ", omega :\n"))
    print(summary(list_hyper$om_vb))
    cat("\n\n")
  }


  out <- locus_dual_info_core_(Y, X, V, list_hyper, vb$gam_vb, vb$mu_beta_vb,
                               vb$sig2_beta_vb, vb$tau_vb, list_struct, tol,
                               maxit, anneal, verbose, bool_eb = TRUE)

  out$s2 <- list_hyper$s2
  out
  
}

