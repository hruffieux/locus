locus_dual_vbem_core_ <- function(Y, X, list_hyper, gam_vb, mu_beta_vb,
                                  sig2_beta_vb, tau_vb, bool_blocks, tol, maxit,
                                  anneal, verbose) {

  maxit_em <- ceiling(maxit / 5)
  converged_em <- FALSE
  it_em <- 0
  tol_em <- 1e-3
  s02_min <- 1e-6
  lb_old <- -Inf
  list_hyper$om_vb <- rep(1/20, ncol(X)) # initial proportion of hotspots, the initial s02 is user specified...
  
  vb <- create_named_list_(gam_vb, mu_beta_vb, sig2_beta_vb, tau_vb)
  
  
  while ((!converged_em) & list_hyper$s02 > s02_min & (it_em < maxit_em)) {
    
    it_em <- it_em + 1
    
    
    if (verbose)
      cat("---------- VB updates ----------\n")
    
    vb <- locus_dual_pleio_core_(Y, X, list_hyper, vb$gam_vb, vb$mu_beta_vb,
                           vb$sig2_beta_vb, vb$tau_vb, eb = TRUE, tol,
                           maxit, anneal, verbose = FALSE, full_output = TRUE)
    
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
    
    converged_em <- (abs(vb$lb_opt - lb_old) < tol_em)
    
    lb_old <- vb$lb_opt
    
  }
  
  
  if (bool_blocks) {
    
    out <- list("s02" = rep(list_hyper$s02, ncol(X)), "om" = list_hyper$om_vb)
    
  } else {
    
    if (verbose) {
      if (converged_em) {
        cat(paste0("Convergence of the EM hyperparameter optimization run obtained after ", format(it_em), " EM iterations. \n\n"))
      } else if (list_hyper$s02 <= s02_min) {
        cat(paste0("EM hyperparameter optimization run stopped after s02 getting below ", s02_min, ". \n\n"))
      } else {
        warning("Maximal number of EM iterations reached before convergence. Exit EM run. \n\n")
      }
      
      cat("======= Final VB run =======\n") 
      cat(paste0("Empirical-Bayes hyperparameters, s02 : ", format(list_hyper$s02, digits = 4), ", omega :\n"))
      print(summary(list_hyper$om_vb))
      cat("\n\n")
    }
    
    
    out <- locus_dual_pleio_core_(Y, X,  list_hyper, vb$gam_vb, vb$mu_beta_vb,
                                  vb$sig2_beta_vb, vb$tau_vb, eb = TRUE, tol,
                                  maxit, anneal, verbose)
    
    out$s02 <- list_hyper$s02
    
  }

  out
  
}

