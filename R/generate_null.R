#' @export
generate_null <- function(n_perm, Y, X, p0_av, Z, list_hyper, list_init,
                          list_blocks, user_seed, tol, maxit, batch, verbose,
                          results_dir, n_cpus) {

  if (!is.null(user_seed)){
    RNGkind("L'Ecuyer-CMRG") # ensure reproducibility when using mclapply
    set.seed(user_seed)
    if (verbose) cat(paste("Seed set to ", user_seed, ". \n", sep = ""))
  }

  n <- nrow(Y)

  permute <- function(i) {
    ind_perm <- sample(1:n) # random permutation
    rownames(Y) <- NULL

    # user_seed must be NULL here otherwise always the same permutation
    res_perm <- locus(Y = Y[ind_perm, ], X = X, p0_av = p0_av, Z = Z,
                      list_hyper = list_hyper, list_init = list_init,
                      list_cv = NULL, list_blocks = list_blocks, user_seed = NULL,
                      tol = tol, maxit = maxit, batch = batch, verbose = verbose)

    om_vb <- res_perm$om_vb
    gam_vb <- res_perm$gam_vb
    rm(res_perm)

    if (is.null(results_dir)) {
      
      create_named_list_(ind_perm, gam_vb, om_vb)
      
    } else {
      
      save(ind_perm, gam_vb, om_vb,
           file = paste(results_dir, "vb_real_data_", i, ".RData", sep=""))
      rm(om_vb)
      rm(gam_vb)
      NULL
      
    }  
  }

  parallel::mclapply(1:n_perm, function(i) permute(i), mc.cores = n_cpus)

}
