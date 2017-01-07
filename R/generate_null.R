#' @export
generate_null <- function(n_perm, Y, X, p0_av, Z, list_hyper, list_init, user_seed,
                          tol, maxit, batch, verbose, rel_results_dir, n_cpus) {

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
    res_perm <- locus(Y[ind_perm, ], X, p0_av, Z, list_hyper, list_init,
                       list_cv = NULL, user_seed = NULL, tol, maxit, batch, verbose)

    om_vb <- res_perm$om_vb
    gam_vb <- res_perm$gam_vb
    rm(res_perm)

    save(ind_perm, gam_vb, om_vb,
         file = paste(rel_results_dir, "vb_real_data_", i, ".RData", sep=""))
    rm(om_vb)
    rm(gam_vb)
    NULL
  }

  parallel::mclapply(1:n_perm, function(i) permute(i), mc.cores = n_cpus)

}
