#' @export
generate_snps <- function(n, p, cor_type, vec_rho, vec_maf = NULL, n_cpus = 1) {

  check_structure_(n, "vector", "numeric", 1)
  check_natural_(n)

  check_structure_(p, "vector", "numeric", 1)
  check_natural_(p)

  stopifnot(cor_type %in% c("autocorrelated", "equicorrelated", "none"))

  if(is.null(vec_maf)) {
    vec_maf <- runif(p, 0.05, 0.5)
  } else {
    check_structure_(vec_maf, "vector", "numeric", p)
    check_zero_one_(vec_maf)
  }


  if (cor_type == "none") {

    snps <- sapply(vec_maf, function(maf) rbinom(n, 2, maf)) # Hardy-Weinberg equilibrium

  } else {

    if(length(vec_rho) > p)
      stop(paste("Provided number of blocks of correlated SNPs, length(vec_rho), ",
                 "must be smaller than the number of SNPs, p: ",
                 length(vec_rho), " > ", p, sep = ""))

    check_structure_(vec_rho, "vector", "numeric")
    if(cor_type == "equicorrelated") check_zero_one_(vec_rho)
    else check_zero_one_(abs(vec_rho))

    check_structure_(n_cpus, "vector", "numeric", 1)
    check_natural_(n_cpus)

    n_bl <- length(vec_rho)
    ind_bl <- make_chunks_(1:p, n_bl)

    snps <- parallel::mclapply(1:n_bl, function(bl) {

      p_bl <- length(ind_bl[[bl]])
      rho_bl <- vec_rho[[bl]]
      R <- matrix(NA, nrow = p_bl, ncol = p_bl)

      if (cor_type == "autocorrelated") {
        for( i in 1:p_bl ){
          for( j in i:p_bl ){
            R[i,j] <- rho_bl^abs(i - j)
            R[j,i] <- R[i,j]
          }
        }
      } else { # equicorrelated
        R[] <- rho_bl
      }

      diag(R) <- 1

      L <- t(chol(R))
      tZ <- matrix(rnorm(n * p_bl), nrow = p_bl, ncol = n)
      X <- t(L %*% tZ)

      snps <- matrix(1, nrow = n, ncol = p_bl)

      for(j in 1:p_bl) {
        maf <- vec_maf[ind_bl[[bl]]][j]
        snps[X[,j] < qnorm((1 - maf)^2), j] <- 0
        snps[X[,j] > qnorm(1 - maf^2), j] <- 2
      }
      snps
    }, mc.cores = n_cpus)

    snps <- do.call(cbind, snps)
  }

  rownames(snps) <- paste("ind_", 1:n, sep = "")
  colnames(snps) <- paste("snp_", 1:p, sep = "")

  vec_maf <- apply(snps, 2, mean) / 2 # empirical maf
  names(vec_maf) <- colnames(snps)

  list_snps <- create_named_list_(snps, vec_maf)
  class(list_snps) <- "sim_snps"
  list_snps
}

#' @export
replicate_real_snps <- function(n, real_snps, bl_lgth, p = NULL, maf_thres = NULL,
                                n_cpus = 1) {

  check_structure_(real_snps, "matrix", "numeric")

  check_structure_(maf_thres, "vector", "numeric", 1, null_ok = T)
  if(!is.null(maf_thres)) check_zero_one_(maf_thres)

  check_structure_(n_cpus, "vector", "numeric", 1)
  check_natural_(n_cpus)

  # n can be larger than the number of available observations for the real snps
  check_structure_(n, "vector", "numeric", 1)
  check_natural_(n)

  check_structure_(p, "vector", "numeric", 1, null_ok = T)
  if(!is.null(p)) check_natural_(p)

  if (!is.null(p)) {
    p_av <- ncol(real_snps)
    if(p > p_av)
      stop(paste("Provided n greater than number of snps available: ",
                 p, " > ", p_av, sep = ""))
    real_snps <- real_snps[, 1:p]
  } else {
    p <- ncol(real_snps)
  }

  check_structure_(bl_lgth, "vector", "numeric", 1)
  check_natural_(bl_lgth)

  if (bl_lgth > p)
    stop(paste("Provided block length must be smaller than the number of SNPs available: ",
               bl_lgth, " > ", p, sep = ""))

  vec_real_maf <- apply(real_snps, 2, mean) / 2

  n_bl <- floor(p / bl_lgth)
  ind_bl <- make_chunks_(1:p, n_bl)

  snps <- parallel::mclapply(1:n_bl, function(bl) {
    R <- cor(real_snps[, ind_bl[[bl]]])
    R <- Matrix::nearPD(R, corr = T, do2eigen = T)$mat
    p_bl <- ncol(R)
    L <- t(chol(R))
    tZ <- matrix(rnorm(n * p_bl), nrow = p_bl, ncol = n)
    X <- t(L %*% tZ) # Gaussian variables

    snps <- matrix(1, nrow = n, ncol = p_bl)
    for(j in 1:p_bl) {
      maf <- vec_real_maf[ind_bl[[bl]]][j]
      snps[X[,j] < qnorm((1 - maf)^2), j] <- 0
      snps[X[,j] > qnorm(1 - maf^2), j] <- 2
    }
    snps
  }, mc.cores = n_cpus)

  snps <- do.call(cbind, snps)
  rownames(snps) <- paste("ind_", 1:n, sep = "")
  colnames(snps) <- paste("snp_", 1:p, sep = "")

  vec_maf <- apply(snps, 2, mean) / 2
  names(vec_maf) <- colnames(snps)

  if (!is.null(maf_thres)) {
    ind_rare <- vec_maf < maf_thres
    snps <- snps[, !ind_rare]
    vec_maf <- vec_maf[!ind_rare]
  }

  list_snps <- create_named_list_(snps, vec_maf)
  class(list_snps) <- "sim_snps"
  list_snps
}

#' @export
generate_phenos <- function(n, d, cor_type, vec_rho, vec_var_err, n_cpus = 1) {

  check_structure_(n, "vector", "numeric", 1)
  check_natural_(n)

  check_structure_(d, "vector", "numeric", 1)
  check_natural_(d)

  stopifnot(cor_type %in% c("autocorrelated", "equicorrelated", "none"))

  check_structure_(vec_var_err, "vector", "numeric", d)
  check_positive_(vec_var_err)

  if (cor_type == "none") {

    phenos <- sapply(vec_var_err, function(var_err) rnorm(n, 0, sqrt(var_err))) # Hardy-Weinberg equilibrium
    ind_bl <- NULL

  } else {

    if(length(vec_rho) > d)
      stop(paste("Provided number of blocks of correlated phenotypes, length(vec_rho), ",
                 "must be smaller than the number of phenotypes, d: ",
                 length(vec_rho), " > ", d, sep = ""))

    check_structure_(vec_rho, "vector", "numeric")
    if(cor_type == "equicorrelated") check_zero_one_(vec_rho)
    else check_zero_one_(abs(vec_rho))

    check_structure_(n_cpus, "vector", "numeric", 1)
    check_natural_(n_cpus)

    n_bl <- length(vec_rho)
    ind_bl <- make_chunks_(1:d, n_bl)

    phenos <- parallel::mclapply(1:n_bl, function(bl) {

      d_bl <- length(ind_bl[[bl]])
      rho_bl <- vec_rho[[bl]]
      R <- matrix(NA, nrow = d_bl, ncol = d_bl)

      if (cor_type == "autocorrelated") {
        for( i in 1:d_bl ){
          for( j in i:d_bl ){
            R[i,j] <- rho_bl^abs(i - j)
            R[j,i] <- R[i,j]
          }
        }
      } else { # equicorrelated
        R[] <- rho_bl
      }
      diag(R) <- 1

      L <- t(chol(R))
      tZ <- matrix(rnorm(n * d_bl), nrow = d_bl, ncol = n)
      as.matrix(t(L %*% tZ))
    }, mc.cores = n_cpus)

    phenos <- do.call(cbind, phenos)
  }

  rownames(phenos) <- paste("ind_", 1:n, sep = "")
  colnames(phenos) <- paste("pheno_", 1:d, sep = "")

  vec_var_err <- apply(phenos, 2, var) # empirical error variance
  names(vec_var_err) <- colnames(phenos)

  list_phenos <- create_named_list_(phenos, vec_var_err, ind_bl)
  class(list_phenos) <- "sim_phenos"
  list_phenos
}

#' @export
replicate_real_phenos <- function(n, real_phenos, bl_lgth = NULL, d = NULL, n_cpus = 1) {

  check_structure_(real_phenos, "matrix", "numeric")

  check_structure_(n_cpus, "vector", "numeric", 1)
  check_natural_(n_cpus)

  check_structure_(n, "vector", "numeric", 1) # can be larger than the number of available observations in real_pheno!
  check_natural_(n)

  check_structure_(d, "vector", "numeric", 1, null_ok = T)
  if(!is.null(d)) check_natural_(d)

  if (!is.null(d)) {
    d_av <- ncol(real_phenos)
    if(d > d_av)
      stop(paste("Provided n greater than number of phenotypes available: ",
                 d, " > ", d_av, sep = ""))
    real_phenos <- real_phenos[, 1:d]
  } else {
    d <- ncol(real_phenos)
  }

  if (is.null(bl_lgth)) {  # bl_lgth = NULL, means one block

    bl_lgth <- d

  } else {

    check_structure_(bl_lgth, "vector", "numeric", 1)
    check_natural_(bl_lgth)

    if (bl_lgth > d)
      stop(paste("Provided block length must be smaller or equal to the number ",
                 "of phenotypes available: ", bl_lgth, " > ", d, sep = ""))
  }

  n_bl <- floor(d / bl_lgth)
  ind_bl <- make_chunks_(1:d, n_bl)

  phenos <- parallel::mclapply(1:n_bl, function(bl) {
    R <- cor(real_phenos[, ind_bl[[bl]]])
    R <- Matrix::nearPD(R, corr = T, do2eigen = T)$mat
    d_bl <- ncol(R)
    L <- t(chol(R))
    tZ <- matrix(rnorm(n * d_bl), nrow = d_bl, ncol = n)
    as.matrix(t(L %*% tZ)) # Gaussian variables
  }, mc.cores = n_cpus)

  phenos <- do.call(cbind, phenos)

  rownames(phenos) <- paste("ind_", 1:n, sep = "")
  colnames(phenos) <- paste("pheno_", 1:d, sep = "")

  vec_var_err <- apply(phenos, 2, var) # empirical error variance
  names(vec_var_err) <- colnames(phenos)

  ind_bl <- NULL # block chuncks do not necessarily represent blocks of correlated phenotypes.

  list_phenos <- create_named_list_(phenos, vec_var_err, ind_bl)
  class(list_phenos) <- "sim_phenos"
  list_phenos
}


set_pattern_ <- function(d, p, ind_d0, ind_p0, vec_prob_sh, chunks_ph){

  if (is.null(chunks_ph)) { # no imposed correlation block structure (either indep
                            # or correlation from real phenotypes). creates artificial chunks.
    n_chunks_ph <- length(vec_prob_sh)
    chunks_ph <- make_chunks_(1:d, n_chunks_ph)
  } else {
    n_chunks_ph <- length(chunks_ph)
  }

  ind_d0 <- unique(ind_d0)
  if(!all(ind_d0 %in% 1:d))
    stop("All indices provided in ind_d0 must be integers between 1 and d.")

  ind_p0 <- unique(ind_p0)
  if(!all(ind_p0 %in% 1:p))
    stop("All indices provided in ind_p0 must be integers between 1 and p.")

  check_structure_(vec_prob_sh, "vector", "numeric")
  check_zero_one_(vec_prob_sh)

  pat <- matrix(F, nrow=p, ncol=d)

  for(ind_j in ind_p0) {

    # random permutation, so that two different SNPs can be associated with a given
    # block with different probabilities
    vec_prob_sh_perm <- sample(vec_prob_sh, n_chunks_ph, replace = T)
    for(ch in 1:n_chunks_ph) {
      ind_ch_d0 <- intersect(ind_d0, chunks_ph[[ch]])
      pat[ind_j, ind_ch_d0] <- sample(c(T, F),
                                      length(ind_ch_d0), replace=T,
                                      prob=c(vec_prob_sh_perm[ch], 1 - vec_prob_sh_perm[ch]))
    }

    if (all(!pat[ind_j, ind_d0]))
      pat[ind_j, ind_d0][sample(1:length(ind_d0), 1)] <- T # each active covariate must be
    # associated with at least one response
  }

  for(ind_k in ind_d0) {
    if (all(!pat[ind_p0, ind_k]))
      # each active covariate must be associated with at least one response
      pat[ind_p0, ind_k][sample(1:length(ind_p0), 1)] <- T
  }
  pat
}


generate_eff_sizes_ <- function(d, p, ind_d0, ind_p0, vec_prob_sh, vec_maf,
                                pve_per_snp, max_tot_pve, vec_var_err, chunks_ph) {

  # pve_per_snp average variance explained per snp
  check_structure_(ind_d0, "vector", "numeric", null_ok = T)
  check_structure_(ind_p0, "vector", "numeric", null_ok = T)

  d0 <- length(ind_d0)
  p0 <- length(ind_p0)

  if (xor(d0 < 1, p0 < 1))
    stop("ind_d0 and ind_p0 must either have both length 0 (no association) or both length >= 1.")

  if (d0 < 1 & p0 < 1) {

    pat <- matrix(F, nrow = p, ncol = d)
    beta <- matrix(0.0, nrow = p, ncol = d)

  } else {

    bool_cst <- vec_maf[ind_p0] %in% c(0.0, 1.0) # there may be remaining constant snps! rm them after.
    if (any(bool_cst)) {
      ind_p0 <- ind_p0[!bool_cst]
      if (length(ind_p0) == 0)
        stop(paste("SNP(s) number ", ind_p0[bool_cst], " constant. Effect(s) on ",
                   "the responses removed.\n No remaining ``active'' snps, change ",
                   "ind_p0 (now empty).\n", sep = ""))
      warning(paste("SNP(s) number ", ind_p0[bool_cst], " constant. Effect(s) on",
                    "the responses removed.", sep = ""))
    }

    pat <- set_pattern_(d, p, ind_d0, ind_p0, vec_prob_sh, chunks_ph)

    check_structure_(vec_maf, "vector", "numeric", p)
    check_zero_one_(vec_maf)

    check_structure_(vec_var_err, "vector", "numeric", d)
    check_positive_(vec_var_err)

    max_per_resp <- max(colSums(pat))
    eps <- .Machine$double.eps^0.75
    if (is.null(pve_per_snp)) {
      # sets pve_per_snp to the max possible so that the tot_pve for all responses are below 1.
      if (is.null(max_tot_pve)) max_tot_pve <- 1 - eps

      pve_per_snp <- max_tot_pve / max_per_resp
    } else {
      check_structure_(pve_per_snp, "vector", "numeric", 1)
      check_zero_one_(pve_per_snp)

      if (max_per_resp * pve_per_snp > 1 - eps)
        stop(paste("Provided average proportion of variance explained per SNP too ",
                   "high, would lead to a total genetic variance explained above ",
                   "100% for at least one response. \n Setting pve_per_snp < 1 / length(ind_p0) ",
                   "will work for any pattern. \n", sep = ""))
    }

    beta <- matrix(0.0, nrow = p, ncol = d)

    beta[, ind_d0] <- sapply(ind_d0, function(k) {

      p0_k <- sum(pat[,k])
      vec_pve_per_snp <- rbeta(p0_k, shape1 = 2, shape2 = 5) # positively skewed Beta distribution,
                                                             # to give more weight to smaller effect sizes
      vec_pve_per_snp <- vec_pve_per_snp / sum(vec_pve_per_snp) * pve_per_snp * p0_k

      tot_var_expl <- pve_per_snp * p0_k * vec_var_err[k] / (1 - pve_per_snp * p0_k)

      vec_maf_act <- vec_maf[pat[,k]]
      vec_var_act <- 2 * vec_maf_act * (1 - vec_maf_act)

      beta_k <- rep(0.0, p)
      beta_k[pat[,k]] <- sqrt((tot_var_expl + vec_var_err[k]) * vec_pve_per_snp / vec_var_act)

      # switches signs with probabilty 0.5
      beta_k[pat[,k]] <- sample(c(1, -1), p0_k, replace = T) * beta_k[pat[,k]]

      beta_k
    })

  }

  create_named_list_(beta, pat, pve_per_snp)

}

#' @export
generate_dependence <- function(list_snps, list_phenos, ind_d0, ind_p0, vec_prob_sh,
                                pve_per_snp = NULL, max_tot_pve = NULL) {

  if (!is.null(pve_per_snp) & !is.null(max_tot_pve))
    stop("Either pve_per_snp or max_tot_pve must be NULL.")

  if (is.null(pve_per_snp) & is.null(max_tot_pve))
    warning(paste("As both pve_per_snp or max_tot_pve were provided as NULL, the ",
                  "pve per SNP was to its maximum value so that the total pve for ",
                  "the responses are all below 1.", sep = ""))

  if (class(list_snps) != "sim_snps")
    stop(paste("The provided list_snps must be an object of class ``sim_snps''. \n",
               "*** You must either use the function generate_snps to simulate snps ",
               "under Hardy-Weinberg equilibrium or the function replicate_real_snps ",
               "to simulate SNPs from real SNP data, by replicating their minor ",
               " allele frequencies and linkage desequilibrium structure. ***",
               sep=""))
  list2env(list_snps, envir=environment())
  rm(list_snps)

  n <- nrow(snps)
  p <- ncol(snps)

  if (class(list_phenos) != "sim_phenos")
    stop(paste("The provided list_phenos must be an object of class ``sim_phenos''. \n",
               "*** You must either use the function generate_phenos to simulate ",
               "phenotypes from (possibly correlated) gaussian variables or the ",
               "function replicate_real_phenos to simulate phenotypes from real ",
               "phenotypic data, by replicating their correlation structure. ***",
               sep=""))
  list2env(list_phenos, envir=environment())
  rm(list_phenos)

  if(n != nrow(phenos))
    stop("The number of observations used for list_snps and for list_phenos does not match.")

  d <- ncol(phenos)

  bool_cst <- is.nan(colSums(snps[, ind_p0]))
  if (any(bool_cst)) {
    ind_p0 <- ind_p0[!bool_cst]
    if (length(ind_p0) == 0)
      stop(paste("SNP(s) number ", ind_p0[bool_cst], " constant. Effect(s) on the ",
                 "responses removed.\n No remaining ``active'' snps, change ",
                 "ind_p0 (now empty).", sep = ""))
  }

  list_eff <- generate_eff_sizes_(d, p, ind_d0, ind_p0, vec_prob_sh, vec_maf,
                                  pve_per_snp, max_tot_pve, vec_var_err,
                                  chunks_ph = ind_bl)

  list2env(list_eff, envir=environment())
  rm(list_eff)

  phenos <- phenos + snps %*% beta

  create_named_list_(phenos, snps, beta, pat, pve_per_snp)
}
