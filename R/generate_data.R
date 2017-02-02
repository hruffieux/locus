#' Generate SNP data with prespecified spatial correlation structure.
#'
#' This function generates SNPs under Hardy-Weinberg equilibrium with specific
#' block correlation structure and minor allele frequencies.
#'
#' @param n Number of observations.
#' @param p Number of SNPs.
#' @param cor_type String describing the type of dependence structure. The SNPs
#'   can \code{autocorrelated}, \code{equicorrelated}. Set to \code{NULL} for
#'   independent SNPs.
#' @param vec_rho Vector of correlation coefficients. Its length determines the
#'   number of blocks of correlated SNPs. Must be smaller than p. Set to
#'   \code{NULL} if independent SNPs.
#' @param vec_maf Vector of size p containing the reference minor allele
#'   frequencies used to draw the SNPs. If \code{NULL}, the minor allele
#'   frequencies drawn uniformly at random between 0.05 and 0.5.
#' @param n_cpus Number of CPUs used when simulating correlated SNP blocks.
#'   Ignored if independent SNPs. Set to 1 for serial execution.
#' @param user_seed Seed set for reproducibility. Default is \code{NULL}, no
#'   seed set.
#'
#' @return An object of class "\code{sim_snps}".
#'  \item{snps}{Matrix containing the generated SNP data.}
#'  \item{vec_maf}{Vector containing the SNP sample minor allele frequencies.}
#'
#' @examples
#' user_seed <- 123; set.seed(user_seed)
#' n <- 500; p <- 10000
#' cor_type <- "autocorrelated"; vec_rho <- runif(100, min = 0.25, max = 0.95)
#' list_snps <- generate_snps(n, p, cor_type, vec_rho, n_cpus = 2,
#'                            user_seed = user_seed)
#'
#' @seealso \code{\link{replicate_real_snps}}, \code{\link{generate_phenos}},
#'   \code{\link{replicate_real_phenos}}, \code{\link{generate_dependence}}
#'
#' @export
#'
generate_snps <- function(n, p, cor_type = NULL, vec_rho = NULL, vec_maf = NULL,
                          n_cpus = 1, user_seed = NULL) {

  check_structure_(user_seed, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(user_seed)){
    RNGkind("L'Ecuyer-CMRG") # ensure reproducibility when using mclapply
    set.seed(user_seed)
  }

  check_structure_(n, "vector", "numeric", 1)
  check_natural_(n)

  check_structure_(p, "vector", "numeric", 1)
  check_natural_(p)

  if (!is.null(cor_type))
    stopifnot(cor_type %in% c("autocorrelated", "equicorrelated"))

  if(is.null(vec_maf)) {
    vec_maf <- runif(p, 0.05, 0.5)
  } else {
    check_structure_(vec_maf, "vector", "numeric", p)
    check_zero_one_(vec_maf)
  }

  if (is.null(cor_type)) {

    if (n_cpus > 1)
      warning("n_cpus is ignored when the SNPs are generated independently of one another.")

    snps <- sapply(vec_maf, function(maf) rbinom(n, 2, maf)) # Hardy-Weinberg equilibrium

  } else {

    check_structure_(vec_rho, "vector", "numeric")
    if(cor_type == "equicorrelated") check_zero_one_(vec_rho)
    else check_zero_one_(abs(vec_rho))

    if(length(vec_rho) > p)
      stop(paste("Provided number of blocks of correlated SNPs, length(vec_rho), ",
                 "must be smaller than the number of SNPs, p: ",
                 length(vec_rho), " > ", p, sep = ""))

    check_structure_(n_cpus, "vector", "numeric", 1)
    check_natural_(n_cpus)

    if (n_cpus > 1) {
      n_cpus_avail <- parallel::detectCores()
      if (n_cpus > n_cpus_avail) {
        n_cpus <- n_cpus_avail
        warning(paste("The number of CPUs specified exceeds the number of CPUs ",
                      "available on the machine. The latter has been used instead.", sep=""))
      }
    }

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
      tZ <- matrix(rnorm(n * p_bl), ncol = n)
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


#' Generate SNPs emulating real SNP data at hand.
#'
#' This function simulates SNPs from real SNP data based on their sample minor
#' allele frequencies and correlation structure.
#'
#' @param n Number of observations.
#' @param real_snps Matrix of real SNPs (rows observations, columns SNP
#'   variables), without missing values. The entries must be 0, 1 or 2.
#' @param bl_lgth Number of variables per block for reproducing the dependence
#'   structure of real SNPs. Must be between 2 and p. Must be small enough
#'   (e.g. 1000) for tractability reasons.
#' @param p Number of SNPs. If \code{NULL}, the total number of SNPs
#'   available in the real_snps matrix is used (default).
#' @param maf_thres Lower bound for sample minor allele frequencies. Simulated
#'   SNPs with lower sample minor allele frequencies are excluded (as a result
#'   the number of simulated SNPs can be smaller than p). Default is \code{NULL}
#'   for no exclusion.
#' @param n_cpus Number of CPUs used when simulating SNP by blocks. Set to 1 for
#'   serial execution.
#' @param user_seed Seed set for reproducibility. Default is \code{NULL}, no
#'   seed set.
#'
#' @return An object of class "\code{sim_snps}".
#'  \item{snps}{Matrix containing the generated SNP data.}
#'  \item{vec_maf}{Vector containing the SNP sample minor allele frequencies.}
#'
#' @examples
#' user_seed <- 123; set.seed(user_seed)
#' n <- 500; p <- 7500
#' cor_type <- "autocorrelated"; vec_rho <- runif(100, min = 0.25, max = 0.95)
#' list_fake_real_snps <- generate_snps(n, p, cor_type, vec_rho, n_cpus = 2,
#'                            user_seed = user_seed)
#' list_snps <- replicate_real_snps(n, list_fake_real_snps$snps, bl_lgth = 100,
#'                                  n_cpus = 2, user_seed = user_seed)
#'
#' @seealso \code{\link{generate_snps}}, \code{\link{generate_phenos}},
#'   \code{\link{replicate_real_phenos}}, \code{\link{generate_dependence}}
#'
#' @export
#'
replicate_real_snps <- function(n, real_snps, bl_lgth, p = NULL, maf_thres = NULL,
                                n_cpus = 1, user_seed = NULL) {


  check_structure_(user_seed, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(user_seed)){
    RNGkind("L'Ecuyer-CMRG") # ensure reproducibility when using mclapply
    set.seed(user_seed)
  }

  check_structure_(real_snps, "matrix", "numeric")

  check_structure_(maf_thres, "vector", "numeric", 1, null_ok = TRUE)
  if(!is.null(maf_thres)) check_zero_one_(maf_thres)

  check_structure_(n_cpus, "vector", "numeric", 1)
  check_natural_(n_cpus)

  if (n_cpus > 1) {
    n_cpus_avail <- parallel::detectCores()
    if (n_cpus > n_cpus_avail) {
      n_cpus <- n_cpus_avail
      warning(paste("The number of CPUs specified exceeds the number of CPUs ",
                    "available on the machine. The latter has been used instead.", sep=""))
    }
  }

  # n can be larger than the number of available observations for the real snps
  check_structure_(n, "vector", "numeric", 1)
  check_natural_(n)

  check_structure_(p, "vector", "numeric", 1, null_ok = TRUE)
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

  if (bl_lgth == 1) stop("Provided block length must be larger than 1")
  if (bl_lgth > p)
    stop(paste("Provided block length must be smaller than the number of SNPs available: ",
               bl_lgth, " > ", p, sep = ""))

  vec_real_maf <- apply(real_snps, 2, mean) / 2

  n_bl <- floor(p / bl_lgth)
  ind_bl <- make_chunks_(1:p, n_bl)

  snps <- parallel::mclapply(1:n_bl, function(bl) {
    R <- cor(real_snps[, ind_bl[[bl]]])
    R <- Matrix::nearPD(R, corr = TRUE, do2eigen = TRUE)$mat
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


#' Generate phenotypic data with prespecified block-wise correlation structure.
#'
#' This function generates Gaussian phenotypes with specific block correlation
#' structure.
#'
#' @param n Number of observations.
#' @param d Number of phenos.
#' @param var_err Vector of length 1 or d containing the variances of the
#'   Gaussian distributions used to draw the phenotytpes. If of length 1, the
#'   value is repeated d times.
#' @param cor_type String describing the type of dependence structure. The
#'   phenotypes can \code{autocorrelated}, \code{equicorrelated}. Set to
#'   \code{NULL} for independent phenotypes.
#' @param vec_rho Vector of correlation coefficients. Its length determines the
#'   number of blocks of correlated phenotypes. Must be smaller than d. Set to
#'   \code{NULL} if independent phenotypes.
#' @param n_cpus Number of CPUs used when simulating correlated phenotype blocks.
#'   Ignored if independent phenotypes. Set to 1 for serial execution.
#' @param user_seed Seed set for reproducibility. Default is \code{NULL}, no
#'   seed set.
#'
#' @return An object of class "\code{sim_phenos}".
#'  \item{phenos}{Matrix containing the generated phenotypic data.}
#'  \item{var_err}{Vector containing the sample phenotypic variances.}
#'  \item{ind_bl}{List of length given by the number of blocks, containing the
#'                indices of the phenotypes in each block. Is \code{NULL} if
#'                \code{cor_type} is \code{NULL}.}
#'
#' @examples
#' user_seed <- 123; set.seed(user_seed)
#' n <- 500; d <- 10000; var_err <- runif(d, min = 0.1, max = 0.4)
#' cor_type <- "equicorrelated"; vec_rho <- runif(100, min = 0.25, max = 0.95)
#' list_phenos <- generate_phenos(n, d, var_err, cor_type, vec_rho, n_cpus = 2,
#'                                user_seed = user_seed)
#'
#' @seealso \code{\link{generate_snps}}, \code{\link{replicate_real_snps}},
#'   \code{\link{replicate_real_phenos}}, \code{\link{generate_dependence}}
#'
#' @export
#'
#' @export
generate_phenos <- function(n, d, var_err, cor_type = NULL, vec_rho = NULL,
                            n_cpus = 1, user_seed = NULL) {


  check_structure_(user_seed, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(user_seed)){
    RNGkind("L'Ecuyer-CMRG") # ensure reproducibility when using mclapply
    set.seed(user_seed)
  }

  check_structure_(n, "vector", "numeric", 1)
  check_natural_(n)

  check_structure_(d, "vector", "numeric", 1)
  check_natural_(d)

  if (!is.null(cor_type))
    stopifnot(cor_type %in% c("autocorrelated", "equicorrelated"))

  check_structure_(var_err, "vector", "numeric", c(1, d))
  check_positive_(var_err)
  if (length(var_err) == 1) var_err <- rep(var_err, d)

  if (is.null(cor_type)) {

    if (n_cpus > 1)
      warning("n_cpus is ignored when the phenotypes are generated independently of one another.")

    phenos <- sapply(var_err, function(var_err) rnorm(n, 0, sqrt(var_err))) # Hardy-Weinberg equilibrium
    ind_bl <- NULL

  } else {

    check_structure_(vec_rho, "vector", "numeric")
    if(cor_type == "equicorrelated") check_zero_one_(vec_rho)
    else check_zero_one_(abs(vec_rho))

    if(length(vec_rho) > d)
      stop(paste("Provided number of blocks of correlated phenotypes, length(vec_rho), ",
                 "must be smaller than the number of phenotypes, d: ",
                 length(vec_rho), " > ", d, sep = ""))

    check_structure_(n_cpus, "vector", "numeric", 1)
    check_natural_(n_cpus)

    if (n_cpus > 1) {
      n_cpus_avail <- parallel::detectCores()
      if (n_cpus > n_cpus_avail) {
        n_cpus <- n_cpus_avail
        warning(paste("The number of CPUs specified exceeds the number of CPUs ",
                      "available on the machine. The latter has been used instead.", sep=""))
      }
    }

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
      tZ <- matrix(sapply(var_err[ind_bl[[bl]]], function(ve) rnorm(n, 0, sqrt(ve))),
                   ncol = n, byrow = TRUE)
      as.matrix(t(L %*% tZ))
    }, mc.cores = n_cpus)

    phenos <- do.call(cbind, phenos)
  }

  rownames(phenos) <- paste("ind_", 1:n, sep = "")
  colnames(phenos) <- paste("pheno_", 1:d, sep = "")

  var_err <- apply(phenos, 2, var) # empirical error variance
  names(var_err) <- colnames(phenos)

  list_phenos <- create_named_list_(phenos, var_err, ind_bl)
  class(list_phenos) <- "sim_phenos"
  list_phenos
}


#' Generate phenotypes emulating real phenotypic data at hand.
#'
#' This function simulates phenotypes from real phenotypic data based on their
#' sample correlation structure.
#'
#' @param n Number of observations.
#' @param real_phenos Matrix of real phenotypes (rows observations, columns
#'   phenotypic variables), without missing values.
#' @param bl_lgth Number of variables per block for reproducing the dependence
#'   structure of real phenotypes. Must be between 2 and d. Must be small enough
#'   (e.g. 1000) for tractability reasons. Default is \code{NULL} for a single
#'   block.
#' @param d Number of phenotypes. If \code{NULL}, the total number of phenotypes
#'   available in the real_phenos matrix is used (default).
#' @param n_cpus Number of CPUs used when simulating phenotypes by blocks. Set
#'   to 1 for serial execution.
#' @param user_seed Seed set for reproducibility. Default is \code{NULL}, no
#'   seed set.
#'
#' @return An object of class "\code{sim_phenos}".
#'  \item{phenos}{Matrix containing the generated phenotypic data.}
#'  \item{var_err}{Vector containing the sample phenotypic variances.}
#'  \item{ind_bl}{List of length given by the number of blocks, containing the
#'                indices of the phenotypes in each block. Is \code{NULL} if
#'                \code{cor_type} is \code{NULL}.}
#'
#' @examples
#' user_seed <- 123; set.seed(user_seed)
#' n <- 500; d <- 1000; var_err <- runif(d, min = 0.1, max = 0.4)
#' cor_type <- "equicorrelated"; vec_rho <- runif(100, min = 0.25, max = 0.95)
#' list_fake_real_phenos <- generate_phenos(n, d, var_err, cor_type, vec_rho,
#'                                          n_cpus = 2, user_seed = user_seed)
#' list_phenos <- replicate_real_phenos(n, list_fake_real_phenos$phenos,
#'                                      bl_lgth = 100, n_cpus = 2,
#'                                      user_seed = user_seed)
#'
#' @seealso \code{\link{generate_snps}}, \code{\link{replicate_real_snps}},
#'   \code{\link{generate_phenos}}, \code{\link{generate_dependence}}
#'
#' @export
#'
replicate_real_phenos <- function(n, real_phenos, bl_lgth = NULL, d = NULL,
                                  n_cpus = 1, user_seed = NULL) {


  check_structure_(user_seed, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(user_seed)){
    RNGkind("L'Ecuyer-CMRG") # ensure reproducibility when using mclapply
    set.seed(user_seed)
  }

  check_structure_(real_phenos, "matrix", "numeric")

  check_structure_(n_cpus, "vector", "numeric", 1)
  check_natural_(n_cpus)

  if (n_cpus > 1) {
    n_cpus_avail <- parallel::detectCores()
    if (n_cpus > n_cpus_avail) {
      n_cpus <- n_cpus_avail
      warning(paste("The number of CPUs specified exceeds the number of CPUs ",
                    "available on the machine. The latter has been used instead.", sep=""))
    }
  }

  check_structure_(n, "vector", "numeric", 1) # can be larger than the number of available observations in real_pheno!
  check_natural_(n)

  check_structure_(d, "vector", "numeric", 1, null_ok = TRUE)
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

    if (bl_lgth == 1) stop("Provided block length must be larger than 1")
    if (bl_lgth > d)
      stop(paste("Provided block length must be smaller or equal to the number ",
                 "of phenotypes available: ", bl_lgth, " > ", d, sep = ""))
  }

  n_bl <- floor(d / bl_lgth)
  ind_bl <- make_chunks_(1:d, n_bl)

  phenos <- parallel::mclapply(1:n_bl, function(bl) {
    R <- cor(real_phenos[, ind_bl[[bl]]])
    R <- Matrix::nearPD(R, corr = TRUE, do2eigen = TRUE)$mat
    d_bl <- ncol(R)
    L <- t(chol(R))
    tZ <- matrix(rnorm(n * d_bl), nrow = d_bl, ncol = n)
    as.matrix(t(L %*% tZ)) # Gaussian variables
  }, mc.cores = n_cpus)

  phenos <- do.call(cbind, phenos)

  rownames(phenos) <- paste("ind_", 1:n, sep = "")
  colnames(phenos) <- paste("pheno_", 1:d, sep = "")

  var_err <- apply(phenos, 2, var) # empirical error variance
  names(var_err) <- colnames(phenos)

  ind_bl <- NULL # block chuncks do not necessarily represent blocks of correlated phenotypes.

  list_phenos <- create_named_list_(phenos, var_err, ind_bl)
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

  pat <- matrix(FALSE, nrow=p, ncol=d)

  for(ind_j in ind_p0) {

    # random permutation, so that two different SNPs can be associated with a given
    # block with different probabilities
    vec_prob_sh_perm <- sample(vec_prob_sh, n_chunks_ph, replace = TRUE)
    for(ch in 1:n_chunks_ph) {
      ind_ch_d0 <- intersect(ind_d0, chunks_ph[[ch]])
      pat[ind_j, ind_ch_d0] <- sample(c(TRUE, FALSE),
                                      length(ind_ch_d0), replace = TRUE,
                                      prob=c(vec_prob_sh_perm[ch], 1 - vec_prob_sh_perm[ch]))
    }

    if (all(!pat[ind_j, ind_d0]))
      pat[ind_j, ind_d0][sample(1:length(ind_d0), 1)] <- TRUE # each active covariate must be
    # associated with at least one response
  }

  for(ind_k in ind_d0) {
    if (all(!pat[ind_p0, ind_k]))
      # each active covariate must be associated with at least one response
      pat[ind_p0, ind_k][sample(1:length(ind_p0), 1)] <- TRUE
  }
  pat
}


generate_eff_sizes_ <- function(d, p, ind_d0, ind_p0, vec_prob_sh, vec_maf,
                                pve_per_snp, max_tot_pve, var_err, chunks_ph) {

  # pve_per_snp average variance explained per snp
  check_structure_(ind_d0, "vector", "numeric", null_ok = TRUE)
  check_structure_(ind_p0, "vector", "numeric", null_ok = TRUE)

  d0 <- length(ind_d0)
  p0 <- length(ind_p0)

  if (xor(d0 < 1, p0 < 1))
    stop("ind_d0 and ind_p0 must either have both length 0 (no association) or both length >= 1.")

  if (d0 < 1 & p0 < 1) {

    pat <- matrix(FALSE, nrow = p, ncol = d)
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

    check_structure_(var_err, "vector", "numeric", d)
    check_positive_(var_err)

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

      tot_var_expl <- pve_per_snp * p0_k * var_err[k] / (1 - pve_per_snp * p0_k)

      vec_maf_act <- vec_maf[pat[,k]]
      vec_var_act <- 2 * vec_maf_act * (1 - vec_maf_act)

      beta_k <- rep(0.0, p)
      beta_k[pat[,k]] <- sqrt((tot_var_expl + var_err[k]) * vec_pve_per_snp / vec_var_act)

      # switches signs with probabilty 0.5
      beta_k[pat[,k]] <- sample(c(1, -1), p0_k, replace = TRUE) * beta_k[pat[,k]]

      beta_k
    })

  }

  create_named_list_(beta, pat, pve_per_snp)

}

#' Generate pleiotropic associations between SNPs and phenotypes.
#'
#' This function sets the association pattern and the effect sizes between SNP
#' and phenotype objects previously obtained from the functions
#' \code{\link{generate_snps}} or \code{\link{replicate_real_snps}}, and
#' \code{\link{generate_phenos}} or \code{\link{replicate_real_phenos}}. It
#' therefore adds a genetic contribution to the phenotypic data.
#'
#' The user can provide using the argument \code{vec_prob_sh} a selection of
#' probabilities describing the propensity with which a given active SNP (i.e.,
#' associated with at least one phenotype) will be associated with active
#' phenotypes (i.e., associated with at least one SNP) of given phenotypic
#' blocks. More precisely, for each active SNP and each phenotypic block, a
#' value from this vector is selected uniformly at random; for instance a large
#' probability implies that the SNPs is highly likely to be associated with
#' each active phenotype in the block. If a single value is provided, all active
#' SNPs will have the same probability to be associated with active phenotypes
#' of all blocks.
#'
#' The user can provide either argument \code{pve_per_snp}, specifying the
#' average proportion of phenotypic variance explained per active SNP for a
#' given active phenotype, or \code{max_tot_pve}, specifying the maximum value
#' for an active phenotype of its proportion of variance explained by the
#' cummulated genetic effects. If both \code{pve_per_snp} and \code{max_tot_pve}
#' are \code{NULL}, the proportion of phenotypic variance explained per SNP is
#' set to its maximum value so that the total proportion of variance explained
#' for the phenotypes are all below 1. Individual proportions of variance
#' explained are drawn from a Beta distribution with shape parameters 2 and 5,
#' putting more weights on smaller effects.
#'
#' @param list_snps An object of class "sim_snps" containing simulated SNP data
#'   and their corresponding sample minor allele frequencies. It must be
#'   obtained from the function \code{\link{generate_snps}} or
#'   \code{\link{replicate_real_snps}}.
#' @param list_phenos An object of class "sim_pheno" containing simulated
#'   phenotypic data, their sample variance and block structure information.
#'   It must be obtained from the function \code{\link{generate_phenos}} or
#'   \code{\link{replicate_real_phenos}}.
#' @param ind_d0 A vector of indices specifying the position of the "active"
#'   phenotypes (i.e., which will be associated with at least one SNP). Must
#'   range between 1 and \code{ncol(list_phenos$phenos)}.
#' @param ind_p0 A vector of indices specifying the position of the "active"
#'   SNPs (i.e., which will be associated with at least one phenotype). Must
#'   range between 1 and \code{ncol(list_snps$snps)}.
#' @param vec_prob_sh Vector providing a set of probabilities with which an
#'   active SNP is associated with an additional active phenotype in a given
#'   phenotypic block. See Details section.
#' @param pve_per_snp Average proportion of phenotypic variance explained by
#'   each active SNP (for an active phenotype). Must be \code{NULL} if
#'   \code{max_tot_pve} is provided. See Details section.
#' @param max_tot_pve Maximum proportion of phenotypic variance explained by the
#'   active SNPs across all phenotypes. Must be \code{NULL} if
#'   \code{pve_per_snp} is provided. See Details section.
#' @param user_seed Seed set for reproducibility. Default is \code{NULL}, no
#'   seed set.
#'
#' @return An object of class "\code{sim_data}".
#'  \item{phenos}{Matrix containing the updated phenotypic data (whose variance
#'                is now partly explained by genetic effects).}
#'  \item{snps}{Matrix containing the original SNPs data.}
#'  \item{beta}{Matrix containing the generated effect sizes between the SNPs
#'             (rows) and phenotypes (columns).}
#'  \item{pat}{Matrix of booleans specifying the generated association pattern
#'             between the SNPs (rows) and phenotypes (columns).}
#'  \item{pve_per_snp}{Average proportion of phenotypic variance explained by
#'                     each active SNP (for an active phenotype).}
#'
#' @seealso \code{\link{generate_snps}}, \code{\link{replicate_real_snps}},
#'   \code{\link{generate_phenos}}, \code{\link{replicate_real_phenos}}
#'
#'
#' @examples
#' user_seed <- 123; set.seed(user_seed)
#' n <- 500; p <- 5000; p0 <- 200; d <- 500; d0 <- 400
#' list_snps <- generate_snps(n = n, p = p)
#' list_phenos <- generate_phenos(n = n, d = d, var_err = 0.25)
#'
#' dat <- generate_dependence(list_snps, list_phenos, ind_d0 = sample(1:d, d0),
#'                            ind_p0 = sample(1:p, p0), vec_prob_sh = 0.05,
#'                            max_tot_pve = 0.5)
#'
#' @export
#'
generate_dependence <- function(list_snps, list_phenos, ind_d0, ind_p0, vec_prob_sh,
                                pve_per_snp = NULL, max_tot_pve = NULL,
                                user_seed = NULL) {

  check_structure_(user_seed, "vector", "numeric", 1, null_ok = TRUE)
  if (!is.null(user_seed)){
    RNGkind("L'Ecuyer-CMRG") # ensure reproducibility when using mclapply
    set.seed(user_seed)
  }

  if (!is.null(pve_per_snp) & !is.null(max_tot_pve))
    stop("Either pve_per_snp or max_tot_pve must be NULL.")

  if (is.null(pve_per_snp) & is.null(max_tot_pve))
    warning(paste("As both pve_per_snp or max_tot_pve were provided as NULL, the ",
                  "pve per SNP was set to its maximum value so that the total ",
                  "pve for the responses are all below 1.", sep = ""))

  if (!inherits(list_snps, "sim_snps"))
    stop(paste("The provided list_snps must be an object of class ``sim_snps''. \n",
               "*** You must either use the function generate_snps to simulate snps ",
               "under Hardy-Weinberg equilibrium or the function replicate_real_snps ",
               "to simulate SNPs from real SNP data, by replicating their minor ",
               "allele frequencies and linkage desequilibrium structure. ***",
               sep=""))


  if (!inherits(list_phenos, "sim_phenos"))
    stop(paste("The provided list_phenos must be an object of class ``sim_phenos''. \n",
               "*** You must either use the function generate_phenos to simulate ",
               "phenotypes from (possibly correlated) gaussian variables or the ",
               "function replicate_real_phenos to simulate phenotypes from real ",
               "phenotypic data, by replicating their correlation structure. ***",
               sep=""))

  with(c(list_snps, list_phenos), {

    n <- nrow(snps)
    p <- ncol(snps)

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
                                    pve_per_snp, max_tot_pve, var_err,
                                    chunks_ph = ind_bl)
    with(list_eff, {
      phenos <- phenos + snps %*% beta

      list_data <- create_named_list_(phenos, snps, beta, pat, pve_per_snp)
      class(list_data) <- "sim_data"
      list_data
    })
  })
}
