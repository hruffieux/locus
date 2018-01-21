# This file is part of the `locus` R package:
#     https://github.com/hruffieux/locus
#

# Diverse utility functions implementing sanity checks, basic preprocessing,
# and ticks to prevent overflow/underflow.
#

check_natural_ <- function(x, eps = .Machine$double.eps^0.75){
  if (any(x < eps | abs(x - round(x)) > eps)) {
    stop(paste(deparse(substitute(x)),
               " must be natural.", sep=""))
  }
}

check_positive_ <- function(x, eps = .Machine$double.eps^0.75){
  if (any(x < eps)) {
    err_mess <- paste(deparse(substitute(x)), " must be positive, greater than ",
                      format(eps, digits = 3), ".", sep="")
    if (length(x) > 1) err_mess <- paste("All entries of ", err_mess, sep="")
    stop(err_mess)
  }
}

check_zero_one_ <- function(x){
  if (any(x < 0) | any(x > 1)) {
    err_mess <- paste(deparse(substitute(x)), " must lie between 0 and 1.", sep="")
    if (length(x) > 1) err_mess <- paste("All entries of ", err_mess, sep="")
    stop(err_mess)
  }
}

check_structure_ <- function(x, struct, type, size = NULL,
                             null_ok = FALSE,  inf_ok = FALSE, na_ok = FALSE) {
  if (type == "double") {
    bool_type <-  is.double(x)
    type_mess <- "a double-precision "
  } else if (type == "integer") {
    bool_type <- is.integer(x)
    type_mess <- "an integer "
  } else if (type == "numeric") {
    bool_type <- is.numeric(x)
    type_mess <- "a numeric "
  } else if (type == "logical") {
    bool_type <- is.logical(x)
    type_mess <- "a boolean "
  } else if (type == "string") {
    bool_type <- is.character(x)
    type_mess <- "string "
  }

  bool_size <- TRUE # for case size = NULL (no assertion on the size/dimension)
  size_mess <- ""
  if (struct == "vector") {
    bool_struct <- is.vector(x) & (length(x) > 0) # not an empty vector
    if (!is.null(size)) {
      bool_size <- length(x) %in% size
      size_mess <- paste(" of length ", paste(size, collapse=" or "), sep = "")
    }
  } else if (struct == "matrix") {
    bool_struct <- is.matrix(x) & (length(x) > 0) # not an empty matrix
    if (!is.null(size)) {
      bool_size <- all(dim(x) == size)
      size_mess <- paste(" of dimension ", size[1], " x ", size[2], sep = "")
    }
  }

  correct_obj <- bool_struct & bool_type & bool_size

  bool_null <- is.null(x)

  if (!is.list(x) & type != "string") {
    na_mess <- ""
    if (!na_ok) {
      if (!bool_null) correct_obj <- correct_obj & !any(is.na(x))
      na_mess <- " without missing value"
    }

    inf_mess <- ""
    if (!inf_ok) {
      if (!bool_null) correct_obj <- correct_obj & all(is.finite(x[!is.na(x)]))
      inf_mess <- ", finite"
    }
  } else {
    na_mess <- ""
    inf_mess <- ""
  }

  null_mess <- ""
  if (null_ok) {
    correct_obj <- correct_obj | bool_null
    null_mess <- " or must be NULL"
  }

  if(!(correct_obj)) {
    stop(paste(deparse(substitute(x)), " must be a non-empty ", type_mess, struct,
               size_mess, inf_mess, na_mess, null_mess, ".", sep = ""))
  }
}


create_named_list_ <- function(...) {
  setNames(list(...), as.character(match.call()[-1]))
}



get_annealing_ladder_ <- function(anneal, verbose) {

  # ladder set following:
  # Importance Tempering, Robert B. Gramacy & Richard J. Samworth, pp.9-10

  k_m <- 1 / anneal[2]
  m <- anneal[3]

  if(anneal[1] == 1) {

    type <- "geometric"

    delta_k <- k_m^(1 / (1 - m)) - 1

    ladder <- (1 + delta_k)^(1 - m:1)

  } else if (anneal[1] == 2) { # harmonic spacing

    type <- "harmonic"

    delta_k <- ( 1 / k_m - 1) / (m - 1)

    ladder <- 1 / (1 + delta_k * (m:1 - 1))

  } else { # linear spacing

    type <- "linear"

    delta_k <- (1 - k_m) / (m - 1)

    ladder <- k_m + delta_k * (1:m - 1)
  }

  if (verbose)
    cat(paste0("** Annealing with ", type," spacing ** \n\n"))

  ladder

}


log_one_plus_exp_ <- function(x) { # computes log(1 + exp(x)) avoiding
                                   # numerical overflow
  m <- x
  m[x < 0] <- 0

  log(exp(x - m) + exp(- m)) + m
}


log_sigmoid_ <- function(chi) {

  - log(1 + exp(- chi)) # chi is always positive so no overflow possible (underflow neither, thanks to the "+1")

}

log_det <- function(list_mat) {

  if (is.list(list_mat)) {
    sapply(list_mat, function(mat) {
      log_det <- determinant(mat, logarithm = TRUE)
      log_det$modulus * log_det$sign
    })
  } else {
    log_det <- determinant(list_mat, logarithm = TRUE)
    log_det$modulus * log_det$sign
  }

}

inv_mills_ratio_ <- function(Y, U) {

  if (is.matrix(U)) m <- matrix(NA, nrow = nrow(U), ncol = ncol(U))
  else m <- rep(NA, length(U))

  U_1 <- U[Y==1]
  m_1 <- exp(dnorm(U_1, log = TRUE) - pnorm(U_1, log.p = TRUE))
  m_1[m_1 < -U_1] <- -U_1

  m[Y==1] <- m_1


  U_0 <- U[Y==0]
  m_0 <- - exp(dnorm(U[Y==0], log = TRUE) - pnorm(U[Y==0], lower.tail = FALSE, log.p = TRUE))
  m_0[m_0 > -U_0] <- -U_0

  m[Y==0] <- m_0

  m

}


log_sum_exp_ <- function(x) {
  # Computes log(sum(exp(x))

  if ( max(abs(x)) > max(x) )
    offset <- min(x)
  else
    offset <- max(x)
  log(sum(exp(x - offset))) + offset

}


# entropy_ <- function(Y, U) {
#
#   log((2 * pi * exp(1))^(1/2) *
#            exp(Y * pnorm(U, log.p = TRUE) +
#                  (1-Y) * pnorm(U, lower.tail = FALSE, log.p = TRUE))) -
#   U * inv_mills_ratio_(Y, U) / 2
#
# }


# Functions for hyperparameter settings in dual_core (similarly to what is done in HESS)
#
E_Phi_X <- function(mu, s2, lower_tail = TRUE) {

  pnorm(mu / sqrt(1 + s2), lower.tail = lower_tail)

}

E_Phi_X_2 <- function(mu, s2) {

  pnorm(mu / sqrt(1 + s2)) -
    2 * PowerTOST::OwensT(mu / sqrt(1 + s2), 1 / sqrt(1 + 2 * s2))

}

get_V_p_t <- function(mu, s2, p) {
  p * (p - 1) * E_Phi_X_2(mu, s2) -
    p^2 * E_Phi_X(mu, s2)^2 +
    p * E_Phi_X(mu, s2)
}


get_mu <- function(E_p_t, s2, p) {

  sqrt(1 + s2) * qnorm(1- E_p_t / p)

}


get_n0_t02 <- function(d, p, p_star) {
  
  E_p_t <- p_star[1]
  V_p_t <- min(p_star[2], floor(2 * p / 3))

  dn <- 1e-6
  up <- 1e5
  
  # Get n0 and t02 similarly as for a_omega_t and b_omega_t in HESS
  # (specify expectation and variance of number of active predictors per response)
  #
  # Look at : gam_st | theta_s = 0
  #
  tryCatch(t02 <- uniroot(function(x)
    get_V_p_t(get_mu(E_p_t, x, p), x, p) - V_p_t,
    interval = c(dn, up))$root,
    error = function(e) {
      stop(paste0("No hyperparameter values matching the expectation and variance ",
                  "of the number of active predictors per responses supplied in p0_av.",
                  "Please change p0_av."))
    })
  
  # n0 sets the level of sparsity.
  n0 <- get_mu(E_p_t, t02, p)
  n0 <- rep(-n0, d)
  
  create_named_list_(n0, t02)
}



rm_constant_ <- function(mat, verbose) {

  bool_cst <- is.nan(colSums(mat))

  if (any(bool_cst)) {

    rmvd_cst <- colnames(mat)[bool_cst]

    if (verbose) {
      if (sum(bool_cst) < 50) {
        cat(paste("Variable(s) ", paste(rmvd_cst, collapse=", "),
                  " constant across subjects. \n",
                  "Removing corresponding column(s) and saving its/their id(s) ",
                  "in the function output ... \n\n",
                  sep=""))
      } else {
        cat(paste(sum(bool_cst), " variables constant across subjects. \n",
                  "Removing corresponding column(s) and saving their ids ",
                  "in the function output ... \n\n",
                  sep=""))
      }
    }

    mat <- mat[, !bool_cst, drop = FALSE]
  } else {
    rmvd_cst <- NULL
  }

  create_named_list_(mat, bool_cst, rmvd_cst)
}

rm_collinear_ <- function(mat, verbose) {

  bool_coll <- duplicated(mat, MARGIN = 2)

  if (any(bool_coll)) {

    mat_coll <- mat[, bool_coll, drop = FALSE]
    rmvd_coll <- colnames(mat_coll)

    if (verbose) {
      if (length(rmvd_coll) < 50) {
        cat(paste("Presence of collinear variable(s). ",
                  paste(rmvd_coll, collapse=", "), " redundant. \n",
                  "Removing corresponding column(s) and saving its/their id(s) ",
                  "in the function output ... \n",
                  sep=""))
      } else {
        cat(paste("Presence of collinear variables. ", length(rmvd_coll),
                  " redundant.\n", "Removing corresponding columns and saving ",
                  "their ids in the function output ... \n",
                  sep=""))
      }
    }

    # associate to each removed replicate the name of the covariate with which
    # it is duplicated and that is kept in the dataset
    bool_with_coll <- duplicated(mat, MARGIN = 2, fromLast = TRUE) & !bool_coll
    mat_with_coll <- mat[, bool_with_coll, drop = FALSE]

    assoc_coll <- colnames(mat_with_coll)[match(data.frame(mat_coll),
                                                data.frame(mat_with_coll))]
    names(rmvd_coll) <- assoc_coll

    mat <- mat[, !bool_coll, drop = FALSE]

  } else {
    rmvd_coll <- NULL
  }

  create_named_list_(mat, bool_coll, rmvd_coll)
}

make_chunks_ <- function(x, n_g) split(x, factor(sort(rank(x) %% n_g)))


Q_approx <- function(x, eps1 = 1e-30, eps2 = 1e-7) {
  
  if(x <= 1) {
    
    gsl::expint_E1(x) * exp(x)
    
  } else {
    
    f_p <- eps1
    C_p <- eps1
    D_p <- 0
    Delta <- 2 + eps2
    j <- 1
    
    while( abs(Delta-1) >= eps2 ) {
      
      j <- j+1
      
      D_c <- x + 2*j - 1 - ((j-1)^{2}) * D_p
      C_c <- x + 2*j - 1 - ((j-1)^{2}) / C_p
      D_c <- 1 / D_c
      
      Delta <- C_c * D_c
      f_c <- f_p * Delta
      f_p <- f_c
      C_p <- C_c
      D_p <- D_c
    }
    
    1/(x + 1 + f_c)
  }
}


compute_integral_hs_ <- function(alpha, beta, m, n, Q_ab) {
  
  # computes int_0^infty x^n (1 + alpha * x)^(-m) * exp(- beta * x) dx
  # for m = n or m = n + 1, n natural n > 0, beta > 0 (if n = 0, then = Q_ab)
  
  # Q_ab = Q_approx(alpha / beta) # precomputed to avoid computing it too many times
  # = exp(beta/alpha) * E_1(beta/alpha)
  
  
  if (m == n) {
    
    out <- alpha^(-n) * beta^(-1) 
    
    if (n == 1) {
      
      out <- out - alpha^(-2) * Q_ab
      
    } else if (n == 2) {
      
      out <- out - 2 * alpha^(-3) * Q_ab + alpha^(-3) - alpha^(-4) * beta * Q_ab
      
    } else if (n == 3) {
      
      out <- out - 3 * alpha^(-4) * Q_ab + 
        3 * alpha^(-5) * (alpha - beta * Q_ab) -
        alpha^(-6) / 2 * (-beta * alpha + alpha^2 + beta^2 * Q_ab)
    
    } else if (n == 4){

      v1 <- c(-n * log(alpha) - log(beta),
              log(4) - 6 * log(alpha) + log(alpha),
              log(4) - 7 * log(alpha) - log(2) + 2 * log(alpha),
              log(4) - 7 * log(alpha) - log(2) + 2 * log(beta) + log(Q_ab),
              -8 * log(alpha) - log(6) + 2 * log(beta) + log(alpha),
              -8 * log(alpha) - log(6) + log(2) + 3 * log(alpha))
      
      v2 <- c(log(4) - 5 * log(alpha) + log(Q_ab),
              log(4) - 6 * log(alpha) + log(beta) + log(Q_ab),
              log(4) - 7 * log(alpha) - log(2)  + log(beta),
              -8 * log(alpha) - log(6) + log(beta) + 2 * log(alpha))
      
      out <- exp(log_sum_exp_(v1)) - exp(log_sum_exp_(v2))
    
    } else {

      v1 <- c(-n * log(alpha) - log(beta), 
              -n * log(alpha) + log(n) + unlist(sapply(2:(n-1), function(k) {
                -k * log(alpha) - lfactorial(k-1) +
                  sapply(seq(1, k-1, by = 2), function(j) {
                    lfactorial(j-1) + (k-j-1) * log(beta) + j * log(alpha)}) })),
              -2*n *log(alpha) - lfactorial(n-1) +
                sapply(seq(1, n-1, by = 2), function(j) {
                  lfactorial(j-1) + (n-j-1) * log(beta) + j * log(alpha)})
      )
      
      v2 <- c(log(n) - (n + 1) * log(alpha) + log(Q_ab), 
              -n * log(alpha) + log(n) + unlist(sapply(3:(n-1), function(k) { # k = 2 doesn't contribute
                -k * log(alpha) - lfactorial(k-1) +
                  sapply(seq(2, k-1, by = 2), function(j) {
                    lfactorial(j-1) * (k-j-1) * log(beta) + j * log(alpha)}) })),
              -n * log(alpha) + log(n) + sapply(2:(n-1), function(k) {
                -k * log(alpha) - lfactorial(k-1) + (k-1) * log(beta) + log(Q_ab)}),
              - 2*n * log(alpha) - lfactorial(n-1) +
                sapply(seq(2, n-1, by = 2), function(j) {
                  lfactorial(j-1) + (n-j-1) * log(beta) + j * log(alpha)}), 
              -2*n * log(alpha) - lfactorial(n-1) + (n-1) * log(beta) + log(Q_ab)
      )
      
      out <- exp(log_sum_exp_(v1)) - exp(log_sum_exp_(v2))

    }
  
    
  } else if (m == n + 1) { # not stable for n >= 4, i.e., can't be used for df = 9 and higher.
    
    if (n == 1) {

      out <- alpha^(-2) * Q_ab - alpha^(-2) +  alpha^(-3) * beta * Q_ab

    } else if (n == 2){

      v1 <- c(-3 * log(alpha) + log(Q_ab),
              -5 * log(alpha) - log(2) + 2 * log(alpha),
              -5 * log(alpha) - log(2) + 2 * log(beta) + log(Q_ab),
              -2 * log(alpha) + log(2) - 2 * log(alpha) + log(beta) + log(Q_ab)
      )

      v2 <- c(-5 * log(alpha) - log(2) + log(beta) + log(alpha),
              -2 * log(alpha) + log(2) - 2 * log(alpha) + log(alpha)
      )


      out <- exp(log_sum_exp_(v1)) - exp(log_sum_exp_(v2))


     } else {
      
       
       v1 <- c(-(n+1) * log(alpha) + log(Q_ab),
               -(2*n+1) * log(alpha) - lfactorial(n) +
                 sapply(seq(2, n, by = 2), function(j) {
                   lfactorial(j-1) + (n-j) * log(beta) + j * log(alpha)}),
               -(2*n+1) * log(alpha) - lfactorial(n) + n * log(beta) + log(Q_ab),
               -n * log(alpha) + log(n) + unlist(sapply(2:(n-1), function(k) {
                 -(1+k) * log(alpha) - lfactorial(k) +
                   sapply(seq(2, k, by = 2), function(j) {
                     lfactorial(j-1) + (k-j) * log(beta) + j * log(alpha) })
               })),
               -n * log(alpha) + log(n) + sapply(1:(n-1), function(k) {
                 -(1+k) * log(alpha) - lfactorial(k) + k * log(beta) + log(Q_ab)
               })
       )

       v2 <- c(-(2*n+1) * log(alpha) - lfactorial(n) +
                 sapply(seq(1, n, by = 2), function(j) {
                   lfactorial(j-1) + (n-j) * log(beta) + j * log(alpha)}),
               -n * log(alpha) + log(n) + unlist(sapply(1:(n-1), function(k) {
                 -(1+k) * log(alpha) - lfactorial(k) +
                   sapply(seq(1, k, by = 2), function(j) {
                     lfactorial(j-1) + (k-j) * log(beta) + j * log(alpha) })
               }))

       )
       
       out <- exp(log_sum_exp_(v1)) - exp(log_sum_exp_(v2))
        
  }

  } else {
    
    stop("Invalid value of m, must be n or n + 1.")
    
  }
 
  out
}


checkpoint_ <- function(it, checkpoint_path, 
                        gam_vb, converged, lb_new, lb_old, b_vb = NULL,
                        mu_rho_vb = NULL, mu_theta_vb = NULL, om_vb = NULL,
                        S0_inv_vb = NULL, rate = 100) {
  
  if (!is.null(checkpoint_path) && it %% rate == 0) {
    
    diff_lb <- abs(lb_new - lb_old)
    
    tmp_vb <- create_named_list_(gam_vb, converged, it, lb_new, diff_lb, 
                                 b_vb, mu_theta_vb, mu_rho_vb, om_vb, S0_inv_vb)
    
    file_save <- paste0(checkpoint_path, "tmp_output_it_", it, ".RData")
    
    save(tmp_vb, file = file_save)
    
    old_file_clean_up <- paste0(checkpoint_path, "tmp_output_it_", it - 2 * rate, ".RData") # keep only the last two for comparison
    
    if (file.exists(old_file_clean_up)) 
      file.remove(old_file_clean_up)

  }
    
}


checkpoint_clean_up_ <- function(checkpoint_path) {
  
  if (!is.null(checkpoint_path)) {
    
    old_files_clean_up <- list.files(path = checkpoint_path, pattern = "tmp_output_it_")
    
    sapply(old_files_clean_up, function(ff) {
      if (file.exists(file.path(checkpoint_path, ff))) 
        file.remove(file.path(checkpoint_path, ff))
    })
    
  } 
  
}
 
