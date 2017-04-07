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
    err_mess <- paste(deparse(substitute(x)), " must be positive.", sep="")
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


log_one_plus_exp_ <- function(x) { # computes log(1 + exp(x)) avoiding
                                     # numerical overflow
  m <- x
  m[x < 0] <- 0

  log(exp(x - m) + exp(- m)) + m
}


inv_mills_ratio_ <- function(Y, U) {

  m <- matrix(NA, nrow = nrow(U), ncol = ncol(U))

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

# entropy_ <- function(Y, U) {
#
#   log((2 * pi * exp(1))^(1/2) *
#            exp(Y * pnorm(U, log.p = TRUE) +
#                  (1-Y) * pnorm(U, lower.tail = FALSE, log.p = TRUE))) -
#   U * inv_mills_ratio_(Y, U) / 2
#
# }

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
