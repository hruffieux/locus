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

  if (type != "string") {
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

log_sum_exp_ <- function(x) { # avoid numerical underflow or overflow
  if ( max(abs(x)) > max(x) )
    offset <- min(x)
  else
    offset <- max(x)
  log(sum(exp(x - offset))) + offset
}

log_sum_exp_mat_ <- function(list_vec) { # avoid numerical underflow or overflow
                                         # for a list of two matrices (or vectors)
  stopifnot(length(list_vec) == 2)

  a <- list_vec[[1]]
  b <- list_vec[[2]]
  rm(list_vec)

  bool_max <- (a - b > 0)
  M <- b
  M[bool_max] <- a[bool_max]
  m <- a
  m[bool_max] <- b[bool_max]

  absa <- abs(a)
  absb <- abs(b)
  bool_absmax <- (absa - absb > 0)
  absM <- absb
  absM[bool_absmax] <- absb[bool_absmax]

  bool_off <- absM > M
  offset <- M
  offset[bool_off] <- m[bool_off]

  log(exp(a - offset) + exp(b - offset)) + offset

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
