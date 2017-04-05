#ifndef LOCUS_UTILS_H_
#define LOCUS_UTILS_H_

#include <RcppEigen.h>
#include "locus_types.h"

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

Arr1D logOnePlusExp(const  Arr1D& x);

Arr2D logOnePlusExpMat(const  Arr2D& x);

#endif // LOCUS_UTILS_H_
