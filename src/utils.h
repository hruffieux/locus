#ifndef LOCUS_UTILS_H_
#define LOCUS_UTILS_H_

#include <RcppEigen.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;

// NOTE: can't be used as not properly copied to RcppExports....
//
// typedef Eigen::ArrayXd Arr1D;
// typedef Eigen::ArrayXXd Arr2D;
// typedef Eigen::Map<Eigen::ArrayXd> MapArr1D;
// typedef Eigen::Map<Eigen::ArrayXXd> MapArr2D;
// typedef Eigen::Map<Eigen::MatrixXd> MapMat;


Eigen::ArrayXd logOnePlusExp(const  Eigen::ArrayXd& x);

#endif // LOCUS_UTILS_H_
