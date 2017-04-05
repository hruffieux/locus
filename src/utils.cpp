#include "utils.h"

Arr1D logOnePlusExp(const Arr1D& x){

  Arr1D m = x;
  for (int i = 0; i < x.size(); ++i) {
    if (x[i] < 0) m[i] = 0;
  }

  return log(exp(x - m) + exp(-m)) + m;

}


Arr2D logOnePlusExpMat(const Arr2D& x){

  Arr2D m = x;
  for (int i = 0; i < x.rows(); ++i) {

    for (int j = 0; j < x.cols(); ++j) {

     if (x(i, j) < 0) m(i, j) = 0;

    }

  }

  return log(exp(x - m) + exp(-m)) + m;

}
