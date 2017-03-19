#include "utils.h"

Arr1D logOnePlusExp(const Arr1D& x){

  Arr1D m = x;
  for (int i = 0; i < x.size(); ++i) {
    if (x[i] < 0) m[i] = 0;
  }

  return log(exp(x - m) + exp(-m)) + m;

}
