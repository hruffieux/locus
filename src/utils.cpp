#include "utils.h"

Eigen::ArrayXd logOnePlusExp(const  Eigen::ArrayXd& x){
  Eigen::ArrayXd m = x;
  for (int i = 0; i < x.size(); ++i) {
    if (x[i] < 0) m[i] = 0;
  }

  return log(exp(x - m) + exp(-m)) + m;

}
