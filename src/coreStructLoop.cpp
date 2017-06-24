/*
 *
 * This file is part of the `locus` R package:
 *     https://github.com/hruffieux/locus
 *
 * Functions for computationally expensive updates for structured regression.
 *
 * These functions use Eigen::Map to pass large matrices by reference from R.
 * Given dimensionalities involved in some applications, copying such matrices
 * would imply a prohibitive RAM overconsumption.
 *
 */

#include "utils.h"

// for locus_struct_core function
// [[Rcpp::export]]
void coreStructLoop(const MapMat X,
                    const MapMat Y,
                    MapArr2D gam_vb,
                    const MapArr1D log_Phi_mu_theta_vb,
                    const MapArr1D log_1_min_Phi_mu_theta_vb,
                    const double log_sig2_inv_vb,
                    const MapArr1D log_tau_vb,
                    MapMat m1_beta,
                    MapMat mat_x_m1,
                    MapArr2D mu_beta_vb,
                    const MapArr1D sig2_beta_vb,
                    const MapArr1D tau_vb,
                    const MapArr1D shuffled_ind) {

  const Arr1D c = -(log_tau_vb + log_sig2_inv_vb + log(sig2_beta_vb) )/ 2;

  for (int i = 0; i < X.cols(); ++i) {

    int j = shuffled_ind[i];

    mat_x_m1.noalias() -= X.col(j) * m1_beta.row(j);

    mu_beta_vb.row(j) = sig2_beta_vb * tau_vb *
      ((Y - mat_x_m1).transpose() * X.col(j)).array();

    gam_vb.row(j) = exp(-logOnePlusExp(log_1_min_Phi_mu_theta_vb(j) - log_Phi_mu_theta_vb(j) -
      mu_beta_vb.row(j).square() / (2 * sig2_beta_vb.transpose()) + c.transpose()));

    m1_beta.row(j) = mu_beta_vb.row(j) * gam_vb.row(j);

    mat_x_m1.noalias() += X.col(j) * m1_beta.row(j);

  }

}

