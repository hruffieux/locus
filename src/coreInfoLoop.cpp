#include "utils.h"

// NOTE: typedefs can't be used as not properly copied to RcppExports...

// [[Rcpp::export]]
void coreInfoLoop(const Eigen::Map<Eigen::MatrixXd> X,
                  const Eigen::Map<Eigen::MatrixXd> Y,
                  Eigen::Map<Eigen::ArrayXXd> gam_vb,
                  const Eigen::Map<Eigen::ArrayXXd> log_Phi_mat_v_mu,
                  const Eigen::Map<Eigen::ArrayXXd> log_1_min_Phi_mat_v_mu,
                  const double log_sig2_inv_vb,
                  const Eigen::Map<Eigen::ArrayXd> log_tau_vb,
                  Eigen::Map<Eigen::MatrixXd> m1_beta,
                  Eigen::Map<Eigen::MatrixXd> mat_x_m1,
                  Eigen::Map<Eigen::ArrayXXd> mu_beta_vb,
                  const Eigen::Map<Eigen::ArrayXd> sig2_beta_vb,
                  const Eigen::Map<Eigen::ArrayXd> tau_vb) {

  const Eigen::ArrayXd c = -(log_tau_vb + log_sig2_inv_vb + log(sig2_beta_vb) )/ 2;

  for (int j = 0; j < X.cols(); ++j) {

    mat_x_m1.noalias() -= X.col(j) * m1_beta.row(j);

    mu_beta_vb.row(j) = sig2_beta_vb * tau_vb *
      ((Y - mat_x_m1).transpose() * X.col(j)).array();

    gam_vb.row(j) = exp(-logOnePlusExp(log_1_min_Phi_mat_v_mu.row(j) - log_Phi_mat_v_mu.row(j) -
      mu_beta_vb.row(j).square() / (2 * sig2_beta_vb.transpose()) + c.transpose()));

    m1_beta.row(j) = mu_beta_vb.row(j) * gam_vb.row(j);

    mat_x_m1.noalias() += X.col(j) * m1_beta.row(j);

  }

}


// for locus_z_info_core and locus_mix_info_core function
// [[Rcpp::export]]
void coreZInfoLoop(const Eigen::Map<Eigen::MatrixXd> X,
                   const Eigen::Map<Eigen::MatrixXd> Y,
                   Eigen::Map<Eigen::ArrayXXd> gam_vb,
                   const Eigen::Map<Eigen::ArrayXXd> log_Phi_mat_v_mu,
                   const Eigen::Map<Eigen::ArrayXXd> log_1_min_Phi_mat_v_mu,
                   const double log_sig2_inv_vb,
                   const Eigen::Map<Eigen::ArrayXd> log_tau_vb,
                   Eigen::Map<Eigen::MatrixXd> m1_beta,
                   Eigen::Map<Eigen::MatrixXd> mat_x_m1,
                   const Eigen::Map<Eigen::MatrixXd> mat_z_mu,
                   Eigen::Map<Eigen::ArrayXXd> mu_beta_vb,
                   const Eigen::Map<Eigen::ArrayXd> sig2_beta_vb,
                   const Eigen::Map<Eigen::ArrayXd> tau_vb) {

  const Eigen::ArrayXd c = -(log_tau_vb + log_sig2_inv_vb + log(sig2_beta_vb) )/ 2;

  for (int j = 0; j < X.cols(); ++j) {

    mat_x_m1.noalias() -= X.col(j) * m1_beta.row(j);

    mu_beta_vb.row(j) = sig2_beta_vb * tau_vb *
      ((Y - mat_x_m1 - mat_z_mu).transpose() * X.col(j)).array();

    gam_vb.row(j) = exp(-logOnePlusExp(log_1_min_Phi_mat_v_mu.row(j) - log_Phi_mat_v_mu.row(j) -
      mu_beta_vb.row(j).square() / (2 * sig2_beta_vb.transpose()) + c.transpose()));

    m1_beta.row(j) = mu_beta_vb.row(j) * gam_vb.row(j);

    mat_x_m1.noalias() += X.col(j) * m1_beta.row(j);

  }

}


// for locus_logit_core function
// [[Rcpp::export]]
void coreLogitInfoLoop(const Eigen::Map<Eigen::MatrixXd> X,
                       const Eigen::Map<Eigen::ArrayXXd> Y,
                       Eigen::Map<Eigen::ArrayXXd> gam_vb,
                       const Eigen::Map<Eigen::ArrayXXd> log_Phi_mat_v_mu,
                       const Eigen::Map<Eigen::ArrayXXd> log_1_min_Phi_mat_v_mu,
                       const double log_sig2_inv_vb,
                       Eigen::Map<Eigen::MatrixXd> m1_beta,
                       Eigen::Map<Eigen::ArrayXXd> mat_x_m1,
                       Eigen::Map<Eigen::ArrayXXd> mat_z_mu,
                       Eigen::Map<Eigen::ArrayXXd> mu_beta_vb,
                       const Eigen::Map<Eigen::ArrayXXd> psi_vb,
                       const Eigen::Map<Eigen::ArrayXXd> sig2_beta_vb) {

  for (int j = 0; j < X.cols(); ++j) {

    mat_x_m1.matrix().noalias() -= X.col(j) * m1_beta.row(j);

    mu_beta_vb.row(j) = sig2_beta_vb.row(j) * (X.col(j).transpose() * (Y - 2 * psi_vb * (mat_x_m1 + mat_z_mu)).matrix()).array();

    gam_vb.row(j) = exp(-logOnePlusExp(log_1_min_Phi_mat_v_mu.row(j) - log_Phi_mat_v_mu.row(j) -
      log_sig2_inv_vb / 2 - mu_beta_vb.row(j).square() / (2 * sig2_beta_vb.row(j)) -
      log(sig2_beta_vb.row(j)) / 2));

    m1_beta.row(j) = mu_beta_vb.row(j) * gam_vb.row(j);

    mat_x_m1.matrix().noalias() += X.col(j) * m1_beta.row(j);

  }

}



// for locus_probit_info_core function
// [[Rcpp::export]]
void coreProbitInfoLoop(const Eigen::Map<Eigen::MatrixXd> X,
                        const Eigen::Map<Eigen::MatrixXd> Wy,
                        Eigen::Map<Eigen::ArrayXXd> gam_vb,
                        const Eigen::Map<Eigen::ArrayXXd> log_Phi_mat_v_mu,
                        const Eigen::Map<Eigen::ArrayXXd> log_1_min_Phi_mat_v_mu,
                        const double log_sig2_inv_vb,
                        Eigen::Map<Eigen::MatrixXd> m1_beta,
                        Eigen::Map<Eigen::MatrixXd> mat_x_m1,
                        Eigen::Map<Eigen::MatrixXd> mat_z_mu,
                        Eigen::Map<Eigen::ArrayXXd> mu_beta_vb,
                        const double sig2_beta_vb) {

  const double c = -(log_sig2_inv_vb + log(sig2_beta_vb) )/ 2;

  for (int j = 0; j < X.cols(); ++j) {

    mat_x_m1.noalias() -= X.col(j) * m1_beta.row(j);

    mu_beta_vb.row(j) = sig2_beta_vb * ((Wy - mat_x_m1 - mat_z_mu).transpose() * X.col(j)).array();

    gam_vb.row(j) = exp(-logOnePlusExp(log_1_min_Phi_mat_v_mu.row(j) - log_Phi_mat_v_mu.row(j) -
      mu_beta_vb.row(j).square() / (2 * sig2_beta_vb) + c));

    m1_beta.row(j) = mu_beta_vb.row(j) * gam_vb.row(j);

    mat_x_m1.noalias() += X.col(j) * m1_beta.row(j);

  }

}

