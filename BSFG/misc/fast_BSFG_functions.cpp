#include <math.h>
#include <iostream>
#include "BSFG_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace RcppParallel;


// [[Rcpp::export()]]
MatrixXd sample_coefs_parallel_sparse_c_2(
    Map<MatrixXd> UtEta,
    Map<MatrixXd> UtW,
    Map<VectorXd> h2,
    Map<VectorXd> tot_Eta_prec,
    Map<VectorXd> s,
    Map<MatrixXd> prior_mean,
    Map<MatrixXd> prior_prec,
    Map<MatrixXd> randn_theta,
    Map<MatrixXd> randn_e,
    int grainSize){

  // Sample regression coefficients
  // columns of matrices are independent
  // each column conditional posterior is a MVN due to conjugacy


  struct sampleColumn : public Worker {
    MatrixXd UtW, UtEta, prior_mean, prior_prec, randn_theta, randn_e;
    VectorXd h2, tot_Eta_prec, s;
    int b,n;
    MatrixXd &coefs;

    sampleColumn(MatrixXd UtW, MatrixXd UtEta, MatrixXd prior_mean, MatrixXd prior_prec, VectorXd h2, VectorXd tot_Eta_prec,
                 MatrixXd randn_theta, MatrixXd randn_e,
                 VectorXd s, int b, int n,
                 MatrixXd &coefs) :
      UtW(UtW), UtEta(UtEta), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
      h2(h2), tot_Eta_prec(tot_Eta_prec),
      s(s), b(b), n(n),
      coefs(coefs) {}

    void operator()(std::size_t begin, std::size_t end) {
      MatrixXd RinvSqUtW, WtURinvy;
      VectorXd R_sq_diag, theta_star, e_star, UtW_theta_star, eta_resid, theta_tilda;
      for(std::size_t j = begin; j < end; j++){
        double tot_Eta_prec_j = tot_Eta_prec[j];
        R_sq_diag = ((h2[j] * s.array() + (1-h2[j]))/tot_Eta_prec_j).sqrt();
        theta_star = randn_theta.col(j).array()/prior_prec.col(j).array().sqrt() + prior_mean.col(j).array();
        e_star = randn_e.col(j).array() * R_sq_diag.array();
        UtW_theta_star = UtW * theta_star;
        eta_resid = UtEta.col(j) - UtW_theta_star - e_star;
        RinvSqUtW = UtW * R_sq_diag.cwiseInverse().asDiagonal();
        VectorXd WtURinvy = RinvSqUtW.transpose() * (eta_resid.array()/R_sq_diag.array()).matrix();

        if(b < n) {
          MatrixXd C = RinvSqUtW.transpose() * RinvSqUtW * tot_Eta_prec;
          C.diagonal() = C.diagonal() + prior_prec;
          theta_tilda = C.householderQr().solve(WtURinvy);
        } else{
          MatrixXd VAi = UtW * prior_prec.col(j).cwiseInverse().asDiagonal();
          MatrixXd inner = VAi*UtW.transpose();
          for(int i = 0; i < n; i++) {
            inner(i,i) += (h2(j) * s(i) + (1-h2(j))) / tot_Eta_prec(j);
          }
          MatrixXd outer = VAi.transpose() * inner.householderQr().solve(VAi);
          theta_tilda = WtURinvy.array() / prior_prec.col(j).array();
          theta_tilda -= outer * WtURinvy;
        }

        coefs.col(j) = theta_tilda + theta_star;
      }
    }
  };

  int p = tot_Eta_prec.size();
  int b = UtW.cols();
  int n = UtW.rows();

  MatrixXd coefs(b,p);

  sampleColumn sampler(UtW, UtEta, prior_mean, prior_prec, h2, tot_Eta_prec, randn_theta,randn_e, s, b, n, coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return(coefs);
}
