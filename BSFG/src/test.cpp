#include <math.h>
#include <iostream>
#include "BSFG_types.h"
// // // // // // // // // #include<Eigen/SparseCholesky>

// [[Rcpp::export()]]
void v1(MSpMat chol_R,MatrixXd W,VectorXd y,double tot_Eta_prec){
  SpMat chol_R_t = chol_R.transpose();
  MatrixXd RinvSqW = chol_R_t.triangularView<Lower>().solve(W);
  VectorXd WtRinvy = RinvSqW.transpose() * chol_R_t.triangularView<Lower>().solve(y) * tot_Eta_prec;
}
// [[Rcpp::export()]]
void v2(ArrayXd  resid_prec,MatrixXd X,VectorXd randn_e){
  VectorXd R_sq_diag = resid_prec.sqrt().inverse();
  MatrixXd RinvSqQtW = R_sq_diag.cwiseInverse().asDiagonal() * X;
  VectorXd eta_std = randn_e.array()/R_sq_diag.array();
  VectorXd WtURinvy = RinvSqQtW.transpose() * eta_std;
}

// [[Rcpp::export()]]
VectorXd sample_MME_single_diagK2(
    VectorXd y,
    MatrixXd W,
    VectorXd prior_mean,
    VectorXd prior_prec,
    MSpMat chol_R,
    double tot_Eta_prec,
    VectorXd randn_theta,
    VectorXd randn_e,
    int b,
    int n
){

  chol_R *= 1/tot_Eta_prec;
  VectorXd theta_star = randn_theta.array() / prior_prec.cwiseSqrt().array();
  theta_star += prior_mean;
  VectorXd e_star = chol_R * randn_e;
  MatrixXd W_theta_star = W * theta_star;
  VectorXd y_resid = y - W_theta_star - e_star;

  MatrixXd RinvSqW = chol_R.transpose().triangularView<Lower>().solve(W);
  VectorXd WtRinvy = RinvSqW.transpose() * chol_R.transpose().triangularView<Lower>().solve(y_resid);

  VectorXd theta_tilda;

  if(b < n) {
    MatrixXd C = RinvSqW.transpose() * RinvSqW;
    C.diagonal() += prior_prec;
    theta_tilda = C.llt().solve(WtRinvy);
  } else{
    MatrixXd VAi = W * prior_prec.cwiseInverse().asDiagonal();
    MatrixXd inner = VAi*W.transpose() + chol_R.transpose() * chol_R;
    VectorXd VAiWtURinvy = VAi * WtRinvy;
    VectorXd outerWtURinvy = VAi.transpose() * inner.ldlt().solve(VAiWtURinvy);
    theta_tilda = WtRinvy.array() / prior_prec.array();
    theta_tilda -= outerWtURinvy;
  }

  VectorXd theta = theta_star + theta_tilda;

  return theta;
}

// [[Rcpp::export()]]
VectorXd sample_coefs_uncorrelated2(  // returns bx1
    VectorXd y,             // nx1
    MatrixXd X,             // nxb
    VectorXd prior_mean,    // bx1
    VectorXd prior_prec,    // bx1
    ArrayXd  resid_prec,    // nx1
    VectorXd randn_theta,   // bx1
    VectorXd randn_e,       // nx1
    int b,
    int n
) {
  VectorXd R_sq_diag = resid_prec.sqrt().inverse();
  VectorXd theta_star = randn_theta.array()/prior_prec.array().sqrt();
  theta_star += prior_mean;
  VectorXd e_star = randn_e.array() * R_sq_diag.array();
  MatrixXd QtW_theta_star = X * theta_star;
  VectorXd eta_resid = y - QtW_theta_star - e_star;
  MatrixXd RinvSqQtW = R_sq_diag.cwiseInverse().asDiagonal() * X;
  VectorXd eta_std = eta_resid.array()/R_sq_diag.array();
  VectorXd WtURinvy = RinvSqQtW.transpose() * eta_std;

  VectorXd theta_tilda;
  if(b < n) {
    MatrixXd C = RinvSqQtW.transpose() * RinvSqQtW;
    C.diagonal() += prior_prec;
    theta_tilda = C.llt().solve(WtURinvy);
  } else{
    MatrixXd VAi = X * prior_prec.cwiseInverse().asDiagonal();
    MatrixXd inner = VAi*X.transpose();
    inner.diagonal() += R_sq_diag.cwiseProduct(R_sq_diag);
    VectorXd VAiWtURinvy = VAi * WtURinvy;
    VectorXd outerWtURinvy = VAi.transpose() * inner.ldlt().solve(VAiWtURinvy);
    theta_tilda = WtURinvy.array() / prior_prec.array();
    theta_tilda -= outerWtURinvy;
  }

  VectorXd coefs = theta_tilda + theta_star;
  return coefs;
}
