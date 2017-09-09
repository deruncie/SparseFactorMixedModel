#include <math.h>
#include <iostream>
#include "BSFG_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace RcppParallel;


// -------------------------- //
// -------------------------- //
// -------------------------- //

// functions to speed up sparse multiplication and conversion to dense matrices

// [[Rcpp::export()]]
MatrixXd SxD(MSpMat X, Map<MatrixXd> Y){
  return(X * Y);
}
// [[Rcpp::export()]]
MatrixXd SxS(MSpMat X, MSpMat Y){
  return(X * Y);
}


/*** R
`%**%` = function(X1,X2){
  if(inherits(X1,'dgCMatrix') && inherits(X2,'matrix')) return(SxD(X1,X2))
    if(inherits(X1,'dgCMatrix') && inherits(X2,'dgCMatrix')) return(SxS(X1,X2))
      if(inherits(X1,'matrix') & inherits(X2,'matrix')) return(X1 %*% X2)
        return(as.matrix(X1 %*% X2))
}
*/


// -------------------------- //
// -------------------------- //
// -------------------------- //

// basic functions for

// draws a sample of the vector b from the model:
// y = X %*% beta + e
// with beta[i] ~ N(prior_mean[i],1/prior_prec[i]), i=1:b
// with e[j] ~ N(0,1/resid_prec[i]), i=1:n
// Uses sampling method from MCMCglmm, which requires draws b+n draws from N(0,1), which are passed as randn_theta and randn_e
// If b >= n, inverts the C matrix using Binomial Inverse Theorem
VectorXd sample_coefs_uncorrelated(
    VectorXd y,
    MatrixXd X,
    VectorXd prior_mean,
    VectorXd prior_prec,
    ArrayXd  resid_prec,
    VectorXd randn_theta,
    VectorXd randn_e,
    int b,
    int n
) {
  VectorXd R_sq_diag = resid_prec.sqrt().inverse();
  VectorXd theta_star = randn_theta.array()/prior_prec.array().sqrt();
  theta_star += prior_mean;
  VectorXd e_star = randn_e.array() * R_sq_diag.array();
  MatrixXd UtW_theta_star = X * theta_star;
  VectorXd eta_resid = y - UtW_theta_star - e_star;
  MatrixXd RinvSqUtW = R_sq_diag.cwiseInverse().asDiagonal() * X;
  VectorXd eta_std = eta_resid.array()/R_sq_diag.array();
  VectorXd WtURinvy = RinvSqUtW.transpose() * eta_std;

  VectorXd theta_tilda;
  if(b < n) {
    MatrixXd C = RinvSqUtW.transpose() * RinvSqUtW;
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
