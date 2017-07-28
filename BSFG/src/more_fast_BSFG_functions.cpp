#include <math.h>
#include <iostream>
#include "BSFG_types.h"
#include<Eigen/SparseCholesky>


VectorXd sample_coefs_single_hierarchical(
  VectorXd UtEta,
  SpMat UtW,
  SpMat UtWX,
  SpMat X,
  VectorXd prior_mean,
  VectorXd prior_prec,
  double h2,
  double tot_Eta_prec,
  VectorXd randn_theta,
  VectorXd randn_e,
  VectorXd s,
  int b,
  int n,
  int r
) {

  VectorXd R_sq_diag = ((h2 * s.array() + (1.0-h2))/tot_Eta_prec).sqrt();
  VectorXd theta_star = randn_theta.array()/prior_prec.array().sqrt();
  theta_star += prior_mean;
  VectorXd e_star = randn_e.array() * R_sq_diag.array();
  MatrixXd UtWX_theta_star = UtWX*theta_star;
  VectorXd eta_resid = UtEta - UtWX_theta_star - e_star;
  SpMat RinvSqUtWX = R_sq_diag.cwiseInverse().asDiagonal() * UtWX;
  VectorXd eta_std = eta_resid.array()/R_sq_diag.array();
  VectorXd XtWtURinvy = RinvSqUtWX.transpose() * eta_std;

  VectorXd theta_tilda;
  if(b < r) {
    MatrixXd C = RinvSqUtWX.transpose() * RinvSqUtWX;
    C.diagonal() = C.diagonal() + prior_prec;
    theta_tilda = C.householderQr().solve(XtWtURinvy);
  } else{
    SpMat B = UtW.transpose() * R_sq_diag.array().pow(2).inverse().matrix().asDiagonal() * UtW;
    SimplicialLLT<SparseMatrix<double> >solver;
    SpMat I(r,r);
    I.setIdentity();
    SpMat Bi = solver.compute(B).solve(I);
    MatrixXd VAi = X * prior_prec.cwiseInverse().asDiagonal();
    MatrixXd inner = VAi * X.transpose() + Bi;
    VectorXd VAiXtWtURinvy = VAi * XtWtURinvy;
    VectorXd outerXtWtURinvy = VAi.transpose() * inner.householderQr().solve(VAiXtWtURinvy);
    theta_tilda = XtWtURinvy.array() / prior_prec.array();
    theta_tilda -= outerXtWtURinvy;
  }
  VectorXd coefs = theta_tilda + theta_star;
  return coefs;
}

// [[Rcpp::export()]]
MatrixXd sample_coefs_hierarchical_parallel_sparse_c_Eigen(
  MSpMat Ut,
  Map<MatrixXd> Eta,
  MSpMat W,
  MSpMat X,
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


  MatrixXd UtEta = Ut*Eta;
  SpMat UtW = Ut*W;
  SpMat UtWX = UtW*X;

  int p = UtEta.cols();
  int b = X.cols();
  int n = UtW.rows();
  int r = X.rows();

  MatrixXd coefs(b,p);

  struct sampleColumn : public RcppParallel::Worker {
    SpMat UtW,UtWX,X;
    MatrixXd UtEta;
    MatrixXd prior_mean, prior_prec, randn_theta, randn_e;
    VectorXd h2, tot_Eta_prec, s;
    int b,n,r;
    MatrixXd &coefs;

    sampleColumn(SpMat UtW, SpMat UtWX, SpMat X,MatrixXd UtEta, MatrixXd prior_mean, MatrixXd prior_prec, VectorXd h2, VectorXd tot_Eta_prec,
                 MatrixXd randn_theta, MatrixXd randn_e,
                 VectorXd s, int b, int n,int r,
                 MatrixXd &coefs) :
      UtW(UtW), UtWX(UtWX), X(X),UtEta(UtEta), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
    h2(h2), tot_Eta_prec(tot_Eta_prec),
    s(s), b(b), n(n),r(r),
    coefs(coefs) {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        // for(int j = 0; j < 1; j++){
          coefs.col(j) = sample_coefs_single_hierarchical(UtEta.col(j), UtW, UtWX, X,prior_mean.col(j), prior_prec.col(j), h2(j), tot_Eta_prec(j), randn_theta.col(j),randn_e.col(j),s,b,n,r);
        }
      }
    };

    sampleColumn sampler(UtW, UtWX, X,UtEta, prior_mean, prior_prec, h2, tot_Eta_prec, randn_theta,randn_e, s, b, n, r,coefs);
    RcppParallel::parallelFor(0,p,sampler,grainSize);

    return(coefs);
  }
