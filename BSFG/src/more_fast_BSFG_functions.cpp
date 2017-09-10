#include <math.h>
#include <iostream>
#include "BSFG_types.h"
#include<Eigen/SparseCholesky>


VectorXd sample_coefs_single_hierarchical(
  VectorXd QtEta,         // nx1
  SpMat QtW,              // nxr
  SpMat QtWX,             // nxp
  SpMat X,                // rxp
  VectorXd prior_mean,    // bx1
  VectorXd prior_prec,    // bx1
  ArrayXd  resid_prec,    // nx1
  VectorXd randn_theta,   // bx1
  VectorXd randn_e,       // nx1
  int b,
  int r
) {

  VectorXd R_sq_diag = resid_prec.sqrt().inverse();
  VectorXd theta_star = randn_theta.array()/prior_prec.array().sqrt();
  theta_star += prior_mean;
  VectorXd e_star = randn_e.array() * R_sq_diag.array();
  MatrixXd QtWX_theta_star = QtWX*theta_star;
  VectorXd eta_resid = QtEta - QtWX_theta_star - e_star;
  SpMat RinvSqQtWX = R_sq_diag.cwiseInverse().asDiagonal() * QtWX;
  VectorXd eta_std = eta_resid.array()/R_sq_diag.array();
  VectorXd XtWtQRinvy = RinvSqQtWX.transpose() * eta_std;

  VectorXd theta_tilda;
  if(b < r) {
    MatrixXd C = RinvSqQtWX.transpose() * RinvSqQtWX;
    C.diagonal() = C.diagonal() + prior_prec;
    theta_tilda = C.householderQr().solve(XtWtQRinvy);
  } else{
    SpMat B = QtW.transpose() * resid_prec.matrix().asDiagonal() * QtW;
    SimplicialLLT<SparseMatrix<double> >solver;
    SpMat I(r,r);
    I.setIdentity();
    SpMat Bi = solver.compute(B).solve(I);
    MatrixXd VAi = X * prior_prec.cwiseInverse().asDiagonal();
    MatrixXd inner = VAi * X.transpose() + Bi;
    VectorXd VAiXtWtQRinvy = VAi * XtWtQRinvy;
    VectorXd outerXtWtQRinvy = VAi.transpose() * inner.householderQr().solve(VAiXtWtQRinvy);
    theta_tilda = XtWtQRinvy.array() / prior_prec.array();
    theta_tilda -= outerXtWtQRinvy;
  }
  VectorXd coefs = theta_tilda + theta_star;
  return coefs;
}

// [[Rcpp::export()]]
MatrixXd sample_coefs_hierarchical_parallel_sparse_c_Eigen(
  MSpMat Qt,                  // nxn dgCMatrix
  Map<MatrixXd> Eta,          // nxp matrix
  MSpMat W,                   // nxr dgCMatrix
  MSpMat X,                   // rxb dgCMatrix
  Map<VectorXd> resid_prec,   // nxp matrix
  Map<MatrixXd> prior_mean,   // bxp matrix
  Map<MatrixXd> prior_prec,   // bxp matrix
  int grainSize){

  // Sample regression coefficients
  // columns of matrices are independent
  // each column conditional posterior is a MVN due to conjugacy


  MatrixXd QtEta = Qt*Eta;
  SpMat QtW = Qt*W;
  SpMat QtWX = QtW*X;

  int p = QtEta.cols();
  int b = X.cols();
  int n = QtW.rows();
  int r = X.rows();

  MatrixXd randn_theta = rstdnorm_mat(b,p);
  MatrixXd randn_e = rstdnorm_mat(n,p);

  MatrixXd coefs(b,p);

  struct sampleColumn : public RcppParallel::Worker {
    SpMat QtW,QtWX,X;
    MatrixXd QtEta;
    MatrixXd prior_mean, prior_prec, resid_prec, randn_theta, randn_e;
    int b,r;
    MatrixXd &coefs;

    sampleColumn(SpMat QtW, SpMat QtWX, SpMat X,MatrixXd QtEta, MatrixXd prior_mean, MatrixXd prior_prec, MatrixXd resid_prec,
                 MatrixXd randn_theta, MatrixXd randn_e,
                 int b, int r,
                 MatrixXd &coefs) :
      QtW(QtW), QtWX(QtWX), X(X),QtEta(QtEta), prior_mean(prior_mean), prior_prec(prior_prec), resid_prec(resid_prec),
      randn_theta(randn_theta), randn_e(randn_e), b(b), r(r),
      coefs(coefs) {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        // for(int j = 0; j < 1; j++){
          coefs.col(j) = sample_coefs_single_hierarchical(QtEta.col(j), QtW, QtWX, X,
                                                          prior_mean.col(j), prior_prec.col(j), resid_prec.col(j),
                                                          randn_theta.col(j),randn_e.col(j),b, r);
        }
      }
    };

    sampleColumn sampler(QtW, QtWX, X,QtEta, prior_mean, prior_prec, resid_prec, randn_theta,randn_e, b, r,coefs);
    RcppParallel::parallelFor(0,p,sampler,grainSize);

    return(coefs);
  }
