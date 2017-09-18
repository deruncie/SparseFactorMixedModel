#include <math.h>
#include <iostream>
#include "BSFG_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace RcppParallel;


// [[Rcpp::export()]]
MatrixXd rstdnorm_mat2(int n,int p) {  // returns nxp matrix
  VectorXd X_vec = as<VectorXd>(rnorm(n*p));
  Map<MatrixXd> X_mat(X_vec.data(),n,p);
  return(X_mat);
}

// [[Rcpp::export()]]
VectorXd sample_MME_single_diagR2(
    VectorXd y,           // nx1
    SpMat ZQt,              // nxr dgCMatrix
    SpMat Q,              // nxr dgCMatrix
    SpMat chol_S,         // rxr dgCMatrix upper triangular: chol(ZtRinvZ + diag(Kinv))
    double pe,            // double
    VectorXd randn_theta // rx1
){
  VectorXd b = ZQt * y * pe;
  b = chol_S.transpose().triangularView<Lower>().solve(b);
  b += randn_theta;
  b = Q * chol_S.triangularView<Upper>().solve(b);
  return(b);
}

// [[Rcpp::export()]]
VectorXd sample_MME_single_diagK2(  // returns b x 1 vector
    VectorXd y,           // nx1
    MatrixXd X,           // nxb
    VectorXd prior_mean,  // bx1
    VectorXd prior_prec,  // bx1
    SpMat chol_R,        // nxn upper triangular Cholesky decomposition of R. Format: dgCMatrix
    VectorXd randn_theta
){
  MatrixXd RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X);
  MatrixXd C = RinvSqX.transpose() * RinvSqX;
  C.diagonal() += prior_prec;

  VectorXd XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y);

  LLT<MatrixXd> C_llt;
  C_llt.compute(C);
  MatrixXd chol_C = C_llt.matrixU();

  VectorXd b = chol_C.transpose().triangularView<Lower>().solve(XtRinvy);
  b += randn_theta;
  b = chol_C.triangularView<Upper>().solve(b);
  return(b);
}


// RcppParallel struct for sampling coefficients from a set of parallel regression problems:
// Y = X %*% B + E, where
// b[i,j] ~ N(prior_mean[i,j],1/prior_prec[i,j])
// E[,j] ~ N(0,t(chol_R) %*% chol(R) / tot_Eta_prec[j])
// where chol_R is selected from a list based on h2s_index[j]
// Y is complete, so everythign has the same dimensions
struct sample_MME_single_diagK_worker2 : public RcppParallel::Worker {
  MatrixXd Y;
  MatrixXd X;
  MatrixXd prior_mean, prior_prec, randn_theta;
  const std::vector<MSpMat> chol_R_list;
  VectorXi h2s_index;
  VectorXd tot_Eta_prec;
  MatrixXd &coefs;

  sample_MME_single_diagK_worker2(
    MatrixXd Y,           // nxp
    MatrixXd X,           // nxb
    MatrixXd prior_mean,  // bxp
    MatrixXd prior_prec,  // bxp
    const std::vector<MSpMat> chol_R_list, // each element a MSpMat nxn upper-triangular cholesky decomposition of R
    VectorXi h2s_index,   // px1, 1-based index
    VectorXd tot_Eta_prec,// px1
    MatrixXd randn_theta, // bxp
    MatrixXd &coefs       // bxp
  ):
    Y(Y), X(X), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta),
    chol_R_list(chol_R_list), h2s_index(h2s_index), tot_Eta_prec(tot_Eta_prec), coefs(coefs) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++){
      int h2_index = h2s_index[j] - 1;
      SpMat chol_R = chol_R_list[h2_index];
      chol_R *= 1/sqrt(tot_Eta_prec[j]);
      coefs.col(j) = sample_MME_single_diagK2(Y.col(j), X, prior_mean.col(j), prior_prec.col(j), chol_R, randn_theta.col(j));
    }
  }
};


// -------------------------------------------------- //
// -- Versions of the independent residuals regression model --- //
// -------------------------------------------------- //


// basic version - all regression have same dimensions
// [[Rcpp::export()]]
MatrixXd sample_MME_fixedEffects_c2(  // returns bxp matrix
    Map<MatrixXd> Y,              // nxp
    Map<MatrixXd> X,              // nxb
    Rcpp::List Sigma_Choleskys,   // list of Cholesky decompositions of residuals. nxn upper triangular, dgCMatrix
    VectorXi h2s_index,           // px1 index of Cholesky matrix for each column
    Map<VectorXd> tot_Eta_prec,   // px1
    Map<MatrixXd> prior_mean,     // bxp
    Map<MatrixXd> prior_prec,     // bxp
    int grainSize) {

  int b = X.cols();
  int p = Y.cols();

  MatrixXd randn_theta = rstdnorm_mat2(b,p);

  std::vector<MSpMat> chol_R_list;
  for(int i = 0; i < h2s_index.maxCoeff(); i++){
    Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);
    chol_R_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
  }

  MatrixXd coefs(b,p);

  sample_MME_single_diagK_worker2 sampler(Y,X,prior_mean,prior_prec,chol_R_list,h2s_index,tot_Eta_prec,randn_theta,coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return(coefs);
}
