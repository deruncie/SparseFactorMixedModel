#include <math.h>
#include <iostream>
#include "BSFG_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace RcppParallel;


// [[Rcpp::export()]]
VectorXd sample_delta_c_Eigen(
    VectorXd delta,
    VectorXd tauh,
    Map<VectorXd> scores,
    double delta_1_rate,
    double delta_2_rate,
    Map<MatrixXd> randg_draws  // all done with rate = 1;
) {
  int times = randg_draws.rows();
  int k = tauh.size();

  double rate,delta_old;
  for(int i = 0; i < times; i++){
    delta_old = delta(0);
    rate = delta_1_rate + (1/delta(0)) * tauh.dot(scores);
    delta(0) = randg_draws(i,0) / rate;
    // tauh = cumprod(delta);
    tauh *= delta(0)/delta_old;   // replaces re-calculating cumprod

    for(int h = 1; h < k; h++) {
      delta_old = delta(h);
      rate = delta_2_rate + (1/delta(h))*tauh.tail(k-h).dot(scores.tail(k-h));
      delta(h) = randg_draws(i,h) / rate;
      // tauh = cumprod(delta);
      tauh.tail(k-h) *= delta(h)/delta_old; // replaces re-calculating cumprod
      // Rcout << (tauh - cumprod(delta)).sum() << std::endl;
    }
  }
  return(delta);
}


// -------------------------- //
// -------------------------- //
// -------------------------- //

// [[Rcpp::export()]]
NumericVector rgig_multiple(
  int n,
  NumericVector lambda,
  NumericVector chi,
  NumericVector psi
){
  NumericVector res(n);
  SEXP (*fun)(SEXP,SEXP,SEXP,SEXP) = NULL;
  if (!fun)
    fun = (SEXP(*)(SEXP,SEXP,SEXP,SEXP)) R_GetCCallable("GIGrvg", "rgig");
  for(int i = 0; i < n; i++){
    res[i] = as<double>(fun(wrap(1.0),wrap(lambda[i]),wrap(chi[i]),wrap(psi[i])));
  }
  return(res);
}


// -------------------------- //
// ----Y = ZXB + E ---------- //
// -------------------------- //

VectorXd sample_MME_single_hierarchical_diagK(  // returns b x 1 vector
    VectorXd y,           // nx1
    MSpMat   Z,           // nxr   // use when r < n < b
    MatrixXd X,           // rxb
    VectorXd prior_mean,  // bx1
    VectorXd prior_prec,  // bx1
    MSpMat chol_R,        // nxn upper triangular Cholesky decomposition of R. Format: dgCMatrix
    double tot_Eta_prec, // double
    VectorXd randn_theta, // bx1
    VectorXd randn_e      // 0x1 or nx1. 0x1 if b<n
){
  VectorXd theta_star = randn_theta.array() / prior_prec.cwiseSqrt().array();
  theta_star += prior_mean;
  VectorXd e_star = chol_R.transpose() * (randn_e / sqrt(tot_Eta_prec));
  MatrixXd ZX = Z*X;
  MatrixXd ZX_theta_star = ZX * theta_star;
  VectorXd y_resid = y - ZX_theta_star - e_star;

  MatrixXd RinvSqZ = chol_R.transpose().triangularView<Lower>().solve(Z.toDense() * sqrt(tot_Eta_prec));
  MatrixXd RinvSqZX = RinvSqZ*X;
  VectorXd XtZtRinvy = RinvSqZX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y_resid * sqrt(tot_Eta_prec));

  VectorXd theta_tilda;
  MatrixXd VAi = X * prior_prec.cwiseInverse().asDiagonal();
  MatrixXd I(Z.cols(),Z.cols());
  I.setIdentity();
  MatrixXd ZtRinvZ = RinvSqZ.transpose()*RinvSqZ;
  MatrixXd inner = VAi*X.transpose() + ZtRinvZ.ldlt().solve(I);
  VectorXd VAiXtURinvy = VAi * XtZtRinvy;
  VectorXd outerXtURinvy = VAi.transpose() * inner.ldlt().solve(VAiXtURinvy);
  theta_tilda = XtZtRinvy.array() / prior_prec.array();
  theta_tilda -= outerXtURinvy;
  VectorXd theta = theta_star + theta_tilda;

  return(theta);
}

struct sample_MME_single_hierarchical_diagK_worker : public RcppParallel::Worker {
  MatrixXd Y;
  MSpMat Z;
  MatrixXd X;
  MatrixXd prior_mean, prior_prec, randn_theta, randn_e;
  const std::vector<MSpMat> chol_R_list;
  VectorXi h2s_index;
  VectorXd tot_Eta_prec;
  MatrixXd &coefs;

  sample_MME_single_hierarchical_diagK_worker(
    MatrixXd Y,           // nxp
    MSpMat   Z,           // nxb
    MatrixXd X,           // nxb
    MatrixXd prior_mean,  // bxp
    MatrixXd prior_prec,  // bxp
    const std::vector<MSpMat> &chol_R_list, // each element a MSpMat nxn upper-triangular cholesky decomposition of R
    VectorXi h2s_index,   // px1, 1-based index
    VectorXd tot_Eta_prec,// px1
    MatrixXd randn_theta, // bxp
    MatrixXd randn_e,     // 0xp or nxp. The former if b<n, the latter if b >= n
    MatrixXd &coefs       // bxp
  ):
    Y(Y), Z(Z), X(X), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
    chol_R_list(chol_R_list), h2s_index(h2s_index), tot_Eta_prec(tot_Eta_prec), coefs(coefs) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++){
      int h2_index = h2s_index[j] - 1;
      MSpMat chol_R = chol_R_list[h2_index];
      coefs.col(j) = sample_MME_single_hierarchical_diagK(Y.col(j), Z, X, prior_mean.col(j), prior_prec.col(j), chol_R, tot_Eta_prec[j], randn_theta.col(j),randn_e.col(j));
    }
  }
};

// Samples from B in model:
// Y = ZXB + E
// [[Rcpp::export()]]
MatrixXd sample_MME_fixedEffects_hierarchical_c(  // returns bxp matrix
    Map<MatrixXd> Y,              // nxp
    MSpMat        Z,              // nxr   // use when r < n < b
    Map<MatrixXd> X,              // rxb
    Rcpp::List Sigma_Choleskys,   // list of Cholesky decompositions of residuals. nxn upper triangular, dgCMatrix
    VectorXi h2s_index,           // px1 index of Cholesky matrix for each column
    Map<VectorXd> tot_Eta_prec,   // px1
    Map<MatrixXd> prior_mean,     // bxp
    Map<MatrixXd> prior_prec,     // bxp
    int grainSize) {

  int b = X.cols();
  int p = Y.cols();
  int n = Y.rows();

  MatrixXd randn_theta = rstdnorm_mat(b,p);
  MatrixXd randn_e = rstdnorm_mat(n,p);

  std::vector<MSpMat> chol_R_list;
  for(int i = 0; i < h2s_index.maxCoeff(); i++){
    Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);
    chol_R_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
  }

  MatrixXd coefs(b,p);

  sample_MME_single_hierarchical_diagK_worker sampler(Y,Z,X,prior_mean,prior_prec,chol_R_list,h2s_index,tot_Eta_prec,randn_theta,randn_e, coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);

  return(coefs);
}


