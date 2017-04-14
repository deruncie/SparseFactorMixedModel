#include <math.h>
#include <iostream>
#include "BSFG_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace RcppParallel;

// -------------------------- //
// -------------------------- //
// -------------------------- //

VectorXd sample_coefs_single(
    VectorXd UtEta,
    MatrixXd UtW,
    VectorXd prior_mean,
    VectorXd prior_prec,
    double h2,
    double tot_Eta_prec,
    VectorXd randn_theta,
    VectorXd randn_e,
    VectorXd s,
    int b,
    int n
) {

  VectorXd R_sq_diag = ((h2 * s.array() + (1.0-h2))/tot_Eta_prec).sqrt();
  VectorXd theta_star = randn_theta.array()/prior_prec.array().sqrt();
  theta_star += prior_mean;
  VectorXd e_star = randn_e.array() * R_sq_diag.array();
  MatrixXd UtW_theta_star = UtW * theta_star;
  VectorXd eta_resid = UtEta - UtW_theta_star - e_star;
  MatrixXd RinvSqUtW = R_sq_diag.cwiseInverse().asDiagonal() * UtW;
  VectorXd eta_std = eta_resid.array()/R_sq_diag.array();
  VectorXd WtURinvy = RinvSqUtW.transpose() * eta_std;

  VectorXd theta_tilda;
  if(b < n) {
    MatrixXd C = RinvSqUtW.transpose() * RinvSqUtW;
    C.diagonal() = C.diagonal() + prior_prec;
    theta_tilda = C.householderQr().solve(WtURinvy);
  } else{
    MatrixXd VAi = UtW * prior_prec.cwiseInverse().asDiagonal();
    MatrixXd inner = VAi*UtW.transpose();
    for(int i = 0; i < n; i++) {
      inner(i,i) += (h2 * s(i) + (1.0-h2)) / tot_Eta_prec;
    }
    MatrixXd outer = VAi.transpose() * inner.householderQr().solve(VAi);
    theta_tilda = WtURinvy.array() / prior_prec.array();
    theta_tilda -= outer * WtURinvy;
  }

  VectorXd coefs = theta_tilda + theta_star;
  return coefs;
}

// [[Rcpp::export()]]
MatrixXd sample_coefs_parallel_sparse_c(
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
      for(std::size_t j = begin; j < end; j++){
        coefs.col(j) = sample_coefs_single(UtEta.col(j), UtW, prior_mean.col(j), prior_prec.col(j), h2(j), tot_Eta_prec(j), randn_theta.col(j),randn_e.col(j),s,b,n);
      }
    }
  };

  int p = UtEta.cols();
  int b = UtW.cols();
  int n = UtW.rows();

  MatrixXd coefs(b,p);

  sampleColumn sampler(UtW, UtEta, prior_mean, prior_prec, h2, tot_Eta_prec, randn_theta,randn_e, s, b, n, coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);

  return(coefs);
}


// [[Rcpp::export()]]
MatrixXd sample_coefs_set_c(
    List model_matrices,
    List randn_draws,
    List s_vectors,
    Map<VectorXd> h2s,
    Map<VectorXd> tot_Eta_prec,
    Map<MatrixXd> prior_mean,
    Map<MatrixXd> prior_prec,
    int n,
    int grainSize){

  struct sampleColumn : public Worker {
    std::vector<VectorXd> UtEta_list;
    std::vector<MatrixXd> UtW_list;
    std::vector<VectorXd> randn_theta_list;
    std::vector<VectorXd> randn_e_list;
    std::vector<VectorXd> s_list;
    VectorXd h2s,tot_Eta_prec;
    MatrixXd prior_mean,prior_prec;
    MatrixXd &coefs;

    sampleColumn(
      std::vector<VectorXd> UtEta_list,
      std::vector<MatrixXd> UtW_list,
      std::vector<VectorXd> randn_theta_list,
      std::vector<VectorXd> randn_e_list,
      std::vector<VectorXd> s_list,
      VectorXd h2s,
      VectorXd tot_Eta_prec,
      MatrixXd prior_mean,
      MatrixXd prior_prec,
      MatrixXd &coefs) :
      UtEta_list(UtEta_list), UtW_list(UtW_list), randn_theta_list(randn_theta_list), randn_e_list(randn_e_list),s_list(s_list),
      h2s(h2s), tot_Eta_prec(tot_Eta_prec),prior_mean(prior_mean), prior_prec(prior_prec),
      coefs(coefs)
    {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        int b = randn_theta_list[j].size();
        int n = randn_e_list[j].size();
        coefs.col(j) = sample_coefs_single(UtEta_list[j], UtW_list[j], prior_mean.col(j),
                  prior_prec.col(j), h2s(j), tot_Eta_prec(j), randn_theta_list[j],
                                                                              randn_e_list[j],s_list[j],b,n);
      }
    }
  };

  std::vector<VectorXd> UtEta_list;
  std::vector<MatrixXd> UtW_list;
  std::vector<VectorXd> randn_theta_list;
  std::vector<VectorXd> randn_e_list;
  std::vector<VectorXd> s_list;
  for(int i = 0; i < n; i++){
    Rcpp::List model_matrix_i = Rcpp::as<Rcpp::List>(model_matrices[i]);
    Rcpp::List randn_draws_i = Rcpp::as<Rcpp::List>(randn_draws[i]);
    UtEta_list.push_back(Rcpp::as<VectorXd>(model_matrix_i["y"]));
    UtW_list.push_back(Rcpp::as<MatrixXd>(model_matrix_i["X"]));
    randn_theta_list.push_back(Rcpp::as<VectorXd>(randn_draws_i["randn_theta"]));
    randn_e_list.push_back(Rcpp::as<VectorXd>(randn_draws_i["randn_e"]));
    s_list.push_back(Rcpp::as<VectorXd>(s_vectors[i]));
  }

  int b = randn_theta_list[0].size();
  int p = UtEta_list.size();


  MatrixXd coefs(b,p);

  sampleColumn sampler(UtEta_list, UtW_list,randn_theta_list,randn_e_list,s_list,h2s,tot_Eta_prec,prior_mean,prior_prec,coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);

  return(coefs);
}



// -------------------------- //
// -------------------------- //
// -------------------------- //

// [[Rcpp::export()]]
VectorXd tot_prec_scores_c (
    Map<MatrixXd> UtEta,
    Map<VectorXd> h2,
    Map<VectorXd> s
) {

  int p = UtEta.cols();

  VectorXd scores(p);

  for(int i = 0; i < p; i++){
    ArrayXXd Sigma_sqrt = sqrt(h2(i) * s.array() + (1.0 - h2(i)));
    VectorXd SiUtEta_i = UtEta.col(i).array() / Sigma_sqrt;
    scores(i) = SiUtEta_i.dot(SiUtEta_i);;
  }
  return scores;
}

