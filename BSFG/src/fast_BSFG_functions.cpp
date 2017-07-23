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
    VectorXd VAiWtURinvy = VAi * WtURinvy;
    VectorXd outerWtURinvy = VAi.transpose() * inner.householderQr().solve(VAiWtURinvy);
    theta_tilda = WtURinvy.array() / prior_prec.array();
    theta_tilda -= outerWtURinvy;
  }

  VectorXd coefs = theta_tilda + theta_star;
  return coefs;
}

// [[Rcpp::export()]]
MatrixXd sample_coefs_parallel_sparse_c_Eigen(
    MSpMat Ut,
    Map<MatrixXd> Eta,
    MSpMat W,
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
    MatrixXd UtW, UtEta;
    MatrixXd prior_mean, prior_prec, randn_theta, randn_e;
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

  MatrixXd UtEta = Ut*Eta;
  MatrixXd UtW = Ut*W;

  int p = UtEta.cols();
  int b = UtW.cols();
  int n = UtW.rows();

  MatrixXd coefs(b,p);

  sampleColumn sampler(UtW, UtEta, prior_mean, prior_prec, h2, tot_Eta_prec, randn_theta,randn_e, s, b, n, coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);

  return(coefs);
}

// [[Rcpp::export()]]
Rcpp::List sample_cis_coefs_parallel_sparse_c_Eigen(
    MSpMat Ut,
    Map<MatrixXd> Eta,
    Map<MatrixXd> W,
    Rcpp::List cis_genotypes,
    Map<VectorXd> h2,
    Map<VectorXd> tot_Eta_prec,
    Map<VectorXd> s,
    Map<MatrixXd> prior_mean,
    Map<MatrixXd> prior_prec,
    Map<MatrixXd> randn_theta,
    Map<MatrixXd> randn_e,
    Map<VectorXd> randn_cis,
    Map<VectorXd> cis_effect_index,
    int grainSize){

  // Sample regression coefficients
  // columns of matrices are independent
  // each column conditional posterior is a MVN due to conjugacy

  MatrixXd UtEta = Ut*Eta;
  MatrixXd UtW = Ut*W;

  int p = UtEta.cols();
  int b = UtW.cols();
  int n = UtW.rows();

  std::vector<MatrixXd> cis_X;
  int length_cis = 0;
  for(int i = 0; i < p; i++){
    MatrixXd cXi = Rcpp::as<MatrixXd>(cis_genotypes[i]);
    cis_X.push_back(cXi);
    length_cis += cXi.cols();
  }

  MatrixXd coefs(b,p);
  VectorXd cis_effects(length_cis);

  struct sampleColumn : public Worker {
    MatrixXd UtW, UtEta;
    SpMat Ut;
    std::vector<MatrixXd> cis_X;
    MatrixXd prior_mean, prior_prec, randn_theta, randn_e;
    VectorXd randn_cis;
    VectorXd h2, tot_Eta_prec, cis_effect_index,s;
    int b,n;
    MatrixXd &coefs;
    VectorXd &cis_effects;

    sampleColumn(MatrixXd UtW, MatrixXd UtEta, SpMat Ut, std::vector<MatrixXd> cis_X,
                 MatrixXd prior_mean, MatrixXd prior_prec, VectorXd h2,
                 VectorXd tot_Eta_prec, VectorXd cis_effect_index,
                 MatrixXd randn_theta, MatrixXd randn_e,VectorXd randn_cis,
                 VectorXd s, int b, int n,
                 MatrixXd &coefs, VectorXd &cis_effects) :
      UtW(UtW), UtEta(UtEta),
      Ut(Ut), cis_X(cis_X),
      prior_mean(prior_mean), prior_prec(prior_prec),
      randn_theta(randn_theta), randn_e(randn_e), randn_cis(randn_cis),
      h2(h2), tot_Eta_prec(tot_Eta_prec), cis_effect_index(cis_effect_index),
      s(s), b(b), n(n),
      coefs(coefs), cis_effects(cis_effects) {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        int b_cis = cis_X[j].cols();
        MatrixXd UtW_cisj(n,b+b_cis);
        UtW_cisj << UtW, Ut*cis_X[j];

        VectorXd prior_mean_j = VectorXd::Zero(b+b_cis);
        prior_mean_j.head(b) = prior_mean.col(j);

        VectorXd prior_prec_j = VectorXd::Constant(b+b_cis,1e-10);
        prior_prec_j.head(b) = prior_prec.col(j);

        VectorXd randn_theta_j(b+b_cis);
        randn_theta_j.head(b) = randn_theta.col(j);
        randn_theta_j.tail(b_cis) = randn_cis.segment(cis_effect_index[j],b_cis);

        VectorXd result = sample_coefs_single(UtEta.col(j), UtW_cisj, prior_mean_j, prior_prec_j, h2(j), tot_Eta_prec(j), randn_theta_j,randn_e.col(j),s,b,n);
        coefs.col(j) = result.head(b);
        cis_effects.segment(cis_effect_index[j],b_cis) = result.tail(b_cis);
      }
    }
  };

  sampleColumn sampler(UtW, UtEta, Ut, cis_X, prior_mean, prior_prec, h2, tot_Eta_prec, cis_effect_index, randn_theta,randn_e, randn_cis, s, b, n, coefs,cis_effects);
  RcppParallel::parallelFor(0,p,sampler,grainSize);

  return(Rcpp::List::create(coefs,cis_effects));
}


// [[Rcpp::export()]]
MatrixXd sample_coefs_parallel_sparse_c_Eigen_group(
    MSpMat Ut,
    Map<MatrixXd> Eta,
    Map<MatrixXd> W,
    Map<MatrixXd> old_coefs,  // previous iteration coefficients (this matrix may change?)
    VectorXi row_groups,      // list of indexes (one-based) giving the beginnings of groups to sample jointly
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

  // Break problem into a set of samplings in a series. may be helpful to break up large matrix inversions

  struct sampleColumn : public Worker {
    MatrixXd UtW, UtEta;
    MatrixXd prior_mean, prior_prec, randn_theta, randn_e;
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

  MatrixXd UtEta = Ut*Eta;
  MatrixXd UtW = Ut*W;

  int p = UtEta.cols();
  int b = UtW.cols();
  int n = UtW.rows();

  MatrixXd coefs = old_coefs;
  int n_groups = row_groups.size();

  if(n_groups == 1){
    sampleColumn sampler(UtW, UtEta, prior_mean, prior_prec, h2, tot_Eta_prec, randn_theta,randn_e, s, b, n, coefs);
    RcppParallel::parallelFor(0,p,sampler,grainSize);
  } else{
    for(int i = 0; i < n_groups; i++){
      int start = row_groups[i]-1;
      int group_size = b - row_groups[i]+1;
      if(i < (n_groups - 1)){
        group_size = row_groups[i+1] - row_groups[i];
      }
      if(group_size <= 0) continue;
      int end = start+group_size;
      MatrixXd UtEta_resid(UtEta);
      if(start > 0){
        UtEta_resid -= UtW.block(0,0,n,start) * coefs.block(0,0,start,p);
      }
      if(end < b){
        UtEta_resid -= UtW.block(0,end,n,b-end) * coefs.block(end,0,b-end,p);
      }
      MatrixXd coefs_block(group_size,p);
      sampleColumn sampler(UtW.block(0,start,n,group_size), UtEta_resid,
                           prior_mean.block(start,0,group_size,p), prior_prec.block(start,0,group_size,p),
                           h2, tot_Eta_prec,
                           randn_theta.block(start,0,group_size,p),randn_e.block(n*i,0,n,p), s, b, n, coefs_block);
      RcppParallel::parallelFor(0,p,sampler,grainSize);
      coefs.block(start,0,group_size,p) = coefs_block;
    }
  }
  return(coefs);
}

// [[Rcpp::export()]]
List sample_coefs_set_c(
    Rcpp::List model_matrices,
    Map<VectorXd> randn_theta_vec,
    Map<VectorXd> randn_e_vec,
    Map<MatrixXd> h2s,
    Map<MatrixXd> tot_Eta_prec,
    Map<MatrixXd> prior_mean,
    Map<MatrixXd> prior_prec,
    int n,
    int grainSize){

  std::vector<MatrixXd> UtEta_list;
  std::vector<MatrixXd> UtW_list;
  std::vector<VectorXd> s_list;
  std::vector<VectorXd> position_list;
  std::vector<MatrixXd> randn_theta_list;
  std::vector<MatrixXd> randn_e_list;
  int randn_theta_index = 0;
  int randn_e_index = 0;
  int total_obs = 0;
  for(int i = 0; i < n; i++){
    Rcpp::List model_matrix_i = Rcpp::as<Rcpp::List>(model_matrices[i]);
    UtEta_list.push_back(Rcpp::as<MatrixXd>(model_matrix_i["y"]));
    UtW_list.push_back(Rcpp::as<MatrixXd>(model_matrix_i["X"]));
    position_list.push_back(Rcpp::as<VectorXd>(model_matrix_i["position"]));
    s_list.push_back(Rcpp::as<VectorXd>(model_matrix_i["s"]));

    int p = UtEta_list[i].cols();
    int n_obs = UtEta_list[i].rows();
    int b = UtW_list[i].cols();

    total_obs += n_obs;

    MatrixXd r_theta = randn_theta_vec.segment(randn_theta_index,b*p);
    Map<MatrixXd> randn_theta(r_theta.data(),b,p);
    randn_theta_list.push_back(randn_theta);
    randn_theta_index += b*p;
    MatrixXd r_e = randn_e_vec.segment(randn_e_index,n_obs*p);
    Map<MatrixXd> randn_e(r_e.data(),n_obs,p);
    randn_e_list.push_back(randn_e);
    randn_e_index += n_obs*p;
  }


  int n_traits = UtEta_list[0].cols();
  int b = randn_theta_list[0].rows()*n_traits;

  MatrixXd coefs(b,n);
  MatrixXd Y_fitted(total_obs,n_traits);

  struct sampleColumn : public RcppParallel::Worker {
    std::vector<MatrixXd> UtEta_list;
    std::vector<MatrixXd> UtW_list;
    std::vector<MatrixXd> randn_theta_list;
    std::vector<MatrixXd> randn_e_list;
    std::vector<VectorXd> position_list;
    std::vector<VectorXd> s_list;
    MatrixXd h2s,tot_Eta_prec;
    MatrixXd prior_mean,prior_prec;
    int n_traits;
    MatrixXd &coefs;
    MatrixXd &Y_fitted;

    sampleColumn(
      std::vector<MatrixXd> UtEta_list,
      std::vector<MatrixXd> UtW_list,
      std::vector<MatrixXd> randn_theta_list,
      std::vector<MatrixXd> randn_e_list,
      std::vector<VectorXd> position_list,
      std::vector<VectorXd> s_list,
      MatrixXd h2s,
      MatrixXd tot_Eta_prec,
      MatrixXd prior_mean,
      MatrixXd prior_prec,
      int n_traits,
      MatrixXd &coefs,
      MatrixXd &Y_fitted) :
      UtEta_list(UtEta_list), UtW_list(UtW_list), randn_theta_list(randn_theta_list), randn_e_list(randn_e_list),position_list(position_list),s_list(s_list),
      h2s(h2s), tot_Eta_prec(tot_Eta_prec),prior_mean(prior_mean), prior_prec(prior_prec),n_traits(n_traits),
      coefs(coefs),Y_fitted(Y_fitted)
    {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        // for(int j = 0; j < n; j++){
        int b = randn_theta_list[j].rows();
        int n_obs = randn_e_list[j].rows();
        for(int t = 0; t < n_traits; t++) {
          coefs.block(t*b,j,b,1) = sample_coefs_single(UtEta_list[j].col(t), UtW_list[j], prior_mean.block(t*b,j,b,1),
                      prior_prec.block(t*b,j,b,1), h2s(j,t), tot_Eta_prec(j,t), randn_theta_list[j].col(t),randn_e_list[j].col(t),s_list[j],b,n_obs);
        }
        Map<MatrixXd> Eta_i(coefs.col(j).data(),b,n_traits);
        MatrixXd Y_fitted_j = UtW_list[j] * Eta_i;
        for(int i = 0; i < position_list[j].size(); i++) {
          Y_fitted.row(position_list[j][i]-1) = Y_fitted_j.row(i);
        }
      }
    }
  };


  sampleColumn sampler(UtEta_list, UtW_list,randn_theta_list,randn_e_list,position_list,s_list,h2s,tot_Eta_prec,prior_mean,prior_prec,n_traits,coefs,Y_fitted);
  RcppParallel::parallelFor(0,n,sampler,grainSize);

  return(Rcpp::List::create(
      Rcpp::Named("coefs") = coefs,
      Rcpp::Named("Y_fitted") = Y_fitted));
}



// -------------------------- //
// -------------------------- //
// -------------------------- //

// [[Rcpp::export()]]
VectorXd tot_prec_scores_c (
    Map<MatrixXd> UtEta,
    Map<VectorXd> h2,
    Map<ArrayXd> s
) {

  int p = UtEta.cols();

  VectorXd scores(p);

  for(int i = 0; i < p; i++){
    ArrayXd Sigma_sqrt = sqrt(h2(i) * s + (1.0 - h2(i)));
    VectorXd SiUtEta_i = UtEta.col(i).array() / Sigma_sqrt;
    scores(i) = SiUtEta_i.dot(SiUtEta_i);
  }
  return scores;
}

// [[Rcpp::export()]]
VectorXd tot_prec_scores_withX_c (
    Map<MatrixXd> UtEta,
    Map<MatrixXd> B_F,
    Map<VectorXd> h2,
    Map<VectorXd> s,
    Map<MatrixXd> prec_B_F
) {

  int p = UtEta.cols();

  VectorXd scores(p);

  for(int i = 0; i < p; i++){
    ArrayXd Sigma_sqrt = sqrt(h2(i) * s.array() + (1.0 - h2(i)));
    VectorXd SiUtEta_i = UtEta.col(i).array() / Sigma_sqrt;
    VectorXd b_std = B_F.col(i).cwiseProduct(prec_B_F.col(i).cwiseSqrt());
    scores(i) = SiUtEta_i.dot(SiUtEta_i) + b_std.dot(b_std);
  }
  return scores;
}


// -------------------------- //
// -------------------------- //
// -------------------------- //

// [[Rcpp::export()]]
MatrixXd log_p_h2s_fast(
    Map<MatrixXd> UtEta,
    Map<VectorXd> tot_Eta_prec,
    Map<VectorXd> discrete_priors,
    Map<VectorXd> s,
    int grainSize)
{

  struct sampleColumn : public Worker {
    MatrixXd std_scores_b2;
    VectorXd discrete_priors;
    VectorXd s;
    MatrixXd &log_ps;

    sampleColumn(MatrixXd std_scores_b2,
                 VectorXd discrete_priors,
                 VectorXd s,
                 MatrixXd &log_ps):
      std_scores_b2(std_scores_b2), discrete_priors(discrete_priors), s(s), log_ps(log_ps) {}

    void operator()(std::size_t begin, std::size_t end) {
      int n = std_scores_b2.rows();
      int b = discrete_priors.size();
      for(std::size_t i = begin; i < end; i++){
        double h2 = double(i)/b;
        VectorXd s2s = h2*s.array() + (1.0-h2);
        MatrixXd std_scores2 = -0.5 * std_scores_b2 * s2s.cwiseInverse().asDiagonal();
        double det = -n/2 * log(2.0*M_PI) - 0.5*s2s.array().log().sum();
        log_ps.row(i) = std_scores2.rowwise().sum().array() + det + log(discrete_priors[i]);
      }
    }
  };

  int b = discrete_priors.size();
  int p = UtEta.cols();

  MatrixXd log_ps(b,p);

  MatrixXd std_scores_b = tot_Eta_prec.cwiseSqrt().asDiagonal() * UtEta.transpose();
  MatrixXd std_scores_b2 = std_scores_b.cwiseProduct(std_scores_b);

  sampleColumn sampler(std_scores_b2,discrete_priors,s, log_ps);
  RcppParallel::parallelFor(0,b,sampler,grainSize);
  return(log_ps);
}

// -------------------------- //
// -------------------------- //
// -------------------------- //



// [[Rcpp::export()]]
MatrixXd sample_randomEffects_parallel_sparse_c_Eigen (
    Map<MatrixXd> Eta,
    MSpMat Z,
    Map<ArrayXd> tot_prec,
    Map<ArrayXd> h2,
    List invert_aZZt_Kinv,
    Map<ArrayXXd> randn_draws,
    int grainSize) {
  //samples genetic effects on factors (F_a) conditional on the factor scores F:
  // F_i = F_{a_i} + E_i, E_i~N(0,s2*(1-h2)*I) for each latent trait i
  // U_i = zeros(r,1) if h2_i = 0
  // it is assumed that s2 = 1 because this scaling factor is absorbed in
  // Lambda
  // invert_aZZt_Kinv has parameters to diagonalize a*Z*Z' + b*K for fast
  // inversion:

  ArrayXd a_prec = tot_prec / h2;
  ArrayXd e_prec = tot_prec / (1.0-h2);

  SpMat U = as<MSpMat>(invert_aZZt_Kinv["U"]);
  VectorXd s1 = as<VectorXd>(invert_aZZt_Kinv["s1"]);
  VectorXd s2 = as<VectorXd>(invert_aZZt_Kinv["s2"]);

  int p = Eta.cols();
  int r = Z.cols();
  // MatrixXd b = U.transpose() * Z.transpose() * (Eta * e_prec.matrix().asDiagonal());
  MatrixXd b = (Z*U).transpose() * (Eta * e_prec.matrix().asDiagonal());

  // MatrixXd z = randn(r,p);

  MatrixXd effects(r,p);

  struct sampleColumn : public Worker {
    VectorXd s1, s2;
    ArrayXd a_prec, e_prec;
    SpMat U;
    MatrixXd b;
    ArrayXXd randn_draws;
    MatrixXd &effects;

    sampleColumn(VectorXd s1, VectorXd s2, ArrayXd a_prec, ArrayXd e_prec, SpMat U, MatrixXd b, ArrayXXd randn_draws, MatrixXd &effects)
      : s1(s1), s2(s2), a_prec(a_prec), e_prec(e_prec), U(U), b(b), randn_draws(randn_draws), effects(effects) {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        ArrayXd d = s2*a_prec(j) + s1*e_prec(j);
        ArrayXd mlam = b.col(j).array() / d;
        effects.col(j) = U * (mlam + randn_draws.col(j) / sqrt(d)).matrix();
      }
    }
  };

  sampleColumn sampler(s1, s2, a_prec, e_prec, U, b, randn_draws, effects);
  RcppParallel::parallelFor(0,p,sampler,grainSize);

  return(effects);
}

// -------------------------- //
// -------------------------- //
// -------------------------- //

// [[Rcpp::export()]]
MatrixXd sample_factors_scores_sparse_c_Eigen(
    Map<MatrixXd> Eta_tilde,
    Map<MatrixXd> prior_mean,
    Map<MatrixXd> Lambda,
    Map<VectorXd> resid_Eta_prec,
    Map<VectorXd> F_e_prec,
    Map<MatrixXd> randn_draws
) {
  //Sample factor scores given factor loadings (F_a), factor residual variances (F_e_prec) and
  //phenotype residuals
  MatrixXd Lmsg = resid_Eta_prec.asDiagonal() * Lambda;
  MatrixXd Sigma = Lambda.transpose() * Lmsg;
  Sigma.diagonal() += F_e_prec;
  Eigen::LLT<MatrixXd> chol_Sigma;
  chol_Sigma.compute(Sigma);
  MatrixXd R = chol_Sigma.matrixU();

  MatrixXd Meta = R.transpose().triangularView<Lower>().solve((Eta_tilde * Lmsg + prior_mean * F_e_prec.asDiagonal()).transpose());

  MatrixXd Ft = R.triangularView<Upper>().solve(Meta + randn_draws.transpose());

  return Ft.transpose();
}

