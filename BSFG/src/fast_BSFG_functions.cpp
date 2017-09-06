#include <math.h>
#include <iostream>
#include "BSFG_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace RcppParallel;

// -------------------------- //
// -------------------------- //
// -------------------------- //

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
  VectorXd theta_tilda2;
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
    theta_tilda2 = WtURinvy.array() / prior_prec.array();
    theta_tilda2 -= outerWtURinvy;
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

  MatrixXd UtEta = Ut*Eta;
  MatrixXd UtW = Ut*W;

  int p = UtEta.cols();
  int b = UtW.cols();
  int n = UtW.rows();

  MatrixXd coefs(b,p);

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
        VectorXd resid_prec = tot_Eta_prec(j) * (h2(j) * s.array() + (1.0-h2(j))).inverse();
        coefs.col(j) = sample_coefs_uncorrelated(UtEta.col(j), UtW, prior_mean.col(j), prior_prec.col(j), resid_prec, randn_theta.col(j),randn_e.col(j),b,n);
      }
    }
  };

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

        VectorXd resid_prec = tot_Eta_prec(j) * (h2(j) * s.array() + (1.0-h2(j))).inverse();
        VectorXd result = sample_coefs_uncorrelated(UtEta.col(j), UtW_cisj, prior_mean_j, prior_prec_j, resid_prec, randn_theta_j,randn_e.col(j),b,n);
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
        VectorXd resid_prec = tot_Eta_prec(j) * (h2(j) * s.array() + (1.0-h2(j))).inverse();
        coefs.col(j) = sample_coefs_uncorrelated(UtEta.col(j), UtW, prior_mean.col(j), prior_prec.col(j), resid_prec, randn_theta.col(j),randn_e.col(j),b,n);
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
          VectorXd resid_prec = tot_Eta_prec(j,t) * (h2s(j,t) * s_list[j].array() + (1.0-h2s(j,t))).inverse();
          coefs.block(t*b,j,b,1) = sample_coefs_uncorrelated(UtEta_list[j].col(t), UtW_list[j], prior_mean.block(t*b,j,b,1),
                      prior_prec.block(t*b,j,b,1), resid_prec, randn_theta_list[j].col(t),randn_e_list[j].col(t),b,n_obs);
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
      _("coefs") = coefs,
      _("Y_fitted") = Y_fitted));
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

  int b = discrete_priors.size();
  int p = UtEta.cols();

  MatrixXd log_ps(b,p);

  MatrixXd std_scores_b = tot_Eta_prec.cwiseSqrt().asDiagonal() * UtEta.transpose();
  MatrixXd std_scores_b2 = std_scores_b.cwiseProduct(std_scores_b);

  struct sampleColumn : public Worker {
    MatrixXd std_scores_b2;
    VectorXd discrete_priors;
    VectorXd tot_Eta_prec;
    VectorXd s;
    MatrixXd &log_ps;

    sampleColumn(MatrixXd std_scores_b2,
                 VectorXd discrete_priors,
                 VectorXd tot_Eta_prec,
                 VectorXd s,
                 MatrixXd &log_ps):
      std_scores_b2(std_scores_b2), discrete_priors(discrete_priors), tot_Eta_prec(tot_Eta_prec),s(s), log_ps(log_ps) {}

    void operator()(std::size_t begin, std::size_t end) {
      int n = std_scores_b2.cols();
      int b = discrete_priors.size();
      for(std::size_t i = begin; i < end; i++){
        double h2 = double(i)/b;
        VectorXd s2s = h2*s.array() + (1.0-h2);
        MatrixXd std_scores2 = -0.5 * std_scores_b2 * s2s.cwiseInverse().asDiagonal();
        double det = -n/2 * log(2.0*M_PI) - 0.5*s2s.array().log().sum();
        log_ps.row(i) = std_scores2.rowwise().sum().array() + det + log(discrete_priors[i]) + n/2*tot_Eta_prec.array().log();
      }
    }
  };

  sampleColumn sampler(std_scores_b2,discrete_priors,tot_Eta_prec,s,log_ps);
  RcppParallel::parallelFor(0,b,sampler,grainSize);
  return(log_ps);
}

double log_prob_h2_fast_c(
    ArrayXd Uty,
    ArrayXd s,
    double h2,
    int n,
    double tot_Eta_prec,
    double discrete_prior
){
  ArrayXd s2s = h2*s + (1.0-h2);
  double score2 = (Uty.pow(2)/s2s).sum() * tot_Eta_prec;
  double log_det_Sigma = s2s.log().sum();

  double log_p = -n/2.0 * log(2*M_PI) - 0.5*(log_det_Sigma - n*log(tot_Eta_prec)) - 0.5 * score2 + log(discrete_prior);
  return log_p;
}


// [[Rcpp::export()]]
VectorXi sample_h2s_discrete_MH_fast_c(
    Map<MatrixXd> UtEta,
    Map<VectorXd> tot_Eta_prec,
    Map<VectorXd> discrete_priors,
    VectorXi h2_index,
    Map<MatrixXd> h2s_matrix,
    Map<VectorXd> s,
    Map<VectorXd> r_draws,
    Map<VectorXd> state_draws,
    double step_size,
    int grainSize
){

  int p = UtEta.cols();

  struct sampleColumn : public RcppParallel::Worker {
    const MatrixXd UtEta;
    const MatrixXd h2s_matrix;
    const VectorXd s;
    const VectorXd tot_Eta_prec;
    const VectorXd discrete_priors;
    const VectorXd r_draws;
    const VectorXd state_draws;
    const VectorXi h2_index;
    const double step_size;
    VectorXi &new_index;

    sampleColumn(const MatrixXd UtEta,
                 const MatrixXd h2s_matrix,
                 const VectorXd s,
                 const VectorXd tot_Eta_prec,
                 const VectorXd discrete_priors,
                 const VectorXd r_draws,
                 const VectorXd state_draws,
                 const VectorXi h2_index,
                 const double step_size,
                 VectorXi &new_index):

      UtEta(UtEta), h2s_matrix(h2s_matrix),s(s),
      tot_Eta_prec(tot_Eta_prec),  discrete_priors(discrete_priors),
      r_draws(r_draws),state_draws(state_draws),h2_index(h2_index),step_size(step_size),new_index(new_index) {}

    void operator()(std::size_t begin, std::size_t end) {
      int n = UtEta.rows();
      for(std::size_t j = begin; j < end; j++){
        int old_state = h2_index[j] - 1;
        double old_h2 = h2s_matrix.coeffRef(0,old_state);
        double old_log_p = log_prob_h2_fast_c(UtEta.col(j),s,old_h2,n,tot_Eta_prec[j],discrete_priors[old_state]);

        VectorXd candidate_new_states = find_candidate_states(h2s_matrix,step_size,old_state);
        int r = state_draws[j] * (candidate_new_states.size());
        int proposed_state = candidate_new_states[r];

        double new_h2 = h2s_matrix.coeffRef(0,proposed_state);
        double new_log_p = log_prob_h2_fast_c(UtEta.col(j),s,new_h2,n,tot_Eta_prec[j],discrete_priors[proposed_state]);

        VectorXd candidate_states_from_new_state = find_candidate_states(h2s_matrix,step_size,proposed_state);

        double forward_prob = 1.0 / candidate_new_states.size();
        double back_prob = 1.0 / candidate_states_from_new_state.size();

        double log_MH_ratio = new_log_p - old_log_p + log(forward_prob) - log(back_prob);

        if(log(r_draws[j]) < log_MH_ratio) {
          new_index[j] = proposed_state;
        } else {
          new_index[j] = old_state;
        }
      }
    }
  };
  VectorXi new_index(p);

  sampleColumn sampler(UtEta,h2s_matrix,s,tot_Eta_prec,discrete_priors,r_draws,state_draws,h2_index,step_size,new_index);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return new_index;
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




// -------------------------- //
// ------missing data model funcitons ------- //
// -------------------------- //

// [[Rcpp::export()]]
MatrixXd sample_coefs_parallel_sparse_missing_c_Eigen2(
    Map<MatrixXd> Eta,
    Map<MatrixXd> W,
    Map<VectorXd> h2,
    Map<VectorXd> tot_Eta_prec,
    Map<MatrixXd> prior_mean,
    Map<MatrixXd> prior_prec,
    Map<MatrixXd> randn_theta,
    Map<VectorXd> randn_e,
    Rcpp::List invert_aI_bZKZ,
    VectorXi Y_obs_index,
    int grainSize){

  // Sample regression coefficients
  // columns of matrices are independent
  // each column conditional posterior is a MVN due to conjugacy

  int p = Eta.cols();
  int b = W.cols();
  MatrixXd coefs(b,p);

  std::vector<MSpMat> Ut_s;
  std::vector<VectorXd> s_s;
  std::vector<VectorXi> Y_obs;
  std::vector<MatrixXd> UtW_s;
  for(int i = 0; i < invert_aI_bZKZ.length(); i++){
    Rcpp::List invert_aI_bZKZ_i = Rcpp::as<Rcpp::List>(invert_aI_bZKZ[i]);
    Ut_s.push_back(as<MSpMat>(invert_aI_bZKZ_i["Ut"]));
    s_s.push_back(as<VectorXd>(invert_aI_bZKZ_i["s"]));
    Y_obs.push_back(as<VectorXi>(invert_aI_bZKZ_i["Y_obs"]));
    MatrixXd UtW(Y_obs[i].size(),b);
    for(int j = 0; j < Y_obs[i].size(); j++) {
      UtW.row(j) = W.row(Y_obs[i][j]-1);
    }
    UtW = Ut_s[i] * UtW;
    UtW_s.push_back(UtW);
  }

  std::vector<VectorXd> randn_e_list;
  int index = 0;
  for(int i = 0; i < Eta.cols(); i++){
    int n_obs = Y_obs[Y_obs_index[i]-1].size();
    randn_e_list.push_back(randn_e.segment(index,n_obs));
    index += n_obs;
  }

  struct sampleColumn : public RcppParallel::Worker {
    MatrixXd Eta;
    MatrixXd prior_mean, prior_prec, randn_theta;
    VectorXd h2, tot_Eta_prec;
    std::vector<MSpMat> Ut_s;
    std::vector<VectorXd> s_s;
    std::vector<VectorXi> Y_obs;
    std::vector<MatrixXd> UtW_s;
    std::vector<VectorXd> randn_e_list;
    VectorXi Y_obs_index;
    int b;
    MatrixXd &coefs;

    sampleColumn(MatrixXd Eta, MatrixXd prior_mean, MatrixXd prior_prec, VectorXd h2, VectorXd tot_Eta_prec,
                 MatrixXd randn_theta,
                 std::vector<MSpMat> Ut_s, std::vector<VectorXd> s_s, std::vector<VectorXi> Y_obs, std::vector<MatrixXd> UtW_s,
                 std::vector<VectorXd> randn_e_list,
                 VectorXi Y_obs_index, int b,
                 MatrixXd &coefs) :
      Eta(Eta), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta),
      h2(h2), tot_Eta_prec(tot_Eta_prec),
      Ut_s(Ut_s), s_s(s_s), Y_obs(Y_obs), UtW_s(UtW_s),
      randn_e_list(randn_e_list),
      Y_obs_index(Y_obs_index), b(b),
      coefs(coefs) {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        VectorXi Y_obs_j = Y_obs[Y_obs_index[j]-1];
        int n_obs = Y_obs_j.size();
        VectorXd UtEta_j(n_obs);
        for(int i = 0; i < n_obs; i++){
          UtEta_j[i] = Eta.coeffRef(Y_obs_j[i]-1,j);
        }
        UtEta_j = Ut_s[Y_obs_index[j]-1] * UtEta_j;

        VectorXd resid_prec = tot_Eta_prec(j) * (h2(j) * s_s[Y_obs_index[j]-1].array() + (1.0-h2(j))).inverse();
        coefs.col(j) = sample_coefs_uncorrelated(UtEta_j, UtW_s[Y_obs_index[j]-1], prior_mean.col(j), prior_prec.col(j), resid_prec,
                  randn_theta.col(j),randn_e_list[j],b,n_obs);
      }
    }
  };

  sampleColumn sampler(Eta, prior_mean, prior_prec, h2, tot_Eta_prec, randn_theta,Ut_s, s_s, Y_obs, UtW_s,randn_e_list,
                       Y_obs_index, b, coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);

  return(coefs);
}


// [[Rcpp::export()]]
VectorXd tot_prec_scores_missing_c (
    Map<MatrixXd> Eta,
    Map<VectorXd> h2,
    Rcpp::List invert_aI_bZKZ,
    VectorXi Y_obs_index
) {

  int p = Eta.cols();

  std::vector<MSpMat> Ut_s;
  std::vector<VectorXd> s_s;
  std::vector<VectorXi> Y_obs;
  for(int i = 0; i < invert_aI_bZKZ.length(); i++){
    Rcpp::List invert_aI_bZKZ_i = Rcpp::as<Rcpp::List>(invert_aI_bZKZ[i]);
    Ut_s.push_back(as<MSpMat>(invert_aI_bZKZ_i["Ut"]));
    s_s.push_back(as<VectorXd>(invert_aI_bZKZ_i["s"]));
    Y_obs.push_back(as<VectorXi>(invert_aI_bZKZ_i["Y_obs"]));
  }

  VectorXd scores(p);

  for(int j = 0; j < p; j++){
    ArrayXd s = s_s[Y_obs_index[j]-1];
    VectorXi Y_obs_j = Y_obs[Y_obs_index[j]-1];
    int n_obs = Y_obs_j.size();
    VectorXd UtEta_j(n_obs);
    for(int i = 0; i < n_obs; i++){
      UtEta_j[i] = Eta.coeffRef(Y_obs_j[i]-1,j);
    }
    UtEta_j = Ut_s[Y_obs_index[j]-1] * UtEta_j;
    ArrayXd Sigma_sqrt = sqrt(h2(j) * s + (1.0 - h2(j)));
    VectorXd SiUtEta_j = UtEta_j.array() / Sigma_sqrt;
    scores(j) = SiUtEta_j.dot(SiUtEta_j);
  }
  return scores;
}

// [[Rcpp::export()]]
MatrixXd log_p_h2s_fast_missing(
    Map<MatrixXd> Eta,
    Map<VectorXd> tot_Eta_prec,
    Map<VectorXd> discrete_priors,
    Rcpp::List invert_aI_bZKZ,
    VectorXi Y_obs_index,
    int grainSize)
{

  int b = discrete_priors.size();
  int p = Eta.cols();

  std::vector<MSpMat> Ut_s;
  std::vector<VectorXd> s_s;
  std::vector<VectorXi> Y_obs;
  for(int i = 0; i < invert_aI_bZKZ.length(); i++){
    Rcpp::List invert_aI_bZKZ_i = Rcpp::as<Rcpp::List>(invert_aI_bZKZ[i]);
    Ut_s.push_back(as<MSpMat>(invert_aI_bZKZ_i["Ut"]));
    s_s.push_back(as<VectorXd>(invert_aI_bZKZ_i["s"]));
    Y_obs.push_back(as<VectorXi>(invert_aI_bZKZ_i["Y_obs"]));
  }

  MatrixXd log_ps(b,p);

  std::vector<VectorXd> std_scores_b2;
  for(int j = 0; j < p; j++) {
    VectorXi Y_obs_j = Y_obs[Y_obs_index[j]-1];
    int n_obs = Y_obs_j.size();
    VectorXd Eta_std_j(n_obs);
    for(int i = 0; i < n_obs; i++){
      Eta_std_j[i] = Eta.coeffRef(Y_obs_j[i]-1,j);
    }
    VectorXd UtEta_std_j = Ut_s[Y_obs_index[j]-1] * Eta_std_j * sqrt(tot_Eta_prec[j]);
    std_scores_b2.push_back(UtEta_std_j.cwiseProduct(UtEta_std_j));
  }


  struct sampleColumn : public Worker {
    std::vector<VectorXd> std_scores_b2;
    VectorXd discrete_priors;
    VectorXd tot_Eta_prec;
    std::vector<MSpMat> Ut_s;
    std::vector<VectorXd> s_s;
    std::vector<VectorXi> Y_obs;
    VectorXi Y_obs_index;
    MatrixXd &log_ps;

    sampleColumn(std::vector<VectorXd> std_scores_b2,
                 VectorXd discrete_priors,
                 VectorXd tot_Eta_prec,
                 std::vector<MSpMat> Ut_s,
                 std::vector<VectorXd> s_s,
                 std::vector<VectorXi> Y_obs,
                 VectorXi Y_obs_index,
                 MatrixXd &log_ps):
      std_scores_b2(std_scores_b2), discrete_priors(discrete_priors), tot_Eta_prec(tot_Eta_prec),Ut_s(Ut_s), s_s(s_s), Y_obs(Y_obs),Y_obs_index(Y_obs_index), log_ps(log_ps) {}

    void operator()(std::size_t begin, std::size_t end) {
      int b = discrete_priors.size();
      int p = std_scores_b2.size();
      for(std::size_t i = begin; i < end; i++){
        double h2 = double(i)/b;
        for(int j = 0; j < p; j++) {
          int index = Y_obs_index[j] - 1;
          int n_obs = Y_obs[index].size();
          VectorXd s2s = h2*s_s[index].array() + (1.0-h2);
          VectorXd std_scores_b2_j = -0.5 * std_scores_b2[j].cwiseProduct(s2s.cwiseInverse());
          double det = -n_obs/2 * log(2.0*M_PI) - 0.5*s2s.array().log().sum();
          log_ps.coeffRef(i,j) = std_scores_b2_j.sum() + det + log(discrete_priors[i]) + n_obs/2*log(tot_Eta_prec[j]);
        }
      }
    }
  };

  sampleColumn sampler(std_scores_b2,discrete_priors,tot_Eta_prec,Ut_s,s_s,Y_obs,Y_obs_index,log_ps);
  RcppParallel::parallelFor(0,b,sampler,grainSize);
  return(log_ps);
}


// [[Rcpp::export()]]
MatrixXd sample_randomEffects_parallel_sparse_missing_c_Eigen (
    Map<MatrixXd> Eta,
    Map<ArrayXd> tot_prec,
    Map<ArrayXd> h2,
    List invert_aZZt_Kinv,
    VectorXi Y_obs_index,
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

  std::vector<MSpMat> U_s;
  std::vector<MSpMat> ZUt_s;
  std::vector<VectorXd> s1_s;
  std::vector<VectorXd> s2_s;
  std::vector<VectorXi> Y_obs;
  for(int i = 0; i < invert_aZZt_Kinv.length(); i++){
    Rcpp::List invert_aZZt_Kinv_i = Rcpp::as<Rcpp::List>(invert_aZZt_Kinv[i]);
    U_s.push_back(as<MSpMat>(invert_aZZt_Kinv_i["U"]));
    ZUt_s.push_back(as<MSpMat>(invert_aZZt_Kinv_i["ZUt"]));
    s1_s.push_back(as<VectorXd>(invert_aZZt_Kinv_i["s1"]));
    s2_s.push_back(as<VectorXd>(invert_aZZt_Kinv_i["s2"]));
    Y_obs.push_back(as<VectorXi>(invert_aZZt_Kinv_i["Y_obs"]));
  }

  int p = Eta.cols();
  int r = randn_draws.rows();
  MatrixXd b(r,p);
  for(int j = 0; j < p; j++){
    VectorXi Y_obs_j = Y_obs[Y_obs_index[j]-1];
    int n_obs = Y_obs_j.size();
    VectorXd Eta_j(n_obs);
    for(int i = 0; i < n_obs; i++){
      Eta_j[i] = Eta.coeffRef(Y_obs_j[i]-1,j);
    }
    b.col(j) = ZUt_s[Y_obs_index[j]-1] * Eta_j * e_prec[j];
  }

  MatrixXd effects(r,p);

  struct sampleColumn : public Worker {
    std::vector<MSpMat> U_s;
    std::vector<VectorXd> s1_s;
    std::vector<VectorXd> s2_s;
    VectorXi Y_obs_index;
    ArrayXd a_prec, e_prec;
    MatrixXd b;
    ArrayXXd randn_draws;
    MatrixXd &effects;

    sampleColumn(std::vector<MSpMat> U_s, std::vector<VectorXd> s1_s, std::vector<VectorXd> s2_s, VectorXi Y_obs_index,
                 ArrayXd a_prec, ArrayXd e_prec, MatrixXd b, ArrayXXd randn_draws, MatrixXd &effects):
      U_s(U_s), s1_s(s1_s), s2_s(s2_s), Y_obs_index(Y_obs_index),
      a_prec(a_prec), e_prec(e_prec), b(b), randn_draws(randn_draws), effects(effects) {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        int index = Y_obs_index[j]-1;
        ArrayXd d = s2_s[index]*a_prec(j) + s1_s[index]*e_prec(j);
        ArrayXd mlam = b.col(j).array() / d;
        effects.col(j) = U_s[index] * (mlam + randn_draws.col(j) / sqrt(d)).matrix();
      }
    }
  };

  sampleColumn sampler(U_s,s1_s,s2_s,Y_obs_index, a_prec, e_prec, b, randn_draws, effects);
  RcppParallel::parallelFor(0,p,sampler,grainSize);

  return(effects);
}



// [[Rcpp::export()]]
MatrixXd sample_factors_scores_sparse_mising_c_Eigen(
    Map<MatrixXd> Eta_tilde,
    Map<MatrixXd> prior_mean,
    Map<MatrixXd> Lambda,
    Map<VectorXd> resid_Eta_prec,
    Map<VectorXd> F_e_prec,
    Map<MatrixXd> randn_draws,
    VectorXi Y_row_obs_index,
    Rcpp::List Y_row_obs_sets
) {
  //Sample factor scores given factor loadings (F_a), factor residual variances (F_e_prec) and
  //phenotype residuals

  int k = randn_draws.rows();
  int n = randn_draws.cols();

  MatrixXd Ft(k,n);
  MatrixXd Lmsg = resid_Eta_prec.asDiagonal() * Lambda;

  std::vector<VectorXi> obs_sets;
  std::vector<MatrixXd> R_s;
  std::vector<MatrixXd> Lmsg_s;
  for(int i = 0; i < Y_row_obs_sets.length(); i++){
    VectorXi obs_set_i = as<VectorXi>(Y_row_obs_sets[i]);
    obs_sets.push_back(obs_set_i);
    int n_traits = obs_set_i.size();
    MatrixXd Lmsg_i(n_traits,k);
    MatrixXd Lambda_i(n_traits,k);
    for(int j = 0; j < n_traits; j++){
      int index = obs_set_i[j]-1;
      Lmsg_i.row(j) = Lmsg.row(index);
      Lambda_i.row(j) = Lambda.row(index);
    }
    Lmsg_s.push_back(Lmsg_i);
    MatrixXd Sigma = Lambda_i.transpose() * Lmsg_i;
    Sigma.diagonal() += F_e_prec;
    Eigen::LLT<MatrixXd> chol_Sigma;
    chol_Sigma.compute(Sigma);
    MatrixXd R = chol_Sigma.matrixU();
    R_s.push_back(R);
  }

  for(int i = 0; i < n; i++){
    int index = Y_row_obs_index[i]-1;
    VectorXi obs_set_i = obs_sets[index];
    int n_traits = obs_set_i.size();
    MatrixXd R_i = R_s[index];
    MatrixXd Lmsg_i = Lmsg_s[index];
    Eigen::RowVectorXd Eta_i(n_traits);
    for(int j = 0; j < n_traits; j++){
      int trait_index = obs_set_i[j]-1;
      Eta_i[j] = Eta_tilde.coeffRef(i,trait_index);
    }
    VectorXd Meta = R_i.transpose().triangularView<Lower>().solve((Eta_i * Lmsg_i + prior_mean.row(i) * F_e_prec.asDiagonal()).transpose());
    Ft.col(i) = R_i.triangularView<Upper>().solve(Meta + randn_draws.col(i));
  }
  return Ft.transpose();
}




/////////////////////

// [[Rcpp::export()]]
MatrixXd sample_coefMat_uncorrelated_parallel_Eigen(
    Map<MatrixXd> Y,
    Map<MatrixXd> X,
    Map<MatrixXd> resid_prec,
    Map<MatrixXd> prior_mean,
    Map<MatrixXd> prior_prec,
    Map<MatrixXd> randn_theta,
    Map<MatrixXd> randn_e,
    int grainSize){

  // Sample regression coefficients
  // columns of matrices are independent
  // each column conditional posterior is a MVN due to conjugacy

  int p = Y.cols();
  int b = X.cols();
  int n = X.rows();

  MatrixXd coefs(b,p);

  struct sampleColumn : public Worker {
    MatrixXd X, Y;
    MatrixXd resid_prec, prior_mean, prior_prec, randn_theta, randn_e;
    int b,n;
    MatrixXd &coefs;

    sampleColumn(MatrixXd X, MatrixXd Y, MatrixXd resid_prec, MatrixXd prior_mean, MatrixXd prior_prec,
                 MatrixXd randn_theta, MatrixXd randn_e,
                 int b, int n,
                 MatrixXd &coefs) :
      X(X), Y(Y), resid_prec(resid_prec), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
      b(b), n(n),
      coefs(coefs) {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        coefs.col(j) = sample_coefs_uncorrelated(Y.col(j), X, prior_mean.col(j), prior_prec.col(j), resid_prec.col(j),
                  randn_theta.col(j),randn_e.col(j),b,n);
      }
    }
  };

  sampleColumn sampler(X, Y, resid_prec, prior_mean, prior_prec, randn_theta,randn_e, b, n, coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);

  return(coefs);
}


// [[Rcpp::export()]]
MatrixXd uncorrelated_prec_mat(
    Map<VectorXd> h2,
    Map<VectorXd> tot_prec,
    Map<VectorXd> s
){
  int n = s.size();
  int p = h2.size();

  MatrixXd Prec = MatrixXd::Constant(n,p,1);
  Prec.rowwise() -= h2.transpose();
  Prec += s*h2.transpose();
  Prec = Prec.cwiseInverse() * tot_prec.asDiagonal();
  return(Prec);
}
