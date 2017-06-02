#include <math.h>
#include <iostream>
#include "BSFG_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace RcppParallel;

// -------------------------------------------- //
// ---------- sample_MME_fixedEffects --------- //
// -------------------------------------------- //

VectorXd sample_MME_single_diagK(
    VectorXd y,
    SpMat W,
    VectorXd prior_mean,
    VectorXd prior_prec,
    SpMat chol_R,
    double tot_Eta_prec,
    VectorXd randn_theta,
    VectorXd randn_e
){

  VectorXd theta_star = randn_theta.array() / prior_prec.cwiseSqrt().array();
  theta_star += prior_mean;
  VectorXd e_star = chol_R * randn_e / sqrt(tot_Eta_prec);
  MatrixXd W_theta_star = W * theta_star;
  VectorXd y_resid = y - W_theta_star - e_star;

  MatrixXd W_mat = W;
  MatrixXd RinvSqW = chol_R.transpose().triangularView<Lower>().solve(W_mat);
  VectorXd WtRinvy = RinvSqW.transpose() * chol_R.transpose().triangularView<Lower>().solve(y_resid) * tot_Eta_prec;

  VectorXd theta_tilda;

  if(W.cols() < W.rows()) {
    MatrixXd C = RinvSqW.transpose() * RinvSqW * tot_Eta_prec;
    C.diagonal() = C.diagonal() + prior_prec;
    theta_tilda = C.householderQr().solve(WtRinvy);
  } else{
    MatrixXd R = chol_R.transpose() * chol_R / tot_Eta_prec;
    MatrixXd AiU = (W * prior_prec.cwiseInverse().asDiagonal()).transpose();
    MatrixXd R_VAiU = R + W * AiU;
    MatrixXd inner = AiU * R_VAiU.householderQr().solve(AiU.transpose());
    theta_tilda = WtRinvy.array()/prior_prec.array();
    theta_tilda -= inner * WtRinvy;
  }

  VectorXd theta = theta_star + theta_tilda;

  return theta;
}

// [[Rcpp::export]]
VectorXd sample_MME_single_diagK_c(
    Map<VectorXd> y,
    MSpMat W,
    Map<VectorXd> prior_mean,
    Map<VectorXd> prior_prec,
    MSpMat chol_R,
    double tot_Eta_prec,
    Map<VectorXd> randn_theta,
    Map<VectorXd> randn_e
){
  return sample_MME_single_diagK(y,W,prior_mean,prior_prec,chol_R,tot_Eta_prec,randn_theta,randn_e);
}

// [[Rcpp::export()]]
MatrixXd sample_MME_fixedEffects_c(
    Map<MatrixXd> Y,
    MSpMat W,
    Rcpp::List Sigma_Choleskys,
    Rcpp::IntegerVector h2s_index,
    Map<VectorXd> tot_Eta_prec,
    Map<MatrixXd> prior_mean,
    Map<MatrixXd> prior_prec,
    Map<MatrixXd> randn_theta,
    Map<MatrixXd> randn_e,
    int grainSize) {

  struct sampleColumn : public Worker {
    MatrixXd Y;
    SpMat W;
    MatrixXd prior_mean, prior_prec, randn_theta, randn_e;
    const std::vector<MSpMat> chol_R_list;
    RVector<int> h2s_index;
    VectorXd tot_Eta_prec;
    MatrixXd &coefs;

    sampleColumn(MatrixXd Y,
                 SpMat W,
                 MatrixXd prior_mean,
                 MatrixXd prior_prec,
                 const std::vector<MSpMat> chol_R_list,
                 const Rcpp::IntegerVector h2s_index,
                 VectorXd tot_Eta_prec,
                 MatrixXd randn_theta,
                 MatrixXd randn_e,
                 MatrixXd &coefs):
      Y(Y), W(W), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
      chol_R_list(chol_R_list), h2s_index(h2s_index), tot_Eta_prec(tot_Eta_prec), coefs(coefs) {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        int h2_index = h2s_index[j] - 1;
        SpMat chol_R = chol_R_list[h2_index];
        coefs.col(j) = sample_MME_single_diagK(Y.col(j), W, prior_mean.col(j), prior_prec.col(j), chol_R, tot_Eta_prec[j], randn_theta.col(j),randn_e.col(j));
      }
    }
  };

  int b = randn_theta.rows();
  int p = randn_theta.cols();

  std::vector<MSpMat> chol_R_list;
  for(int i = 0; i < max(h2s_index); i++){
    Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);
    chol_R_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
  }

  MatrixXd coefs(b,p);

  sampleColumn sampler(Y,W,prior_mean,prior_prec,chol_R_list,h2s_index,tot_Eta_prec,randn_theta,randn_e, coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return(coefs);
}

// -------------------------------------------- //
// ------------ sample_MME_ZKZts -------------- //
// -------------------------------------------- //

// [[Rcpp::export]]
VectorXd sample_MME_single_diagR(
    VectorXd Y,
    SpMat W,
    SpMat chol_C,
    double pe,
    SpMat chol_K_inv,
    double tot_Eta_prec,
    VectorXd randn_theta,
    VectorXd randn_e
){
  VectorXd theta_star = chol_K_inv.triangularView<Upper>().solve(randn_theta);
  VectorXd e_star = randn_e / sqrt(pe);
  MatrixXd W_theta_star = W * theta_star;

  VectorXd Y_resid = Y - W_theta_star - e_star;
  VectorXd WtRiy = W.transpose() * (Y_resid * pe);

  VectorXd theta_tilda = chol_C.triangularView<Upper>().solve(chol_C.transpose().triangularView<Lower>().solve(WtRiy));

  VectorXd theta = theta_tilda / tot_Eta_prec + theta_star;

  return theta;
}

// // [[Rcpp::export]]
// VectorXd sample_MME_single_diagR_c(
//     Map<VectorXd> Y,
//     MSpMat W,
//     MSpMat chol_C,
//     double pe,
//     MSpMat chol_K_inv,
//     double tot_Eta_prec,
//     Map<VectorXd> randn_theta,
//     Map<VectorXd> randn_e
// ){
//   return sample_MME_single_diagR(Y,W,chol_C,pe,chol_K_inv,tot_Eta_prec,randn_theta,randn_e);
// }

// [[Rcpp::export()]]
MatrixXd sample_MME_ZKZts_c(
    Map<MatrixXd> Y,
    MSpMat W,
    Map<VectorXd> tot_Eta_prec,
    Rcpp::List randomEffect_C_Choleskys,
    Map<MatrixXd> h2s,
    Rcpp::IntegerVector h2s_index,
    Map<MatrixXd> randn_theta,
    Map<MatrixXd> randn_e,
    int grainSize) {

  struct sampleColumn : public Worker {
    MatrixXd Y;
    SpMat W;
    const std::vector<MSpMat> chol_C_list,chol_K_inv_list;
    VectorXd pes,tot_Eta_prec;
    RVector<int> h2s_index;
    MatrixXd randn_theta, randn_e;
    MatrixXd &coefs;

    sampleColumn(MatrixXd Y,
                 SpMat W,
                 const std::vector<MSpMat> chol_C_list,
                 const std::vector<MSpMat> chol_K_inv_list,
                 VectorXd pes,
                 VectorXd tot_Eta_prec,
                 Rcpp::IntegerVector h2s_index,
                 MatrixXd randn_theta, MatrixXd randn_e,
                 MatrixXd &coefs):
      Y(Y), W(W),
      chol_C_list(chol_C_list), chol_K_inv_list(chol_K_inv_list),
      pes(pes), tot_Eta_prec(tot_Eta_prec),h2s_index(h2s_index),
      randn_theta(randn_theta), randn_e(randn_e),
      coefs(coefs) {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        int h2_index = h2s_index[j] - 1;
        SpMat chol_C = chol_C_list[h2_index];
        SpMat chol_K_inv = chol_K_inv_list[h2_index];
        chol_K_inv *= sqrt(tot_Eta_prec[j]);
        coefs.col(j) = sample_MME_single_diagR(Y.col(j), W, chol_C, pes[j],chol_K_inv, tot_Eta_prec[j], randn_theta.col(j),randn_e.col(j));
      }
    }
  };
  int b = randn_theta.rows();
  int p = randn_theta.cols();

  std::vector<MSpMat> chol_C_list,chol_K_inv_list;
  for(int i = 0; i < max(h2s_index); i++){
    Rcpp::List randomEffect_C_Cholesky_i = Rcpp::as<Rcpp::List>(randomEffect_C_Choleskys[i]);
    chol_C_list.push_back(Rcpp::as<MSpMat>(randomEffect_C_Cholesky_i["chol_C"]));
    chol_K_inv_list.push_back(Rcpp::as<MSpMat>(randomEffect_C_Cholesky_i["chol_K_inv"]));
  }

  MatrixXd coefs(b,p);
  VectorXd h2_e = 1.0 - h2s.colwise().sum().array();
  VectorXd pes = tot_Eta_prec.array() / h2_e.array();

  sampleColumn sampler(Y,W,chol_C_list,chol_K_inv_list,pes,tot_Eta_prec,h2s_index,randn_theta,randn_e,coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return(coefs);
}

// -------------------------------------------- //
// -------------- tot_prec_scores ------------- //
// -------------------------------------------- //

// [[Rcpp::export()]]
VectorXd tot_prec_scores(
    Map<MatrixXd> Y,
    Rcpp::List Sigma_Choleskys,
    Rcpp::IntegerVector h2s_index,
    int grainSize)
{

  struct sampleColumn : public Worker {
    MatrixXd Y;
    SpMat W;
    const std::vector<MSpMat> chol_Sigma_list;
    RVector<int> h2s_index;
    VectorXd &scores;

    sampleColumn(MatrixXd Y,
                 const std::vector<MSpMat> chol_Sigma_list,
                 Rcpp::IntegerVector h2s_index,
                 VectorXd &scores):
      Y(Y), chol_Sigma_list(chol_Sigma_list), h2s_index(h2s_index),
      scores(scores) {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        int h2_index = h2s_index[j] - 1;
        SpMat chol_Sigma = chol_Sigma_list[h2_index];
        VectorXd score = chol_Sigma.transpose().triangularView<Lower>().solve(Y.col(j));
        scores[j] = score.dot(score);
      }
    }
  };

  int p = Y.cols();
  VectorXd scores(p);

  std::vector<MSpMat> chol_Sigma_list;
  for(int i = 0; i < max(h2s_index); i++){
    Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);
    chol_Sigma_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
  }

  sampleColumn sampler(Y,chol_Sigma_list,h2s_index,scores);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return scores;
}

// -------------------------------------------- //
// ---------------- sample h2s ---------------- //
// -------------------------------------------- //

// [[Rcpp::export()]]
MatrixXd log_p_h2s(
    Map<MatrixXd> Y,
    Map<VectorXd> tot_Eta_prec,
    Rcpp::List Sigma_Choleskys,
    Map<VectorXd> discrete_priors,
    int grainSize)
{

  struct sampleColumn : public Worker {
    MatrixXd Y;
    VectorXd tot_Eta_prec;
    const std::vector<MSpMat> chol_Sigma_list;
    VectorXd log_det_Sigmas;
    VectorXd discrete_priors;
    MatrixXd &log_ps;

    sampleColumn(MatrixXd Y,
                 VectorXd tot_Eta_prec,
                 const std::vector<MSpMat> chol_Sigma_list,
                 VectorXd log_det_Sigmas,
                 VectorXd discrete_priors,
                 MatrixXd &log_ps):
      Y(Y), tot_Eta_prec(tot_Eta_prec), chol_Sigma_list(chol_Sigma_list),
      log_det_Sigmas(log_det_Sigmas), discrete_priors(discrete_priors), log_ps(log_ps) {}

    void operator()(std::size_t begin, std::size_t end) {
      int p = Y.cols();
      int n = Y.rows();
      for(std::size_t i = begin; i < end; i++){
        SpMat chol_Sigma = chol_Sigma_list[i];
        VectorXd scores2(p);
        for(int j = 0; j < p; j++){
          VectorXd x_std = chol_Sigma.transpose().triangularView<Lower>().solve(Y.col(j));
          scores2[j] = tot_Eta_prec[j] * x_std.dot(x_std);
        }
        log_ps.row(i) = (-n/2.0 * log(2*M_PI) - 0.5 * (log_det_Sigmas[i] - n*tot_Eta_prec.array().log()) -
          0.5 * scores2.array() + log(discrete_priors[i]));
      }
    }
  };

  int b = discrete_priors.size();
  int p = Y.cols();

  std::vector<MSpMat> chol_Sigma_list;
  VectorXd log_det_Sigmas(b);
  for(int i = 0; i < b; i++){
    Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);
    chol_Sigma_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
    log_det_Sigmas[i] = Rcpp::as<double>(Sigma_Choleskys_i["log_det"]);
  }

  MatrixXd log_ps(b,p);

  sampleColumn sampler(Y,tot_Eta_prec,chol_Sigma_list,log_det_Sigmas,discrete_priors,log_ps);
  RcppParallel::parallelFor(0,b,sampler,grainSize);
  return(log_ps);
}


// [[Rcpp::export()]]
Rcpp::IntegerVector sample_h2s(
    Map<ArrayXXd> log_ps,
    Map<VectorXd> rs,
    int grainSize
)
{

  struct sampleColumn : public Worker {
    ArrayXXd log_ps;
    VectorXd rs;
    RVector<int> h2s_index;

    sampleColumn(ArrayXXd log_ps,
                 VectorXd rs,
                 Rcpp::IntegerVector h2s_index):
      log_ps(log_ps), rs(rs), h2s_index(h2s_index) {}

    void operator()(std::size_t begin, std::size_t end) {
      int b = log_ps.rows();
      for(std::size_t j = begin; j < end; j++){
        double max_col = log_ps.col(j).maxCoeff();
        double norm_factor = max_col + log((log_ps.col(j) - max_col).exp().sum());
        VectorXd ps_j = (log_ps.col(j) - norm_factor).exp();
        h2s_index[j] = 1;
        double cumsum = 0;
        for(int i = 0; i < b; i++){
          cumsum += ps_j[i];
          if(rs[j] > cumsum) {
            h2s_index[j] ++;
          }
        }
      }
    }
  };
  int p = log_ps.cols();

  Rcpp::IntegerVector h2s_index(p);

  sampleColumn sampler(log_ps,rs,h2s_index);
  RcppParallel::parallelFor(0,p,sampler,grainSize);

  return h2s_index;
}


// -------------------------------------------- //
// --------------- sample h2s MH -------------- //
// -------------------------------------------- //

double log_prob_h2_c(
    VectorXd y,
    SpMat chol_Sigma,
    double log_det_Sigma,
    int n,
    double tot_Eta_prec,
    double discrete_prior
){
  VectorXd x_std = chol_Sigma.transpose().triangularView<Lower>().solve(y);
  double score2 = tot_Eta_prec * x_std.dot(x_std);

  double log_p = -n/2.0 * log(2*M_PI) - 0.5*(log_det_Sigma - n*log(tot_Eta_prec)) - 0.5 * score2 + log(discrete_prior);
  return log_p;
}

// [[Rcpp::export()]]
VectorXd find_candidate_states(
    MatrixXd h2s_matrix,
    double step_size,
    int old_state
) {
  VectorXd dists = (h2s_matrix.colwise() - h2s_matrix.col(old_state)).cwiseAbs().colwise().sum();
  VectorXd indices(dists.size());
  int count = 0;
  for(int i = 0; i < dists.size(); i++){
    if(dists[i] < step_size & dists[i] > 0) {
      indices[count] = i;
      count++;
    }
  }
  if(count == 0) {  // return all indices as candidates
    for(int i = 0; i < dists.size(); i++){
      indices[count] = i;
      count++;
    }
  }
  return indices.head(count);
}

// [[Rcpp::export()]]
VectorXi sample_h2s_discrete_MH_c(
    Map<MatrixXd> Y,
    Map<VectorXd> tot_Eta_prec,
    Map<VectorXd> discrete_priors,
    VectorXi h2_index,
    Map<MatrixXd> h2s_matrix,
    Rcpp::List Sigma_Choleskys,
    Map<VectorXd> r_draws,
    Map<VectorXd> state_draws,
    double step_size,
    int grainSize
){

  struct sampleColumn : public Worker {
    const MatrixXd Y;
    const MatrixXd h2s_matrix;
    const std::vector<MSpMat> chol_Sigma_list;
    const VectorXd log_det_Sigmas;
    const VectorXd tot_Eta_prec;
    const VectorXd discrete_priors;
    const VectorXd r_draws;
    const VectorXd state_draws;
    const VectorXi h2_index;
    const double step_size;
    VectorXi &new_index;

    sampleColumn(const MatrixXd Y,
                 const MatrixXd h2s_matrix,
                 const std::vector<MSpMat> chol_Sigma_list,
                 const VectorXd log_det_Sigmas,
                 const VectorXd tot_Eta_prec,
                 const VectorXd discrete_priors,
                 const VectorXd r_draws,
                 const VectorXd state_draws,
                 const VectorXi h2_index,
                 const double step_size,
                 VectorXi &new_index):

      Y(Y), h2s_matrix(h2s_matrix),
      chol_Sigma_list(chol_Sigma_list), log_det_Sigmas(log_det_Sigmas),
      tot_Eta_prec(tot_Eta_prec),  discrete_priors(discrete_priors),
      r_draws(r_draws),state_draws(state_draws),h2_index(h2_index),step_size(step_size),new_index(new_index) {}

    void operator()(std::size_t begin, std::size_t end) {
      int n = Y.rows();
      for(std::size_t j = begin; j < end; j++){
        int old_state = h2_index[j] - 1;
        double old_log_p = log_prob_h2_c(Y.col(j),chol_Sigma_list[old_state],log_det_Sigmas[old_state],n,tot_Eta_prec[j],discrete_priors[old_state]);

        VectorXd candidate_new_states = find_candidate_states(h2s_matrix,step_size,old_state);
        int r = state_draws[j] * (candidate_new_states.size());
        int proposed_state = candidate_new_states[r];

        double new_log_p = log_prob_h2_c(Y.col(j),chol_Sigma_list[proposed_state],log_det_Sigmas[proposed_state],n,tot_Eta_prec[j],discrete_priors[proposed_state]);

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

  int p = Y.cols();
  int b = discrete_priors.size();

  std::vector<MSpMat> chol_Sigma_list;
  VectorXd log_det_Sigmas(b);
  for(int i = 0; i < b; i++){
    Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);

    chol_Sigma_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
    log_det_Sigmas[i] = Rcpp::as<double>(Sigma_Choleskys_i["log_det"]);
  }

  VectorXi new_index(p);

  sampleColumn sampler(Y,h2s_matrix,chol_Sigma_list,log_det_Sigmas,tot_Eta_prec,discrete_priors,r_draws,state_draws,h2_index,step_size,new_index);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return new_index;
}
