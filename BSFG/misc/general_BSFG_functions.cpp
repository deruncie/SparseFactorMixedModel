
// #include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <RcppParallel.h>
#include <math.h>

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace RcppParallel;

using Eigen::Map;               	      // 'Eigen::Maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::ArrayXXd;                  // variable size matrix, double precision
using Eigen::Upper;
using Eigen::Lower;
typedef Eigen::MappedSparseMatrix<double> MSpMat;
typedef Eigen::SparseMatrix<double> Eigen_SpMat;


// -------------------------------------------- //
// ---------- sample_MME_fixedEffects --------- //
// -------------------------------------------- //

Eigen::VectorXd sample_MME_single_diagK(
    Eigen::VectorXd y,
    Eigen::MatrixXd W,
    Eigen::VectorXd prior_mean,
    Eigen::VectorXd prior_prec,
    Eigen_SpMat chol_R,
    double tot_Eta_prec,
    Eigen::VectorXd randn_theta,
    Eigen::VectorXd randn_e
){

  Eigen::VectorXd theta_star = randn_theta.array() / prior_prec.cwiseSqrt().array();
  theta_star += prior_mean;
  Eigen::VectorXd e_star = chol_R * randn_e / sqrt(tot_Eta_prec);
  Eigen::MatrixXd W_theta_star = W * theta_star;
  Eigen::VectorXd y_resid = y - W_theta_star - e_star;

  Eigen::MatrixXd RinvSqW = chol_R.triangularView<Lower>().solve(W);
  Eigen::VectorXd WtRinvy = RinvSqW.transpose() * chol_R.triangularView<Lower>().solve(y_resid) * tot_Eta_prec;

  Eigen::VectorXd theta_tilda;

  if(W.cols() < W.rows()) {
    Eigen::MatrixXd C = RinvSqW.transpose() * RinvSqW * tot_Eta_prec;
    C.diagonal() = C.diagonal() + prior_prec;
    theta_tilda = C.householderQr().solve(WtRinvy);
  } else{
    Eigen::MatrixXd R = chol_R * chol_R.transpose() / tot_Eta_prec;
    Eigen::MatrixXd AiU = (W * prior_prec.cwiseInverse().asDiagonal()).transpose();
    Eigen::MatrixXd R_VAiU = R + W * AiU;
    Eigen::MatrixXd inner = AiU * R_VAiU.householderQr().solve(AiU.transpose());
    theta_tilda = WtRinvy.array()/prior_prec.array();
    theta_tilda -= inner * WtRinvy;
  }

  Eigen::VectorXd theta = theta_star + theta_tilda;

  return theta;
}

// [[Rcpp::export]]
Eigen::VectorXd sample_MME_single_diagK_c(
    Eigen::Map<Eigen::VectorXd> y,
    Eigen::Map<Eigen::MatrixXd> W,
    Eigen::Map<Eigen::VectorXd> prior_mean,
    Eigen::Map<Eigen::VectorXd> prior_prec,
    MSpMat chol_R,
    double tot_Eta_prec,
    Eigen::Map<Eigen::VectorXd> randn_theta,
    Eigen::Map<Eigen::VectorXd> randn_e
){
  return sample_MME_single_diagK(y,W,prior_mean,prior_prec,chol_R,tot_Eta_prec,randn_theta,randn_e);
}

// [[Rcpp::export()]]
Eigen::MatrixXd sample_MME_fixedEffects_c(
    Eigen::Map<Eigen::MatrixXd> Y,
    Eigen::Map<Eigen::MatrixXd> W,
    Rcpp::List chol_Rs,
    Rcpp::IntegerVector h2s_index,
    Eigen::Map<Eigen::VectorXd> tot_Eta_prec,
    Eigen::Map<Eigen::MatrixXd> prior_mean,
    Eigen::Map<Eigen::MatrixXd> prior_prec,
    Eigen::Map<Eigen::MatrixXd> randn_theta,
    Eigen::Map<Eigen::MatrixXd> randn_e,
    int grainSize) {

  struct sampleColumn : public Worker {
    Eigen::MatrixXd Y, W, prior_mean, prior_prec, randn_theta, randn_e;
    const std::vector<Eigen_SpMat> chol_R_list;
    RVector<int> h2s_index;
    Eigen::VectorXd tot_Eta_prec;
    Eigen::MatrixXd &coefs;

    sampleColumn(Eigen::MatrixXd Y,
                 Eigen::MatrixXd W,
                 Eigen::MatrixXd prior_mean,
                 Eigen::MatrixXd prior_prec,
                 const std::vector<Eigen_SpMat> chol_R_list,
                 const Rcpp::IntegerVector h2s_index,
                 Eigen::VectorXd tot_Eta_prec,
                 Eigen::MatrixXd randn_theta,
                 Eigen::MatrixXd randn_e,
                 Eigen::MatrixXd &coefs):
      Y(Y), W(W), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
      chol_R_list(chol_R_list), h2s_index(h2s_index), tot_Eta_prec(tot_Eta_prec), coefs(coefs) {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        int h2_index = h2s_index[j] - 1;
        Eigen_SpMat chol_R = chol_R_list[h2_index];
        coefs.col(j) = sample_MME_single_diagK(Y.col(j), W, prior_mean.col(j), prior_prec.col(j), chol_R, tot_Eta_prec[j], randn_theta.col(j),randn_e.col(j));
      }
    }
  };

  int b = randn_theta.rows();
  int p = randn_theta.cols();

  std::vector<Eigen_SpMat> chol_R_list;
  for(int i = 0; i < max(h2s_index); i++){
    chol_R_list.push_back(Rcpp::as<MSpMat>(chol_Rs[i]));
  }

  Eigen::MatrixXd coefs(b,p);

  sampleColumn sampler(Y,W,prior_mean,prior_prec,chol_R_list,h2s_index,tot_Eta_prec,randn_theta,randn_e, coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return(coefs);
}



// -------------------------------------------- //
// ------------ sample_MME_ZKZts -------------- //
// -------------------------------------------- //

Eigen::VectorXd sample_MME_single_diagR(
    Eigen::VectorXd Y,
    Eigen_SpMat W,
    Eigen_SpMat chol_C,
    double pe,
    Eigen_SpMat chol_K_inv,
    double tot_Eta_prec,
    Eigen::VectorXd randn_theta,
    Eigen::VectorXd randn_e
){
  Eigen::VectorXd theta_star = chol_K_inv.triangularView<Eigen::Upper>().solve(randn_theta);
  Eigen::VectorXd e_star = randn_e / sqrt(pe);
  Eigen::MatrixXd W_theta_star = W * theta_star;

  Eigen::VectorXd Y_resid = Y - W_theta_star - e_star;
  Eigen::VectorXd WtRiy = W.transpose() * (Y_resid * pe);

  Eigen::VectorXd theta_tilda = chol_C.transpose().triangularView<Eigen::Upper>().solve(chol_C.triangularView<Eigen::Lower>().solve(WtRiy));

  Eigen::VectorXd theta = theta_tilda / tot_Eta_prec + theta_star;

  return theta;
}

// [[Rcpp::export]]
Eigen::VectorXd sample_MME_single_diagR_c(
    Eigen::Map<Eigen::VectorXd> Y,
    MSpMat W,
    MSpMat chol_C,
    double pe,
    MSpMat chol_K_inv,
    double tot_Eta_prec,
    Eigen::Map<Eigen::VectorXd> randn_theta,
    Eigen::Map<Eigen::VectorXd> randn_e
){
  return sample_MME_single_diagR(Y,W,chol_C,pe,chol_K_inv,tot_Eta_prec,randn_theta,randn_e);
}

// [[Rcpp::export()]]
Eigen::MatrixXd sample_MME_ZKZts_c(
    Eigen::Map<Eigen::MatrixXd> Y,
    MSpMat W,
    Eigen::Map<Eigen::VectorXd> tot_Eta_prec,
    Rcpp::List chol_Cs,
    Rcpp::List chol_K_invs,
    Eigen::Map<Eigen::MatrixXd> h2s,
    Rcpp::IntegerVector h2s_index,
    Eigen::Map<Eigen::MatrixXd> randn_theta,
    Eigen::Map<Eigen::MatrixXd> randn_e,
    int grainSize) {

  struct sampleColumn : public Worker {
    Eigen::MatrixXd Y;
    Eigen_SpMat W;
    const std::vector<Eigen_SpMat> chol_C_list,chol_K_inv_list;
    Eigen::VectorXd pes,tot_Eta_prec;
    RVector<int> h2s_index;
    Eigen::MatrixXd randn_theta, randn_e;
    Eigen::MatrixXd &coefs;

    sampleColumn(Eigen::MatrixXd Y,
                 Eigen_SpMat W,
                 const std::vector<Eigen_SpMat> chol_C_list,
                 const std::vector<Eigen_SpMat> chol_K_inv_list,
                 Eigen::VectorXd pes,
                 Eigen::VectorXd tot_Eta_prec,
                 Rcpp::IntegerVector h2s_index,
                 Eigen::MatrixXd randn_theta, Eigen::MatrixXd randn_e,
                 Eigen::MatrixXd &coefs):
      Y(Y), W(W),
      chol_C_list(chol_C_list), chol_K_inv_list(chol_K_inv_list),
      pes(pes), tot_Eta_prec(tot_Eta_prec),h2s_index(h2s_index),
      randn_theta(randn_theta), randn_e(randn_e),
      coefs(coefs) {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        int h2_index = h2s_index[j] - 1;
        Eigen_SpMat chol_C = chol_C_list[h2_index];
        Eigen_SpMat chol_K_inv = chol_K_inv_list[h2_index];
        chol_K_inv *= sqrt(tot_Eta_prec[j]);
        coefs.col(j) = sample_MME_single_diagR(Y.col(j), W, chol_C, pes[j],chol_K_inv, tot_Eta_prec[j], randn_theta.col(j),randn_e.col(j));
      }
    }
  };
  int b = randn_theta.rows();
  int p = randn_theta.cols();

  std::vector<Eigen_SpMat> chol_C_list,chol_K_inv_list;
  for(int i = 0; i < max(h2s_index); i++){
    chol_C_list.push_back(Rcpp::as<MSpMat>(chol_Cs[i]));
    chol_K_inv_list.push_back(Rcpp::as<MSpMat>(chol_K_invs[i]));
  }

  Eigen::MatrixXd coefs(b,p);
  Eigen::VectorXd h2_e = 1.0 - h2s.colwise().sum().array();
  Eigen::VectorXd pes = tot_Eta_prec.array() / h2_e.array();

  sampleColumn sampler(Y,W,chol_C_list,chol_K_inv_list,pes,tot_Eta_prec,h2s_index,randn_theta,randn_e,coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return(coefs);
}




// -------------------------------------------- //
// -------------- tot_prec_scores ------------- //
// -------------------------------------------- //
// [[Rcpp::export()]]
Eigen::VectorXd tot_prec_scores(
    Eigen::Map<Eigen::MatrixXd> Y,
    Rcpp::List chol_Sigmas,
    Rcpp::IntegerVector h2s_index,
    int grainSize)
{

  struct sampleColumn : public Worker {
    Eigen::MatrixXd Y;
    Eigen_SpMat W;
    const std::vector<Eigen_SpMat> chol_Sigma_list;
    RVector<int> h2s_index;
    Eigen::VectorXd &scores;

    sampleColumn(Eigen::MatrixXd Y,
                 const std::vector<Eigen_SpMat> chol_Sigma_list,
                 Rcpp::IntegerVector h2s_index,
                 Eigen::VectorXd &scores):
      Y(Y), chol_Sigma_list(chol_Sigma_list), h2s_index(h2s_index),
      scores(scores) {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        int h2_index = h2s_index[j] - 1;
        Eigen_SpMat chol_Sigma = chol_Sigma_list[h2_index];
        Eigen::VectorXd score = chol_Sigma.triangularView<Eigen::Lower>().solve(Y.col(j));
        scores[j] = score.dot(score);
      }
    }
  };

  int p = Y.cols();
  Eigen::VectorXd scores(p);

  std::vector<Eigen_SpMat> chol_Sigma_list;
  for(int i = 0; i < max(h2s_index); i++){
    chol_Sigma_list.push_back(Rcpp::as<MSpMat>(chol_Sigmas[i]));
  }

  sampleColumn sampler(Y,chol_Sigma_list,h2s_index,scores);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return scores;
}



// -------------------------------------------- //
// ---------------- sample h2s ---------------- //
// -------------------------------------------- //

// [[Rcpp::export()]]
Eigen::VectorXd log_p_h2s(
    Eigen::Map<Eigen::MatrixXd> Y,
    Eigen::Map<Eigen::VectorXd> tot_Eta_prec,
    Rcpp::List chol_Sigmas,
    Eigen::Map<Eigen::VectorXd> log_det_Sigmas,
    Eigen::Map<Eigen::VectorXd> discrete_priors,
    int grainSize)
{

  struct sampleColumn : public Worker {
    Eigen::MatrixXd Y;
    Eigen::VectorXd tot_Eta_prec;
    const std::vector<Eigen_SpMat> chol_Sigma_list;
    Eigen::VectorXd log_det_Sigmas;
    Eigen::VectorXd discrete_priors;
    Eigen::MatrixXd &log_ps;

    sampleColumn(Eigen::MatrixXd Y,
                 Eigen::VectorXd tot_Eta_prec,
                 const std::vector<Eigen_SpMat> chol_Sigma_list,
                 Eigen::VectorXd log_det_Sigmas,
                 Eigen::VectorXd discrete_priors,
                 Eigen::MatrixXd &log_ps):
      Y(Y), tot_Eta_prec(tot_Eta_prec), chol_Sigma_list(chol_Sigma_list),
      log_det_Sigmas(log_det_Sigmas), discrete_priors(discrete_priors), log_ps(log_ps) {}

    void operator()(std::size_t begin, std::size_t end) {
      // int b = discrete_priors.size();
      int p = Y.cols();
      int n = Y.rows();
      for(std::size_t i = begin; i < end; i++){
        Eigen_SpMat chol_Sigma = chol_Sigma_list[i];
        Eigen::VectorXd scores2(p);
        for(int j = 0; j < p; j++){
          Eigen::VectorXd x_std = chol_Sigma.triangularView<Eigen::Lower>().solve(Y.col(j));
          scores2[j] = tot_Eta_prec[j] * x_std.dot(x_std);
        }
        log_ps.row(i) = (-n/2.0 * log(2*M_PI) - 0.5 * (log_det_Sigmas[i] - n*tot_Eta_prec.array().log()) -
          0.5 * scores2.array() + log(discrete_priors[i]));
      }
    }
  };

  int b = discrete_priors.size();
  int p = Y.cols();
  // int n = Y.rows();

  std::vector<Eigen_SpMat> chol_Sigma_list;
  for(int i = 0; i < log_det_Sigmas.size(); i++){
    chol_Sigma_list.push_back(Rcpp::as<MSpMat>(chol_Sigmas[i]));
  }

  Eigen::MatrixXd log_ps(b,p);

  sampleColumn sampler(Y,tot_Eta_prec,chol_Sigma_list,log_det_Sigmas,discrete_priors,log_ps);
  RcppParallel::parallelFor(0,b,sampler,grainSize);
  return(log_ps);
}

// [[Rcpp::export()]]
Rcpp::IntegerVector sample_h2s(
    Eigen::Map<Eigen::ArrayXXd> log_ps,
    Eigen::Map<Eigen::VectorXd> rs,
    int grainSize
)
{

  struct sampleColumn : public Worker {
    Eigen::ArrayXXd log_ps;
    Eigen::VectorXd rs;
    RVector<int> h2s_index;

    sampleColumn(Eigen::ArrayXXd log_ps,
                 Eigen::VectorXd rs,
                 Rcpp::IntegerVector h2s_index):
      log_ps(log_ps), rs(rs), h2s_index(h2s_index) {}

    void operator()(std::size_t begin, std::size_t end) {
      int b = log_ps.rows();
      for(std::size_t j = begin; j < end; j++){
        double max_col = log_ps.col(j).maxCoeff();
        double norm_factor = max_col + log((log_ps.col(j) - max_col).exp().sum());
        Eigen::VectorXd ps_j = (log_ps.col(j) - norm_factor).exp();
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
  // int b = log_ps.rows();
  Rcpp::IntegerVector h2s_index(p);

  sampleColumn sampler(log_ps,rs,h2s_index);
  RcppParallel::parallelFor(0,p,sampler,grainSize);

  return h2s_index;
}


// -------------------------------------------- //
// --------------- sample h2s MH -------------- //
// -------------------------------------------- //

double log_prob_h2_c(
    Eigen::VectorXd y,
    Eigen_SpMat chol_Sigma,
    double log_det_Sigma,
    int n,
    double tot_Eta_prec,
    double discrete_prior
){
  Eigen::VectorXd x_std = chol_Sigma.triangularView<Eigen::Lower>().solve(y);
  double score2 = tot_Eta_prec * x_std.dot(x_std);

  double log_p = -n/2.0 * log(2*M_PI) - 0.5*(log_det_Sigma - n*log(tot_Eta_prec)) - 0.5 * score2 + log(discrete_prior);
  return log_p;
}

Eigen::VectorXd find_candidate_states(
    Eigen::ArrayXXd h2s_matrix,
    double step_size,
    int old_state
) {
  Eigen::VectorXd dists = (h2s_matrix.colwise() - h2s_matrix.col(old_state)).abs().matrix().colwise().sum();
  Eigen::VectorXd indices(dists.size());
  int count = 0;
  for(int i = 0; i < dists.size(); i++){
    if(dists[i] < step_size & dists[i] > 0) {
      indices[count] = i;
      count++;
    }
  }
  return indices.head(count);
}

int randomInt(int max){
  Rcpp::NumericVector x = Rcpp::runif(1)*max;
  int i = x[0];
  return i;
}

// [[Rcpp::export()]]
Rcpp::IntegerVector sample_h2s_discrete_MH_c(
    Eigen::Map<Eigen::MatrixXd> Y,
    Eigen::Map<Eigen::MatrixXd> h2s_matrix,
    Rcpp::List chol_Sigmas,
    Eigen::Map<Eigen::VectorXd> log_det_Sigmas,
    Eigen::Map<Eigen::VectorXd> tot_Eta_prec,
    Eigen::Map<Eigen::VectorXd> discrete_priors,
    Eigen::Map<Eigen::VectorXd> r_draws,
    double step_size,
    Rcpp::IntegerVector h2_index,    // Rcpp::NumericVector state_draws
    Rcpp::NumericVector state_draws,
    int grainSize
){

  struct sampleColumn : public Worker {
    const Eigen::MatrixXd Y;
    const Eigen::MatrixXd h2s_matrix;
    const std::vector<Eigen_SpMat> chol_Sigma_list;
    const Eigen::VectorXd log_det_Sigmas;
    const Eigen::VectorXd tot_Eta_prec;
    const Eigen::VectorXd discrete_priors;
    const Eigen::VectorXd r_draws;
    const double step_size;
    const RVector<int> h2_index;
    RVector<int> new_index;

    sampleColumn(const Eigen::MatrixXd Y,
                 const Eigen::MatrixXd h2s_matrix,
                 const std::vector<Eigen_SpMat> chol_Sigma_list,
                 const Eigen::VectorXd log_det_Sigmas,
                 const Eigen::VectorXd tot_Eta_prec,
                 const Eigen::VectorXd discrete_priors,
                 const Eigen::VectorXd r_draws,
                 const double step_size,
                 const Rcpp::IntegerVector h2_index,
                 Rcpp::IntegerVector new_index):

      Y(Y), h2s_matrix(h2s_matrix),
      chol_Sigma_list(chol_Sigma_list), log_det_Sigmas(log_det_Sigmas),
      tot_Eta_prec(tot_Eta_prec),  discrete_priors(discrete_priors),
      r_draws(r_draws),step_size(step_size),h2_index(h2_index),new_index(new_index) {}

    void operator()(std::size_t begin, std::size_t end) {
      int n = Y.rows();
      for(std::size_t j = begin; j < end; j++){
        int old_state = h2_index[j] - 1;
        Eigen_SpMat chol_Sigma_old = chol_Sigma_list[old_state];
        double old_log_p = log_prob_h2_c(Y.col(j),chol_Sigma_old,log_det_Sigmas[old_state],n,tot_Eta_prec[j],discrete_priors[old_state]);

        Eigen::VectorXd candidate_new_states = find_candidate_states(h2s_matrix,step_size,old_state);
        int proposed_state = candidate_new_states[randomInt(candidate_new_states.size())];

        Eigen_SpMat chol_Sigma_new = chol_Sigma_list[proposed_state];
        double new_log_p = log_prob_h2_c(Y.col(j),chol_Sigma_new,log_det_Sigmas[proposed_state],n,tot_Eta_prec[j],discrete_priors[proposed_state]);

        Eigen::VectorXd candidate_states_from_new_state = find_candidate_states(h2s_matrix,step_size,proposed_state);

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

  // int n = Y.rows();
  int p = Y.cols();

  std::vector<Eigen_SpMat> chol_Sigma_list;
  for(int i = 0; i < log_det_Sigmas.size(); i++){
    chol_Sigma_list.push_back(Rcpp::as<MSpMat>(chol_Sigmas[i]));
  }

  Rcpp::IntegerVector new_index(p);

  sampleColumn sampler(Y,h2s_matrix,chol_Sigma_list,log_det_Sigmas,tot_Eta_prec,discrete_priors,r_draws,step_size,h2_index,new_index);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return new_index;
}

