#include <RcppEigen.h>
#include <RcppParallel.h>
#include <math.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace RcppParallel;
using Eigen::Map;               	// 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::ArrayXXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
// using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers
using Eigen::Upper;
using Eigen::Lower;
typedef Eigen::MappedSparseMatrix<double> MSpMat;
typedef Eigen::SparseMatrix<double> SpMat;


// //[[Rcpp::export()]]
double log_prob_h2_c(
                      VectorXd y,
                      SpMat chol_Sigma,
                      double log_det_Sigma,
                      int n,
                      double tot_Eta_prec,
                      double discrete_prior
                      ){
  VectorXd x_std = chol_Sigma.triangularView<Lower>().solve(y);
  double score2 = tot_Eta_prec * x_std.dot(x_std);

  double log_p = -n/2.0 * log(2*M_PI) - 0.5*(log_det_Sigma - n*log(tot_Eta_prec)) - 0.5 * score2 + log(discrete_prior);
  return log_p;
}

// // [[Rcpp::export()]]
VectorXd find_candidate_states(
                                      ArrayXXd h2s_matrix,
                                      double step_size,
                                      int old_state
) {
  VectorXd dists = (h2s_matrix.colwise() - h2s_matrix.col(old_state)).abs().matrix().colwise().sum();
  VectorXd indices(dists.size());
  int count = 0;
  for(int i = 0; i < dists.size(); i++){
    if(dists[i] < step_size & dists[i] > 0) {
      indices[count] = i;
      count++;
    }
  }
  return indices.head(count);
}

// [[Rcpp::export()]]
int randomInt(int max){
  Rcpp::NumericVector x = Rcpp::runif(1)*max;
  int i = x[0];
  return i;
}


// [[Rcpp::export()]]
Rcpp::IntegerVector sample_h2s_discrete_MH_c(
    Map<MatrixXd> Y,
    Map<MatrixXd> h2s_matrix,
    Rcpp::List chol_Sigmas,
    Map<VectorXd> log_det_Sigmas,
    Map<VectorXd> tot_Eta_prec,
    Map<VectorXd> discrete_priors,
    Map<VectorXd> r_draws,
    double step_size,
    Rcpp::IntegerVector h2_index,
    Rcpp::NumericVector state_draws
){
  int n = Y.rows();
  int p = Y.cols();

  std::vector<SpMat> chol_Sigma_list;
  for(int i = 0; i < log_det_Sigmas.size(); i++){
    chol_Sigma_list.push_back(Rcpp::as<MSpMat>(chol_Sigmas[i]));
  }

  Rcpp::IntegerVector new_index(p);
  for(int j = 0; j < p; j++){
    int old_state = h2_index[j] - 1;
    SpMat chol_Sigma_old = chol_Sigma_list[old_state];
    double old_log_p = log_prob_h2_c(Y.col(j),chol_Sigma_old,log_det_Sigmas[old_state],n,tot_Eta_prec[j],discrete_priors[old_state]);

    VectorXd candidate_new_states = find_candidate_states(h2s_matrix,step_size,old_state);
    // int r = state_draws[j] * (candidate_new_states.size());
    // int proposed_state = candidate_new_states[r];
    int proposed_state = candidate_new_states[randomInt(candidate_new_states.size())];

    SpMat chol_Sigma_new = chol_Sigma_list[proposed_state];
    double new_log_p = log_prob_h2_c(Y.col(j),chol_Sigma_new,log_det_Sigmas[proposed_state],n,tot_Eta_prec[j],discrete_priors[proposed_state]);

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
  return new_index;
}

// [[Rcpp::export()]]
Rcpp::IntegerVector sample_h2s_discrete_MH_c2(
    Map<MatrixXd> Y,
    Map<MatrixXd> h2s_matrix,
    Rcpp::List chol_Sigmas,
    Map<VectorXd> log_det_Sigmas,
    Map<VectorXd> tot_Eta_prec,
    Map<VectorXd> discrete_priors,
    Map<VectorXd> r_draws,
    double step_size,
    Rcpp::IntegerVector h2_index,    // Rcpp::NumericVector state_draws
    Rcpp::NumericVector state_draws,
    int grainSize
){

  struct sampleColumn : public Worker {
    const MatrixXd Y;
    const MatrixXd h2s_matrix;
    const std::vector<SpMat> chol_Sigma_list;
    const VectorXd log_det_Sigmas;
    const VectorXd tot_Eta_prec;
    const VectorXd discrete_priors;
    const VectorXd r_draws;
    const double step_size;
    const RVector<int> h2_index;
    const RVector<double> state_draws;
    RVector<int> new_index;

    sampleColumn(const MatrixXd Y,
                 const MatrixXd h2s_matrix,
                 const std::vector<SpMat> chol_Sigma_list,
                 const VectorXd log_det_Sigmas,
                 const VectorXd tot_Eta_prec,
                 const VectorXd discrete_priors,
                 const VectorXd r_draws,
                 const double step_size,
                 const Rcpp::IntegerVector h2_index,
                 const Rcpp::NumericVector state_draws,
                 Rcpp::IntegerVector new_index):

      Y(Y), h2s_matrix(h2s_matrix),
      chol_Sigma_list(chol_Sigma_list), log_det_Sigmas(log_det_Sigmas),
      tot_Eta_prec(tot_Eta_prec),  discrete_priors(discrete_priors),
      r_draws(r_draws),step_size(step_size),h2_index(h2_index),state_draws(state_draws),new_index(new_index) {}

    void operator()(std::size_t begin, std::size_t end) {
      int n = Y.rows();
      for(std::size_t j = begin; j < end; j++){
        int old_state = h2_index[j] - 1;
        SpMat chol_Sigma_old = chol_Sigma_list[old_state];
        double old_log_p = log_prob_h2_c(Y.col(j),chol_Sigma_old,log_det_Sigmas[old_state],n,tot_Eta_prec[j],discrete_priors[old_state]);

        VectorXd candidate_new_states = find_candidate_states(h2s_matrix,step_size,old_state);
        // int r = state_draws[j] * (candidate_new_states.size());
        // int proposed_state = candidate_new_states[r];
        int proposed_state = candidate_new_states[randomInt(candidate_new_states.size())];

        SpMat chol_Sigma_new = chol_Sigma_list[proposed_state];
        double new_log_p = log_prob_h2_c(Y.col(j),chol_Sigma_new,log_det_Sigmas[proposed_state],n,tot_Eta_prec[j],discrete_priors[proposed_state]);

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

  int n = Y.rows();
  int p = Y.cols();

  std::vector<SpMat> chol_Sigma_list;
  for(int i = 0; i < log_det_Sigmas.size(); i++){
    chol_Sigma_list.push_back(Rcpp::as<MSpMat>(chol_Sigmas[i]));
  }

  Rcpp::IntegerVector new_index(p);

  sampleColumn sampler(Y,h2s_matrix,chol_Sigma_list,log_det_Sigmas,tot_Eta_prec,discrete_priors,r_draws,step_size,h2_index,state_draws,new_index);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return new_index;
}

// [[Rcpp::export()]]
VectorXd log_p_h2s(
    Map<MatrixXd> Y,
    Map<VectorXd> tot_Eta_prec,
    Rcpp::List chol_Sigmas,
    Map<VectorXd> log_det_Sigmas,
    Map<VectorXd> discrete_priors,
    int grainSize)
  {

  struct sampleColumn : public Worker {
    MatrixXd Y;
    VectorXd tot_Eta_prec;
    const std::vector<SpMat> chol_Sigma_list;
    VectorXd log_det_Sigmas;
    VectorXd discrete_priors;
    MatrixXd &log_ps;

    sampleColumn(MatrixXd Y,
                 VectorXd tot_Eta_prec,
                 const std::vector<SpMat> chol_Sigma_list,
                 VectorXd log_det_Sigmas,
                 VectorXd discrete_priors,
                 MatrixXd &log_ps):
      Y(Y), tot_Eta_prec(tot_Eta_prec), chol_Sigma_list(chol_Sigma_list),
      log_det_Sigmas(log_det_Sigmas), discrete_priors(discrete_priors), log_ps(log_ps) {}

    void operator()(std::size_t begin, std::size_t end) {
      // int b = discrete_priors.size();
      int p = Y.cols();
      int n = Y.rows();
      for(std::size_t i = begin; i < end; i++){
        SpMat chol_Sigma = chol_Sigma_list[i];
        VectorXd scores2(p);
        for(int j = 0; j < p; j++){
          VectorXd x_std = chol_Sigma.triangularView<Lower>().solve(Y.col(j));
          scores2[j] = tot_Eta_prec[j] * x_std.dot(x_std);
        }
        log_ps.row(i) = (-n/2.0 * log(2*M_PI) - 0.5 * (log_det_Sigmas[i] - n*tot_Eta_prec.array().log()) -
          0.5 * scores2.array() + log(discrete_priors[i]));
      }
    }
  };

  int b = discrete_priors.size();
  int p = Y.cols();
  int n = Y.rows();

  std::vector<SpMat> chol_Sigma_list;
  for(int i = 0; i < log_det_Sigmas.size(); i++){
    chol_Sigma_list.push_back(Rcpp::as<MSpMat>(chol_Sigmas[i]));
  }

  MatrixXd log_ps(b,p);

  sampleColumn sampler(Y,tot_Eta_prec,chol_Sigma_list,log_det_Sigmas,discrete_priors,log_ps);
  RcppParallel::parallelFor(0,b,sampler,grainSize);
  // for(int i = 0; i < b; i++){
  //   SpMat chol_Sigma(Rcpp::as<MSpMat>(chol_Sigmas[i]));
  //   VectorXd scores2(p);
  //   for(int j = 0; j < p; j++){
  //     VectorXd x_std = chol_Sigma.triangularView<Lower>().solve(Y.col(j));
  //     scores2[j] = tot_Eta_prec[j] * x_std.dot(x_std);
  //   }
  //   log_ps.row(i) = (-n/2.0 * log(2*M_PI) - 0.5 * (log_det_Sigmas[i] - n*tot_Eta_prec.array().log()) -
  //                     0.5 * scores2.array() + log(discrete_priors[i]));
  // }
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
  int b = log_ps.rows();
  Rcpp::IntegerVector h2s_index(p);

  sampleColumn sampler(log_ps,rs,h2s_index);
  RcppParallel::parallelFor(0,p,sampler,grainSize);

  // for(int j = 0; j < p; j++){
  //   double max_col = log_ps.col(j).maxCoeff();
  //   double norm_factor = max_col + log((log_ps.col(j) - max_col).exp().sum());
  //   VectorXd ps_j = (log_ps.col(j) - norm_factor).exp();
  //   h2s_index[j] = 1;
  //   double cumsum = 0;
  //   for(int i = 0; i < b; i++){
  //     cumsum += ps_j[i];
  //     if(rs[j] > cumsum) {
  //       h2s_index[j] ++;
  //     }
  //   }
  // }
  return h2s_index;
}


// [[Rcpp::export()]]
VectorXd tot_prec_scores(
    Map<MatrixXd> Y,
    Rcpp::List chol_Sigmas,
    Rcpp::IntegerVector h2s_index,
    int grainSize)
  {

  struct sampleColumn : public Worker {
    MatrixXd Y;
    SpMat W;
    const std::vector<SpMat> chol_Sigma_list;
    RVector<int> h2s_index;
    VectorXd &scores;

    sampleColumn(MatrixXd Y,
                 const std::vector<SpMat> chol_Sigma_list,
                 Rcpp::IntegerVector h2s_index,
                 VectorXd &scores):
      Y(Y), chol_Sigma_list(chol_Sigma_list), h2s_index(h2s_index),
      scores(scores) {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        int h2_index = h2s_index[j] - 1;
        SpMat chol_Sigma = chol_Sigma_list[h2_index];
        VectorXd score = chol_Sigma.triangularView<Lower>().solve(Y.col(j));
        scores[j] = score.dot(score);
      }
    }
  };

  int p = Y.cols();
  VectorXd scores(p);

  std::vector<SpMat> chol_Sigma_list;
  for(int i = 0; i < max(h2s_index); i++){
    chol_Sigma_list.push_back(Rcpp::as<MSpMat>(chol_Sigmas[i]));
  }

  sampleColumn sampler(Y,chol_Sigma_list,h2s_index,scores);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return scores;
  }
