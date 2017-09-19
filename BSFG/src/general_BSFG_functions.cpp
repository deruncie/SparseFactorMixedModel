#include <math.h>
#include <iostream>
#include "BSFG_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace RcppParallel;


// -------------------------------------------- //
// ---------- helper functions --------- //
// -------------------------------------------- //
// functions to speed up sparse multiplication and conversion to dense matrices
// [[Rcpp::export()]]
MatrixXd SxD(MSpMat X, Map<MatrixXd> Y){
  return(X * Y);
}
// [[Rcpp::export()]]
MatrixXd SxS(MSpMat X, MSpMat Y){
  return(X * Y);
}

// [[Rcpp::export()]]
MatrixXd rstdnorm_mat(int n,int p) {  // returns nxp matrix
  VectorXd X_vec = as<VectorXd>(rnorm(n*p));
  Map<MatrixXd> X_mat(X_vec.data(),n,p);
  return(X_mat);
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

// -------------------------------------------- //
// ---------- sample_MME_fixedEffects --------- //
// -------------------------------------------- //

// [[Rcpp::export()]]
VectorXd sample_MME_single_diagK(  // returns b x 1 vector
    VectorXd y,           // nx1
    MatrixXd X,           // nxb
    VectorXd prior_mean,  // bx1
    VectorXd prior_prec,  // bx1
    SpMat chol_R,        // nxn upper triangular Cholesky decomposition of R. Format: dgCMatrix
    VectorXd randn_theta, // bx1
    VectorXd randn_e      // 0x1 or nx1. 0x1 if b<n
){
  if(randn_e.size() == 0){
    MatrixXd RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X);
    VectorXd XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y);
    VectorXd XtRinvy_std_mu = XtRinvy + prior_prec.asDiagonal()*prior_mean;
    MatrixXd C = RinvSqX.transpose() * RinvSqX;
    C.diagonal() += prior_prec;
    LLT<MatrixXd> C_llt;
    C_llt.compute(C);
    MatrixXd chol_C = C_llt.matrixU();

    VectorXd b = chol_C.transpose().triangularView<Lower>().solve(XtRinvy_std_mu);
    b += randn_theta;
    b = chol_C.triangularView<Upper>().solve(b);
    return(b);
  } else {
    // should check that randn_e.size() == n
    VectorXd theta_star = randn_theta.array() / prior_prec.cwiseSqrt().array();
    theta_star += prior_mean;
    VectorXd e_star = chol_R * randn_e;
    MatrixXd X_theta_star = X * theta_star;
    VectorXd y_resid = y - X_theta_star - e_star;

    MatrixXd RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X);
    VectorXd XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y_resid);

    VectorXd theta_tilda;
    MatrixXd VAi = X * prior_prec.cwiseInverse().asDiagonal();
    MatrixXd inner = VAi*X.transpose() + chol_R.transpose() * chol_R;
    VectorXd VAiXtURinvy = VAi * XtRinvy;
    VectorXd outerXtURinvy = VAi.transpose() * inner.ldlt().solve(VAiXtURinvy);
    theta_tilda = XtRinvy.array() / prior_prec.array();
    theta_tilda -= outerXtURinvy;
    VectorXd theta = theta_star + theta_tilda;

    return(theta);
  }
}


// RcppParallel struct for sampling coefficients from a set of parallel regression problems:
// Y = X %*% B + E, where
// b[i,j] ~ N(prior_mean[i,j],1/prior_prec[i,j])
// E[,j] ~ N(0,t(chol_R) %*% chol(R) / tot_Eta_prec[j])
// where chol_R is selected from a list based on h2s_index[j]
// Y is complete, so everythign has the same dimensions
struct sample_MME_single_diagK_worker : public RcppParallel::Worker {
  MatrixXd Y;
  MatrixXd X;
  MatrixXd prior_mean, prior_prec, randn_theta, randn_e;
  const std::vector<MSpMat> chol_R_list;
  VectorXi h2s_index;
  VectorXd tot_Eta_prec;
  MatrixXd &coefs;

  sample_MME_single_diagK_worker(
    MatrixXd Y,           // nxp
    MatrixXd X,           // nxb
    MatrixXd prior_mean,  // bxp
    MatrixXd prior_prec,  // bxp
    const std::vector<MSpMat> chol_R_list, // each element a MSpMat nxn upper-triangular cholesky decomposition of R
    VectorXi h2s_index,   // px1, 1-based index
    VectorXd tot_Eta_prec,// px1
    MatrixXd randn_theta, // bxp
    MatrixXd randn_e,     // 0xp or nxp. The former if b<n, the latter if b >= n
    MatrixXd &coefs       // bxp
  ):
    Y(Y), X(X), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
    chol_R_list(chol_R_list), h2s_index(h2s_index), tot_Eta_prec(tot_Eta_prec), coefs(coefs) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++){
      int h2_index = h2s_index[j] - 1;
      SpMat chol_R = chol_R_list[h2_index];
      chol_R *= 1/sqrt(tot_Eta_prec[j]);
      coefs.col(j) = sample_MME_single_diagK(Y.col(j), X, prior_mean.col(j), prior_prec.col(j), chol_R, randn_theta.col(j),randn_e.col(j));
    }
  }
};


// -------------------------------------------------- //
// -- Versions of the independent residuals regression model --- //
// -------------------------------------------------- //


// basic version - all regression have same dimensions
// [[Rcpp::export()]]
MatrixXd sample_MME_fixedEffects_c(  // returns bxp matrix
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
  int n = Y.rows();

  MatrixXd randn_theta = rstdnorm_mat(b,p);
  MatrixXd randn_e;
  if(b < n) {
    randn_e = rstdnorm_mat(0,p);
  } else{
    randn_e = rstdnorm_mat(n,p);
  }

  std::vector<MSpMat> chol_R_list;
  for(int i = 0; i < h2s_index.maxCoeff(); i++){
    Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);
    chol_R_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
  }

  MatrixXd coefs(b,p);

  sample_MME_single_diagK_worker sampler(Y,X,prior_mean,prior_prec,chol_R_list,h2s_index,tot_Eta_prec,randn_theta,randn_e, coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return(coefs);
}


// [[Rcpp::export()]]
MatrixXd sample_coefs_set_c(    // return pxn matrix
    Rcpp::List model_matrices,  // List. Each element contains: y (n_i x t), X (n_i x p), position (n_i x 1)
    Map<VectorXd> tot_Y_prec,   // tx1
    Map<MatrixXd> prior_mean,   // pxn
    Map<MatrixXd> prior_prec,   // pxn
    int grainSize){

  int n = model_matrices.size();

  std::vector<MatrixXd> y_list;
  std::vector<MatrixXd> X_list;
  std::vector<MatrixXd> randn_theta_list;
  for(int i = 0; i < n; i++){
    Rcpp::List model_matrix_i = Rcpp::as<Rcpp::List>(model_matrices[i]);
    y_list.push_back(Rcpp::as<MatrixXd>(model_matrix_i["y"]));
    X_list.push_back(Rcpp::as<MatrixXd>(model_matrix_i["X"]));

    int t = y_list[i].cols();
    int b = X_list[i].cols();

    MatrixXd randn_theta = rstdnorm_mat(b,t);
    randn_theta_list.push_back(randn_theta);
  }

  int n_traits = y_list[0].cols();
  int p = X_list[0].cols()*n_traits;

  MatrixXd coefs(p,n);

  struct sampleColumn : public RcppParallel::Worker {
    std::vector<MatrixXd> y_list;
    std::vector<MatrixXd> X_list;
    std::vector<MatrixXd> randn_theta_list;
    VectorXd tot_Y_prec;
    MatrixXd prior_mean,prior_prec;
    int n_traits;
    MatrixXd &coefs;

    sampleColumn(
      std::vector<MatrixXd> y_list,
      std::vector<MatrixXd> X_list,
      std::vector<MatrixXd> randn_theta_list,
      VectorXd tot_Y_prec,
      MatrixXd prior_mean,
      MatrixXd prior_prec,
      int n_traits,
      MatrixXd &coefs) :
      y_list(y_list), X_list(X_list), randn_theta_list(randn_theta_list),
      tot_Y_prec(tot_Y_prec),prior_mean(prior_mean), prior_prec(prior_prec),n_traits(n_traits),
      coefs(coefs)
    {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        MatrixXd Y = y_list[j];
        MatrixXd X = X_list[j];
        int b = X.cols();
        int n_obs = Y.rows();
        SpMat I = MatrixXd::Identity(n_obs,n_obs).sparseView();
        MatrixXd randn_e = MatrixXd::Zero(0,n_traits);
        for(int t = 0; t < n_traits; t++) {
          SpMat chol_R = I * 1/sqrt(tot_Y_prec[t]);
          coefs.block(t*b,j,b,1) = sample_MME_single_diagK(Y.col(t), X,
                                                           prior_mean.block(t*b,j,b,1), prior_prec.block(t*b,j,b,1),
                                                           chol_R,
                                                           randn_theta_list[j].col(t),randn_e.col(t));
        }
      }
    }
  };


  sampleColumn sampler(y_list, X_list,randn_theta_list,tot_Y_prec,prior_mean,prior_prec,n_traits,coefs);
  RcppParallel::parallelFor(0,n,sampler,grainSize);

  return(coefs);
}

// [[Rcpp::export()]]
MatrixXd get_fitted_set_c(  // returns n_tot x p matrix in same order as data
    Rcpp::List model_matrices,  // List. Each element contains: y (n_i x t), X (n_i x b), s (n_i x 1), position (n_i x 1)
    Map<MatrixXd> coefs,  // b x n matrix
    int grainSize){

  std::vector<MatrixXd> X_list;
  std::vector<VectorXd> position_list;
  int total_obs = 0;
  int n = model_matrices.size();
  for(int i = 0; i < n; i++){
    Rcpp::List model_matrix_i = Rcpp::as<Rcpp::List>(model_matrices[i]);
    X_list.push_back(Rcpp::as<MatrixXd>(model_matrix_i["X"]));
    position_list.push_back(Rcpp::as<VectorXd>(model_matrix_i["position"]));

    int n_obs = X_list[i].rows();
    total_obs += n_obs;
  }

  Rcpp::List model_matrix_1 = Rcpp::as<Rcpp::List>(model_matrices[0]);
  MatrixXd Y = Rcpp::as<MatrixXd>(model_matrix_1["y"]);
  int n_traits = Y.cols();

  MatrixXd Y_fitted(total_obs,n_traits);

  struct sampleColumn : public RcppParallel::Worker {
    std::vector<MatrixXd> X_list;
    std::vector<VectorXd> position_list;
    std::vector<VectorXd> s_list;
    int n_traits;
    MatrixXd coefs;
    MatrixXd &Y_fitted;

    sampleColumn(
      std::vector<MatrixXd> X_list,
      std::vector<VectorXd> position_list,
      int n_traits,
      MatrixXd coefs,
      MatrixXd &Y_fitted) :
      X_list(X_list), position_list(position_list),
      n_traits(n_traits),
      coefs(coefs),Y_fitted(Y_fitted)
    {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        int b = X_list[j].cols();
        Map<MatrixXd> Eta_i(coefs.col(j).data(),b,n_traits);
        MatrixXd Y_fitted_j = X_list[j] * Eta_i;
        for(int i = 0; i < position_list[j].size(); i++) {
          Y_fitted.row(position_list[j][i]-1) = Y_fitted_j.row(i);
        }
      }
    }
  };

  sampleColumn sampler(X_list,position_list,n_traits,coefs,Y_fitted);
  RcppParallel::parallelFor(0,n,sampler,grainSize);

  return(Y_fitted);
}

// -------------------------------------------- //
// ------------ sample_MME_ZKZts -------------- //
// -------------------------------------------- //

// Samples from model:
// y = Zu + e
// u ~ N(0,K); solve(K) = t(chol_K_inv) %*% chol_K_inv
// e[i] ~ N(0,1/tot_Eta_prec)
// C = ZtRinvZ + diag(Kinv)
// [[Rcpp::export()]]
VectorXd sample_MME_single_diagR(
    VectorXd y,           // nx1
    SpMat Zt,             // nxr dgCMatrix
    SpMat chol_C,         // rxr dgCMatrix upper triangular: chol(ZtRinvZ + diag(Kinv))
    double pe,            // double
    VectorXd randn_theta  // rx1
){
  VectorXd b = Zt * y * pe;
  b = chol_C.transpose().triangularView<Lower>().solve(b);
  b += randn_theta;
  b = chol_C.triangularView<Upper>().solve(b);
  return(b);
}


struct sample_MME_single_diagR_worker : public RcppParallel::Worker {
  MatrixXd Y;
  SpMat Zt;
  const std::vector<MSpMat> chol_C_list;
  ArrayXd pes;
  VectorXd tot_Eta_prec;
  VectorXi h2s_index;
  MatrixXd randn_theta;
  MatrixXd &coefs;

  sample_MME_single_diagR_worker(MatrixXd Y,                                // nxp
                                  SpMat Zt,                                  // nxr
                                  const std::vector<MSpMat> chol_C_list,     // std::vector of rxr SpMat upper-triangle
                                  ArrayXd pes,                               // px1
                                  VectorXd tot_Eta_prec,                     // px1
                                  VectorXi h2s_index,                        // px1, 1-based index
                                  MatrixXd randn_theta,                      // rxp
                                  MatrixXd &coefs                            // rxp
  ):
    Y(Y), Zt(Zt),
    chol_C_list(chol_C_list),
    pes(pes), tot_Eta_prec(tot_Eta_prec),h2s_index(h2s_index),
    randn_theta(randn_theta),
    coefs(coefs) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++){
      int h2_index = h2s_index[j] - 1;
      SpMat chol_C = chol_C_list[h2_index] * sqrt(tot_Eta_prec[j]);  // scale C by tot_Eta_prec[j]. C = Zt*Rinv*Z + Kinv, where Rinv and Kinv are scaled by h2.
      coefs.col(j) = sample_MME_single_diagR(Y.col(j), Zt, chol_C, pes[j],randn_theta.col(j));
    }
  }
};


// samples random effects from model:
// Y = ZU + E
// U[,j] ~ N(0,1/tot_Eta_prec[j] * h2[j] * K)
// E[,j] ~ N(0,1/tot_Eta_prec[j] * (1-h2[j]) * I_n)
// For complete data, ie no missing obs.
// [[Rcpp::export()]]
MatrixXd sample_MME_ZKZts_c(
    Map<MatrixXd> Y,                    // nxp
    MSpMat Z,                           // nxr
    Map<VectorXd> tot_Eta_prec,         // px1
    Rcpp::List randomEffect_C_Choleskys, // List. Each element contains chol_C and chol_K_inv (both rxr dgCMatrix upper-triangle)
    Map<MatrixXd> h2s,                  // n_RE x p
    VectorXi h2s_index,                 // px1
    int grainSize) {

  int p = Y.cols();
  int r = Z.cols();

  MatrixXd randn_theta = rstdnorm_mat(r,p);

  std::vector<MSpMat> chol_C_list;
  for(int i = 0; i < h2s_index.maxCoeff(); i++){
    Rcpp::List randomEffect_C_Cholesky_i = Rcpp::as<Rcpp::List>(randomEffect_C_Choleskys[i]);
    chol_C_list.push_back(Rcpp::as<MSpMat>(randomEffect_C_Cholesky_i["chol_C"]));
  }

  MatrixXd U(r,p);
  ArrayXd h2_e = 1.0 - h2s.colwise().sum().array();
  ArrayXd pes = tot_Eta_prec.array() / h2_e.array();

  sample_MME_single_diagR_worker sampler(Y,Z.transpose(),chol_C_list,pes,tot_Eta_prec,h2s_index,randn_theta,U);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return(U);
}


// -------------------------------------------- //
// -------------- tot_prec_scores ------------- //
// -------------------------------------------- //


struct tot_prec_scores_worker : public RcppParallel::Worker {
  MatrixXd Y;
  SpMat W;
  const std::vector<MSpMat> chol_R_list;
  VectorXi h2s_index;
  VectorXd &scores;

  tot_prec_scores_worker(MatrixXd Y,                                // nxp
                         const std::vector<MSpMat> chol_R_list,     // std::vector of nxn SpMat upper-triangule
                         VectorXi h2s_index,                        // px1, 1-based index
                         VectorXd &scores                           // px1
                           ):
    Y(Y), chol_R_list(chol_R_list), h2s_index(h2s_index),
    scores(scores) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++){
      int h2_index = h2s_index[j] - 1;
      SpMat chol_R = chol_R_list[h2_index];
      VectorXd score = chol_R.transpose().triangularView<Lower>().solve(Y.col(j));
      scores[j] = score.dot(score);
    }
  }
};

// calculates normal scores for:
// Y_resid %*% Sigma_inv %*% t(Y_resid)
// assumes complete data Y
// [[Rcpp::export()]]
VectorXd tot_prec_scores(
    Map<MatrixXd> Y,              // nxp
    Rcpp::List Sigma_Choleskys,   // List of nxn dgCMatrix upper-triangule
    VectorXi h2s_index,           // px1
    int grainSize)
{


  int p = Y.cols();
  VectorXd scores(p);

  std::vector<MSpMat> chol_R_list;
  for(int i = 0; i < h2s_index.maxCoeff(); i++){
    Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);
    chol_R_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
  }

  tot_prec_scores_worker sampler(Y,chol_R_list,h2s_index,scores);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return scores;
}

// -------------------------------------------- //
// ---------------- sample h2s ---------------- //
// -------------------------------------------- //

struct log_ps_worker : public RcppParallel::Worker {
  MatrixXd Y;
  VectorXd tot_Eta_prec;
  const std::vector<MSpMat> chol_R_list;
  VectorXd log_det_Sigmas;
  VectorXd discrete_priors;
  MatrixXd &log_ps;

  log_ps_worker(MatrixXd Y,                           // nxp
               VectorXd tot_Eta_prec,                 // px1
               const std::vector<MSpMat> chol_R_list, // std::vector of nxn SpMat, upper-triangular
               VectorXd log_det_Sigmas,               // n_h2 x 1
               VectorXd discrete_priors,              // n_h2 x 1
               MatrixXd &log_ps                       // n_h2 x p
    ):
    Y(Y), tot_Eta_prec(tot_Eta_prec), chol_R_list(chol_R_list),
    log_det_Sigmas(log_det_Sigmas), discrete_priors(discrete_priors), log_ps(log_ps) {}

  void operator()(std::size_t begin, std::size_t end) {
    int p = Y.cols();
    int n = Y.rows();
    for(std::size_t i = begin; i < end; i++){
      SpMat chol_R = chol_R_list[i];
      VectorXd scores2(p);
      for(int j = 0; j < p; j++){
        VectorXd x_std = chol_R.transpose().triangularView<Lower>().solve(Y.col(j));
        scores2[j] = tot_Eta_prec[j] * x_std.dot(x_std);
      }
      log_ps.row(i) = (-n/2.0 * log(2*M_PI) - 0.5 * (log_det_Sigmas[i] - n*tot_Eta_prec.array().log()) -
        0.5 * scores2.array() + log(discrete_priors[i]));
    }
  }
};

// [[Rcpp::export()]]
MatrixXd log_p_h2s(
    Map<MatrixXd> Y,              // nxp
    Map<VectorXd> tot_Eta_prec,   // px1
    Rcpp::List Sigma_Choleskys,   // List. Each element contains: chol_Sigma: nxn dgCMatrix, upper-triangular, log_det
    Map<VectorXd> discrete_priors,// n_h2 x 1
    int grainSize)
{
  int b = discrete_priors.size();
  int p = Y.cols();

  std::vector<MSpMat> chol_R_list;
  VectorXd log_det_Sigmas(b);
  for(int i = 0; i < b; i++){
    Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);
    chol_R_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
    log_det_Sigmas[i] = Rcpp::as<double>(Sigma_Choleskys_i["log_det"]);
  }

  MatrixXd log_ps(b,p);

  log_ps_worker sampler(Y,tot_Eta_prec,chol_R_list,log_det_Sigmas,discrete_priors,log_ps);
  RcppParallel::parallelFor(0,b,sampler,grainSize);
  return(log_ps);
}

// [[Rcpp::export()]]
VectorXi sample_h2s(
    Map<ArrayXXd> log_ps,
    int grainSize
)
{

  int p = log_ps.cols();
  VectorXd rs = as<VectorXd>(runif(p));

  VectorXi h2s_index(p);


  struct sampleColumn : public RcppParallel::Worker {
    ArrayXXd log_ps;
    VectorXd rs;
    VectorXi &h2s_index;

    sampleColumn(ArrayXXd log_ps,
                 VectorXd rs,
                 VectorXi &h2s_index):
      log_ps(log_ps), rs(rs), h2s_index(h2s_index) {}

    void operator()(std::size_t begin, std::size_t end) {
      int b = log_ps.rows();
      for(std::size_t j = begin; j < end; j++){
        // for(int j = 0; j < p; j++){
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
  sampleColumn sampler(log_ps,rs,h2s_index);
  RcppParallel::parallelFor(0,p,sampler,grainSize);

  return h2s_index; // 1-based index
}


// -------------------------------------------- //
// --------------- sample h2s MH -------------- //
// -------------------------------------------- //

double log_prob_h2_c(
    VectorXd y,           // nx1
    SpMat chol_R,         // nxn dgCMatrix upper-triangular
    double log_det_Sigma, // double
    int n,                // int
    double tot_Eta_prec,  // double
    double discrete_prior // double
){
  VectorXd x_std = chol_R.transpose().triangularView<Lower>().solve(y);
  double score2 = tot_Eta_prec * x_std.dot(x_std);

  double log_p = -n/2.0 * log(2*M_PI) - 0.5*(log_det_Sigma - n*log(tot_Eta_prec)) - 0.5 * score2 + log(discrete_prior);
  return log_p;
}

struct sample_h2s_discrete_MH_worker : public RcppParallel::Worker {
  const MatrixXd Y;
  const MatrixXd h2s_matrix;
  const std::vector<MSpMat> chol_R_list;
  const VectorXd log_det_Sigmas;
  const VectorXd tot_Eta_prec;
  const VectorXd discrete_priors;
  const VectorXd r_draws;
  const VectorXd state_draws;
  const VectorXi h2_index;
  const double step_size;
  VectorXi &new_index;

  sample_h2s_discrete_MH_worker(const MatrixXd Y,                           // nxp
                                const MatrixXd h2s_matrix,                  // n_RE x n_h2
                                const std::vector<MSpMat> chol_R_list,      // std::vector of nxn SpMat, upper-triangular
                                const VectorXd log_det_Sigmas,              // n_h2 x 1
                                const VectorXd tot_Eta_prec,                // p x 1
                                const VectorXd discrete_priors,             // n_h2 x 1
                                const VectorXd r_draws,                     // px1
                                const VectorXd state_draws,                 // px1
                                const VectorXi h2_index,                    // px1 - 1-based index
                                const double step_size,                     // double
                                VectorXi &new_index                         // px1 - 1-based index
                                  ):

    Y(Y), h2s_matrix(h2s_matrix),
    chol_R_list(chol_R_list), log_det_Sigmas(log_det_Sigmas),
    tot_Eta_prec(tot_Eta_prec),  discrete_priors(discrete_priors),
    r_draws(r_draws),state_draws(state_draws),h2_index(h2_index),step_size(step_size),new_index(new_index) {}

  void operator()(std::size_t begin, std::size_t end) {
    int n = Y.rows();
    for(std::size_t j = begin; j < end; j++){
      int old_state = h2_index[j] - 1;
      double old_log_p = log_prob_h2_c(Y.col(j),chol_R_list[old_state],log_det_Sigmas[old_state],n,tot_Eta_prec[j],discrete_priors[old_state]);

      VectorXd candidate_new_states = find_candidate_states(h2s_matrix,step_size,old_state);
      int r = state_draws[j] * (candidate_new_states.size());
      int proposed_state = candidate_new_states[r];

      double new_log_p = log_prob_h2_c(Y.col(j),chol_R_list[proposed_state],log_det_Sigmas[proposed_state],n,tot_Eta_prec[j],discrete_priors[proposed_state]);

      VectorXd candidate_states_from_new_state = find_candidate_states(h2s_matrix,step_size,proposed_state);

      double forward_prob = 1.0 / candidate_new_states.size();
      double back_prob = 1.0 / candidate_states_from_new_state.size();

      double log_MH_ratio = new_log_p - old_log_p + log(forward_prob) - log(back_prob);

      if(log(r_draws[j]) < log_MH_ratio) {
        new_index[j] = proposed_state;
      } else {
        new_index[j] = old_state;
      }
      new_index[j] += 1;  // convert to 1-based index
    }
  }
};

// [[Rcpp::export()]]
VectorXi sample_h2s_discrete_MH_c(
    Map<MatrixXd> Y,                // nxp
    Map<VectorXd> tot_Eta_prec,     // px1
    Map<VectorXd> discrete_priors,  // n_h2 x 1
    VectorXi h2s_index,              // px1
    Map<MatrixXd> h2s_matrix,       // n_RE x n_h2
    Rcpp::List Sigma_Choleskys,     // List of nxn dgCMatrices, upper-triangular
    double step_size,               // double
    int grainSize
){

  int p = Y.cols();
  int b = discrete_priors.size();
  VectorXd r_draws = as<VectorXd>(runif(p));
  VectorXd state_draws = as<VectorXd>(runif(p));

  VectorXi new_index(p);

  std::vector<MSpMat> chol_R_list;
  VectorXd log_det_Sigmas(b);
  for(int i = 0; i < b; i++){
    Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);

    chol_R_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
    log_det_Sigmas[i] = Rcpp::as<double>(Sigma_Choleskys_i["log_det"]);
  }


  sample_h2s_discrete_MH_worker sampler(Y,h2s_matrix,chol_R_list,log_det_Sigmas,tot_Eta_prec,discrete_priors,r_draws,state_draws,h2s_index,step_size,new_index);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return new_index;
}


// -------------------------------------------------- //
// -- Sample factor scores --- //
// -------------------------------------------------- //

// Sample factor scores given factor loadings (F_a), factor residual variances (F_e_prec) and
// phenotype residuals
// Y - ZU = F * Lambda + E
// F[,k] ~ N(XFBF[,k] + ZUF[,k],1/F_e_prec[,j])
// E[,j] ~ N(0,1/resid_Eta_prec[,j])
// Sampling is done separately for each block of rows with the same pattern of missing observations
// [[Rcpp::export()]]
MatrixXd sample_factors_scores_c( // returns nxk matrix
    Map<MatrixXd> Eta_tilde,      // nxp
    Map<MatrixXd> prior_mean,     // nxk
    Map<MatrixXd> Lambda,         // nxp
    Map<VectorXd> resid_Eta_prec, // px1
    Map<VectorXd> F_e_prec        // kx1
) {
  int n = Eta_tilde.rows();
  int k = Lambda.cols();
  MatrixXd randn_draws = rstdnorm_mat(n,k);

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
