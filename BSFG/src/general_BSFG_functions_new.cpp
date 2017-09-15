#include <math.h>
#include <iostream>
#include "BSFG_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace RcppParallel;


// -------------------------------------------- //
// ---------- helper functions --------- //
// -------------------------------------------- //

MatrixXd get_block_matrix(MatrixXd Y, VectorXi Y_obs, VectorXi Y_cols){
  int n_obs = Y_obs.size();
  int n_cols = Y_cols.size();
  MatrixXd Y_set(n_obs,n_cols);
  for(int i = 0; i < n_obs; i++){
    int row_i = Y_obs[i]-1;
    for(int j = 0; j < n_cols; j++){
      int col_j = Y_cols[j]-1;
      Y_set.coeffRef(i,j) = Y.coeffRef(row_i,col_j);
    }
  }
  return(Y_set);
}

// -------------------------------------------- //
// ---------- sample_MME_fixedEffects --------- //
// -------------------------------------------- //

// draws a sample of the vector b from the model:
// y = X %*% beta + e
// with beta[i] ~ N(prior_mean[i],1/prior_prec[i]), i=1:b
// with e ~ N(0,t(chol_R) %*% chol(R)), i=1:n
// Uses sampling method from MCMCglmm, which requires draws b+n draws from N(0,1), which are passed as randn_theta and randn_e
// If b >= n, inverts the C matrix using Binomial Inverse Theorem
// [[Rcpp::export()]]
VectorXd sample_MME_single_diagK(  // returns b x 1 vector
    VectorXd y,           // nx1
    MatrixXd X,           // nxb
    VectorXd prior_mean,  // bx1
    VectorXd prior_prec,  // bx1
    SpMat chol_R,        // nxn upper triangular Cholesky decomposition of R. Format: dgCMatrix
    VectorXd randn_theta, // bx1
    VectorXd randn_e      // nx1
){
  int n = y.size();
  int b = X.cols();

  VectorXd theta_star = randn_theta.array() / prior_prec.cwiseSqrt().array();
  theta_star += prior_mean;
  VectorXd e_star = chol_R * randn_e;
  MatrixXd X_theta_star = X * theta_star;
  VectorXd y_resid = y - X_theta_star - e_star;

  MatrixXd RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X);
  VectorXd XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y_resid);

  VectorXd theta_tilda;

  if(b < n) {
    MatrixXd C = RinvSqX.transpose() * RinvSqX;
    C.diagonal() += prior_prec;
    theta_tilda = C.llt().solve(XtRinvy);
  } else{
    MatrixXd VAi = X * prior_prec.cwiseInverse().asDiagonal();
    MatrixXd inner = VAi*X.transpose() + chol_R.transpose() * chol_R;
    VectorXd VAiXtURinvy = VAi * XtRinvy;
    VectorXd outerXtURinvy = VAi.transpose() * inner.ldlt().solve(VAiXtURinvy);
    theta_tilda = XtRinvy.array() / prior_prec.array();
    theta_tilda -= outerXtURinvy;
  }

  VectorXd theta = theta_star + theta_tilda;

  return theta;
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
              VectorXi h2s_index,   // px1
              VectorXd tot_Eta_prec,// px1
              MatrixXd randn_theta, // bxp
              MatrixXd randn_e,     // nxp
              MatrixXd &coefs       // bxp
            ):
    Y(Y), X(X), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
    chol_R_list(chol_R_list), h2s_index(h2s_index), tot_Eta_prec(tot_Eta_prec), coefs(coefs) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++){
      int h2_index = h2s_index[j] - 1;
      SpMat chol_R = chol_R_list[h2_index];
      chol_R *= 1/sqrt(tot_Eta_prec[j]);
      coefs.col(j) = sample_MME_single_diagK(Y.col(j), X, prior_mean.col(j), prior_prec.col(j), chol_R_list[h2_index], randn_theta.col(j),randn_e.col(j));
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

  int n = X.rows();
  int b = X.cols();
  int p = Y.cols();

  MatrixXd randn_theta = rstdnorm_mat(b,p);
  MatrixXd randn_e = rstdnorm_mat(n,p);

  std::vector<MSpMat> chol_R_list;
  for(int i = 0; i < h2s_index.maxCoeff(); i++){
    Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);
    chol_R_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
  }

  MatrixXd coefs(b,p);

  sample_MME_single_diagK_worker sampler(Y,X,prior_mean,prior_prec,chol_R_list,h2s_index,tot_Eta_prec,randn_theta,randn_e,coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return(coefs);
}


// missing data version - sets of columns of Y have different patterns of missing data
// each set can be sampled at once using sample_MME_single_diagK_worker
// [[Rcpp::export()]]
MatrixXd sample_MME_fixedEffects_missing_c(  // returns bxp matrix
    Map<MatrixXd> Y,                // nxp
    Map<MatrixXd> X,                // nxb
    Rcpp::List Sigma_Cholesky_sets, // List. Each element contains: List of Cholesky decompositions of residuals. nxn upper triangular, dgCMatrix.
    Rcpp::List Missing_data_map,    // List. Each element containts Y_obs, Y_cols. The length should correspond to the length of Sigma_Choleksy_sets
    VectorXi h2s_index,             // px1 index of Cholesky matrix for each column
    Map<VectorXd> tot_Eta_prec,     // px1
    Map<MatrixXd> prior_mean,       // bxp
    Map<MatrixXd> prior_prec,       // bxp
    int grainSize) {

  int n = X.rows();
  int b = X.cols();
  int p = Y.cols();

  MatrixXd coefs(b,p);

  // sample sets of columns with same pattern of missing data as a block
  for(int col_set = 0; col_set < Missing_data_map.length(); col_set++){
    Rcpp::List Missing_data_map_set = Rcpp::as<Rcpp::List>(Missing_data_map[col_set]);

    VectorXi Y_obs = as<Map<VectorXi>>(Missing_data_map_set["Y_obs"]);
    VectorXi Y_cols = as<Map<VectorXi>>(Missing_data_map_set["Y_cols"]);

    int n_obs = Y_obs.size();
    int n_cols = Y_cols.size();

    // pull out rows of X corresponding to observed individuals
    MatrixXd X_set(n_obs,b);
    for(int i = 0; i < n_obs; i++) {
      X_set.row(i) = X.row(Y_obs[i]-1);
    }

    // pull out elements of Y corresponding to observed individuals
    MatrixXd Y_set = get_block_matrix(Y,Y_obs,Y_cols);

    // pull out prior_mean, prior_prec, randn_theta, randn_e for col_set
    MatrixXd prior_mean_set(b,n_cols);
    MatrixXd prior_prec_set(b,n_cols);
    for(int j = 0; j < n_cols; j++){
      int col_j = Y_cols[j]-1;
      prior_mean_set.col(j) = prior_mean.col(col_j);
      prior_prec_set.col(j) = prior_prec.col(col_j);
    }

    // pull out h2s_index_set and tot_Eta_prec_set
    VectorXd h2s_index_set(n_cols);
    VectorXd tot_Eta_prec_set(n_cols);
    for(int j = 0; j < n_cols; j++){
      int col_j = Y_cols[j]-1;
      h2s_index_set[j] = h2s_index[col_j];
      tot_Eta_prec_set[j] = tot_Eta_prec[col_j];
    }

    // pull out necessary chol_R matrices for set
    Rcpp::List Sigma_Choleskys_set = Rcpp::as<Rcpp::List>(Sigma_Cholesky_sets[col_set]);
    std::vector<MSpMat> chol_R_list_set;
    for(int i = 0; i < h2s_index_set.maxCoeff(); i++){
      Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys_set[i]);
      chol_R_list_set.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
    }

    // draw standard normals
    MatrixXd randn_theta_set = rstdnorm_mat(b,n_cols);
    MatrixXd randn_e_set = rstdnorm_mat(n_obs,n_cols);

    MatrixXd coefs_set(b,n_cols);

    // sample coefs_set as a block
    sample_MME_single_diagK_worker sampler(Y_set,X_set,prior_mean_set,prior_prec_set,chol_R_list_set,h2s_index_set,tot_Eta_prec_set,randn_theta_set,randn_e_set,coefs_set);
    RcppParallel::parallelFor(0,n_cols,sampler,grainSize);

    // copy new coefs into full matrix
    for(int j = 0; j < n_cols; j++){
      coefs.col(Y_cols[j]-1) = coefs_set.col(j);
    }
  }
  return(coefs);
}


// [[Rcpp::export()]]
MatrixXd sample_coefs_set_c(    // return nxp matrix
    Rcpp::List model_matrices,  // List. Each element contains: y (n_i x t), X (n_i x b), s (n_i x 1), position (n_i x 1)
    Map<VectorXd> tot_Y_prec, // px1
    Map<MatrixXd> prior_mean,   // nxp
    Map<MatrixXd> prior_prec,   // nxp
    int grainSize){

  int n = model_matrices.size();

  std::vector<MatrixXd> y_list;
  std::vector<MatrixXd> X_list;
  std::vector<VectorXd> s_list;
  std::vector<MatrixXd> randn_theta_list;
  std::vector<MatrixXd> randn_e_list;
  int randn_theta_index = 0;
  int randn_e_index = 0;
  int total_obs = 0;
  for(int i = 0; i < n; i++){
    Rcpp::List model_matrix_i = Rcpp::as<Rcpp::List>(model_matrices[i]);
    y_list.push_back(Rcpp::as<MatrixXd>(model_matrix_i["y"]));
    X_list.push_back(Rcpp::as<MatrixXd>(model_matrix_i["X"]));
    s_list.push_back(Rcpp::as<VectorXd>(model_matrix_i["s"]));

    int p = y_list[i].cols();
    int n_obs = y_list[i].rows();
    int b = X_list[i].cols();

    total_obs += n_obs;

    MatrixXd randn_theta = rstdnorm_mat(b,p);
    randn_theta_list.push_back(randn_theta);
    randn_theta_index += b*p;
    MatrixXd randn_e = rstdnorm_mat(n_obs,p);
    randn_e_list.push_back(randn_e);
    randn_e_index += n_obs*p;
  }


  int n_traits = y_list[0].cols();
  int b = randn_theta_list[0].rows()*n_traits;

  MatrixXd coefs(b,n);

  struct sampleColumn : public RcppParallel::Worker {
    std::vector<MatrixXd> y_list;
    std::vector<MatrixXd> X_list;
    std::vector<MatrixXd> randn_theta_list;
    std::vector<MatrixXd> randn_e_list;
    std::vector<VectorXd> s_list;
    VectorXd tot_Y_prec;
    MatrixXd prior_mean,prior_prec;
    int n_traits;
    MatrixXd &coefs;

    sampleColumn(
      std::vector<MatrixXd> y_list,
      std::vector<MatrixXd> X_list,
      std::vector<MatrixXd> randn_theta_list,
      std::vector<MatrixXd> randn_e_list,
      std::vector<VectorXd> s_list,
      VectorXd tot_Y_prec,
      MatrixXd prior_mean,
      MatrixXd prior_prec,
      int n_traits,
      MatrixXd &coefs) :
      y_list(y_list), X_list(X_list), randn_theta_list(randn_theta_list), randn_e_list(randn_e_list),s_list(s_list),
      tot_Y_prec(tot_Y_prec),prior_mean(prior_mean), prior_prec(prior_prec),n_traits(n_traits),
      coefs(coefs)
    {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        MatrixXd Y = y_list[j];
        MatrixXd X = X_list[j];
        int b = randn_theta_list[j].rows();
        int n_obs = randn_e_list[j].rows();
        VectorXd s = s_list[j];
        SpMat I = MatrixXd::Identity(n_obs,n_obs).sparseView();
        for(int t = 0; t < n_traits; t++) {
          SpMat chol_R = I * 1/sqrt(tot_Y_prec[t]);
          coefs.block(t*b,j,b,1) = sample_MME_single_diagK(Y.col(t), X,
                                                           prior_mean.block(t*b,j,b,1), prior_prec.block(t*b,j,b,1),
                                                           chol_R,
                                                           randn_theta_list[j].col(t),randn_e_list[j].col(t));
        }
      }
    }
  };


  sampleColumn sampler(y_list, X_list,randn_theta_list,randn_e_list,s_list,tot_Y_prec,prior_mean,prior_prec,n_traits,coefs);
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
// [[Rcpp::export]]
VectorXd sample_MME_single_diagR(
    VectorXd y,           // nx1
    SpMat Z,              // nxr dgCMatrix
    SpMat chol_C,         // rxr dgCMatrix upper triangular: chol(ZtRinvZ + diag(Kinv))
    double pe,            // double
    SpMat chol_K_inv,     // rxr dgCMatrix upper triangular: chol(solve(K))
    VectorXd randn_theta, // rx1
    VectorXd randn_e      // nx1
){
  VectorXd theta_star = chol_K_inv.triangularView<Upper>().solve(randn_theta);
  VectorXd e_star = randn_e / sqrt(pe);
  MatrixXd Z_theta_star = Z * theta_star;

  VectorXd y_resid = y - Z_theta_star - e_star;
  VectorXd ZtRiy_resid = Z.transpose() * (y_resid * pe);

  VectorXd theta_tilda = chol_C.triangularView<Upper>().solve(chol_C.transpose().triangularView<Lower>().solve(ZtRiy_resid));

  VectorXd theta = theta_tilda + theta_star;

  return theta;
}

struct sample_MME_single_diagR_worker : public RcppParallel::Worker {
  MatrixXd Y;
  SpMat Z;
  const std::vector<MSpMat> chol_C_list,chol_K_inv_list;
  VectorXd pes,tot_Eta_prec;
  VectorXi h2s_index;
  MatrixXd randn_theta, randn_e;
  MatrixXd &coefs;

  sample_MME_single_diagR_worker(MatrixXd Y,              // nxp
               SpMat Z,                                   // nxr
               const std::vector<MSpMat> chol_C_list,     // std::vector of rxr SpMat upper-triangle
               const std::vector<MSpMat> chol_K_inv_list, // std::vector of rxr SpMat upper-triangle
               VectorXd pes,                              // px1
               VectorXd tot_Eta_prec,                     // px1
               VectorXi h2s_index,                        // px1
               MatrixXd randn_theta,                      // rxp
               MatrixXd randn_e,                          // nxp
               MatrixXd &coefs                            // rxp
  ):
    Y(Y), Z(Z),
    chol_C_list(chol_C_list), chol_K_inv_list(chol_K_inv_list),
    pes(pes), tot_Eta_prec(tot_Eta_prec),h2s_index(h2s_index),
    randn_theta(randn_theta), randn_e(randn_e),
    coefs(coefs) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++){
      int h2_index = h2s_index[j] - 1;
      SpMat chol_C = chol_C_list[h2_index] * sqrt(tot_Eta_prec[j]);  // scale C by tot_Eta_prec[j]. C = Zt*Rinv*Z + Kinv, where Rinv and Kinv are scaled by h2.
      SpMat chol_K_inv = chol_K_inv_list[h2_index];
      chol_K_inv *= sqrt(tot_Eta_prec[j]);
      coefs.col(j) = sample_MME_single_diagR(Y.col(j), Z, chol_C, pes[j],chol_K_inv, randn_theta.col(j),randn_e.col(j));
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
    Map<VectorXi> h2s_index,            // px1
    int grainSize) {

  int n = Y.rows();
  int p = Y.cols();
  int r = Z.cols();

  MatrixXd randn_theta = rstdnorm_mat(r,p);
  MatrixXd randn_e = rstdnorm_mat(n,p);

  std::vector<MSpMat> chol_C_list,chol_K_inv_list;
  for(int i = 0; i < h2s_index.maxCoeff(); i++){
    Rcpp::List randomEffect_C_Cholesky_i = Rcpp::as<Rcpp::List>(randomEffect_C_Choleskys[i]);
    chol_C_list.push_back(Rcpp::as<MSpMat>(randomEffect_C_Cholesky_i["chol_C"]));
    chol_K_inv_list.push_back(Rcpp::as<MSpMat>(randomEffect_C_Cholesky_i["chol_K_inv"]));
  }

  MatrixXd U(r,p);
  VectorXd h2_e = 1.0 - h2s.colwise().sum().array();
  VectorXd pes = tot_Eta_prec.array() / h2_e.array();

  sample_MME_single_diagR_worker sampler(Y,Z,chol_C_list,chol_K_inv_list,pes,tot_Eta_prec,h2s_index,randn_theta,randn_e,U);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return(U);
}


// samples random effects from model:
// Y = ZU + E
// U[,j] ~ N(0,1/a_prec[j] * K)
// E[,j] ~ N(0,1/e_prec[j] * I_n)
// Where matrix Q and vectors s1 and s2 diagonalize ZZt + K
// invert_aZZt_Kinv has parameters to diagonalize a*Z*Z' + b*K for fast inversion
// with missing observations; Samples sets of columns with same patterns of missingness
// [[Rcpp::export()]]
MatrixXd sample_MME_ZKZts_missing_c (  // returns rxp matrix
    Map<MatrixXd> Y,                    // nxp
    MSpMat Z,                           // nxr
    Map<VectorXd> tot_Eta_prec,         // px1
    Rcpp::List randomEffect_C_Cholesky_sets, // List. Each element contains: List of pairs chol_C and chol_K_inv (both rxr dgCMatrix upper-triangle)
    Rcpp::List Missing_data_map,        // List. Each element containts Y_obs, Y_cols. The length should correspond to the length of randomEffect_C_Cholesky_sets
    Map<MatrixXd> h2s,                  // n_RE x p
    Map<VectorXi> h2s_index,            // px1
    int grainSize) {

  int n = Y.rows();
  int p = Y.cols();
  int r = Z.cols();

  MatrixXd U(r,p);

  // process sets of columns with identical patterns of missing observations
  for(int col_set = 0; col_set < Missing_data_map.length(); col_set++){
    Rcpp::List Missing_data_map_set = Rcpp::as<Rcpp::List>(Missing_data_map[col_set]);
    Rcpp::List randomEffect_C_Choleskys_set = Rcpp::as<Rcpp::List>(randomEffect_C_Cholesky_sets[col_set]);

    VectorXi Y_obs = as<Map<VectorXi>>(Missing_data_map_set["Y_obs"]);
    VectorXi Y_cols = as<Map<VectorXi>>(Missing_data_map_set["Y_cols"]);

    int n_obs = Y_obs.size();
    int n_cols = Y_cols.size();

    // pull out elements of Y corresponding to observed individuals
    MatrixXd Y_set = get_block_matrix(Y,Y_obs,Y_cols);

    // draw random numbers
    MatrixXd randn_theta_set = rstdnorm_mat(r,n_cols);
    MatrixXd randn_e_set = rstdnorm_mat(n_obs,n_cols);

    VectorXi h2s_index_set(n_cols);
    VectorXd tot_Eta_prec_set(n_cols);
    MatrixXd h2s_set(h2s.rows(),n_cols);
    for(int j = 0; j < n_cols; j++){
      int col_j = Y_cols[j]-1;
      h2s_index_set[j] = h2s_index[col_j];
      tot_Eta_prec_set[j] = tot_Eta_prec[col_j];
      h2s_set.col(j) = h2s.col(col_j);
    }
    VectorXd h2_e_set = 1.0 - h2s_set.colwise().sum().array();
    VectorXd pes_set = tot_Eta_prec_set.array() / h2_e_set.array();

    std::vector<MSpMat> chol_C_list_set,chol_K_inv_list_set;
    for(int i = 0; i < h2s_index_set.maxCoeff(); i++){
      Rcpp::List randomEffect_C_Cholesky_i = Rcpp::as<Rcpp::List>(randomEffect_C_Choleskys_set[i]);
      chol_C_list_set.push_back(Rcpp::as<MSpMat>(randomEffect_C_Cholesky_i["chol_C"]));
      chol_K_inv_list_set.push_back(Rcpp::as<MSpMat>(randomEffect_C_Cholesky_i["chol_K_inv"]));
    }

    // sample U values for this set
    MatrixXd U_set(r,n_cols);
    sample_MME_single_diagR_worker sampler(Y_set,Z,chol_C_list_set,chol_K_inv_list_set,pes_set,tot_Eta_prec_set,h2s_index_set,randn_theta_set,randn_e_set,U_set);
    RcppParallel::parallelFor(0,p,sampler,grainSize);

    // copy new U values into full matrix
    for(int j = 0; j < n_cols; j++){
      U.col(Y_cols[j]-1) = U_set.col(j);
    }
  }
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
                         VectorXi h2s_index,                        // px1
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
    Map<VectorXi> h2s_index,      // px1
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

// calculates normal scores for:
// Y_resid %*% Sigma_inv %*% t(Y_resid)
// works on sets of columns of Y with same patterns of missingness
// [[Rcpp::export()]]
VectorXd tot_prec_scores_missing_c (
    Map<MatrixXd> Y,                // nxp
    Map<VectorXi> h2s_index,        // px1
    Rcpp::List Sigma_Cholesky_sets, // List. Each element contains: List of Cholesky decompositions of residuals. nxn upper triangular, dgCMatrix.
    Rcpp::List Missing_data_map,    // List. Each element containts Y_obs, Y_cols. The length should correspond to the length of Sigma_Choleksy_sets
    int grainSize
) {

  int p = Y.cols();
  VectorXd scores(p);


  for(int col_set = 0; col_set < Missing_data_map.size(); col_set++){
    Rcpp::List Missing_data_map_set = Rcpp::as<Rcpp::List>(Missing_data_map[col_set]);

    VectorXi Y_obs = as<Map<VectorXi>>(Missing_data_map_set["Y_obs"]);
    VectorXi Y_cols = as<Map<VectorXi>>(Missing_data_map_set["Y_cols"]);

    int n_obs = Y_obs.size();
    int n_cols = Y_cols.size();

    // pull out elements of Y corresponding to observed individuals
    MatrixXd Y_set = get_block_matrix(Y,Y_obs,Y_cols);

    // pull out h2s_index_set and tot_Eta_prec_set
    VectorXd h2s_index_set(n_cols);
    for(int j = 0; j < n_cols; j++){
      int col_j = Y_cols[j]-1;
      h2s_index_set[j] = h2s_index[col_j];
    }

    // pull out necessary chol_R matrices for set
    Rcpp::List Sigma_Choleskys_set = Rcpp::as<Rcpp::List>(Sigma_Cholesky_sets[col_set]);
    std::vector<MSpMat> chol_R_list_set;
    for(int i = 0; i < h2s_index_set.maxCoeff(); i++){
      Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys_set[i]);
      chol_R_list_set.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
    }

    VectorXd scores_set(n_cols);

    tot_prec_scores_worker sampler(Y_set,chol_R_list_set,h2s_index_set,scores_set);
    RcppParallel::parallelFor(0,n_cols,sampler,grainSize);

    // copy scores_set into scores
    for(int j = 0; j < n_cols; j++){
      scores(Y_cols[j]-1) = scores_set(j);
    }

  }
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
MatrixXd log_p_h2s_missing(
    Map<MatrixXd> Y,              // nxp
    Map<VectorXd> tot_Eta_prec,   // px1
    Rcpp::List Sigma_Cholesky_sets, // List. Each element contains: List. Each element contains: chol_Sigma:. nxn upper triangular, dgCMatrix, log_det_Sigma
    Rcpp::List Missing_data_map,    // List. Each element containts Y_obs, Y_cols. The length should correspond to the length of Sigma_Choleksy_sets
    Map<VectorXd> discrete_priors,// n_h2 x 1
    int grainSize)
{

  int b = discrete_priors.size();
  int p = Y.cols();
  MatrixXd log_ps(b,p);

  for(int col_set = 0; col_set < Missing_data_map.size(); col_set++){
    Rcpp::List Missing_data_map_set = Rcpp::as<Rcpp::List>(Missing_data_map[col_set]);

    VectorXi Y_obs = as<Map<VectorXi>>(Missing_data_map_set["Y_obs"]);
    VectorXi Y_cols = as<Map<VectorXi>>(Missing_data_map_set["Y_cols"]);

    int n_obs = Y_obs.size();
    int n_cols = Y_cols.size();

    // pull out elements of Y corresponding to observed individuals
    MatrixXd Y_set = get_block_matrix(Y,Y_obs,Y_cols);

    VectorXd tot_Eta_prec_set(n_cols);
    for(int j = 0; j < n_cols; j++){
      int col_j = Y_cols[j]-1;
      tot_Eta_prec_set[j] = tot_Eta_prec[col_j];
    }

    Rcpp::List Sigma_Choleskys_set = Rcpp::as<Rcpp::List>(Sigma_Cholesky_sets[col_set]);
    std::vector<MSpMat> chol_R_list_set;
    VectorXd log_det_Sigmas_set(b);
    for(int i = 0; i < b; i++){
      Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys_set[i]);
      chol_R_list_set.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
      log_det_Sigmas_set[i] = Rcpp::as<double>(Sigma_Choleskys_i["log_det"]);
    }

    MatrixXd log_ps_set(b,n_cols);

    log_ps_worker sampler(Y_set,tot_Eta_prec_set,chol_R_list_set,log_det_Sigmas_set,discrete_priors,log_ps_set);
    RcppParallel::parallelFor(0,b,sampler,grainSize);

    // copy new log_ps_set values into full matrix
    for(int j = 0; j < n_cols; j++){
      log_ps.col(Y_cols[j]-1) = log_ps_set.col(j);
    }
  }
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

  return h2s_index;
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
                                const VectorXi h2_index,                    // px1
                                const double step_size,                     // double
                                VectorXi &new_index                         // px1
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
    }
  }
};

// [[Rcpp::export()]]
VectorXi sample_h2s_discrete_MH_c(
    Map<MatrixXd> Y,                // nxp
    Map<VectorXd> tot_Eta_prec,     // px1
    Map<VectorXd> discrete_priors,  // n_h2 x 1
    VectorXi h2_index,              // px1
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


  sample_h2s_discrete_MH_worker sampler(Y,h2s_matrix,chol_R_list,log_det_Sigmas,tot_Eta_prec,discrete_priors,r_draws,state_draws,h2_index,step_size,new_index);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return new_index;
}

// [[Rcpp::export()]]
VectorXi sample_h2s_discrete_MH_missing_c(
    Map<MatrixXd> Y,                // nxp
    Map<VectorXd> tot_Eta_prec,     // px1
    Map<VectorXd> discrete_priors,  // n_h2 x 1
    VectorXi h2_index,              // px1
    Map<MatrixXd> h2s_matrix,       // n_RE x n_h2
    Rcpp::List Sigma_Cholesky_sets, // List. Each element contains: List. Each element contains: chol_Sigma:. nxn upper triangular, dgCMatrix, log_det_Sigma
    Rcpp::List Missing_data_map,    // List. Each element containts Y_obs, Y_cols. The length should correspond to the length of Sigma_Choleksy_sets
    double step_size,               // double
    int grainSize
){

  int p = Y.cols();
  int b = discrete_priors.size();
  VectorXi new_index(p);

  for(int col_set = 0; col_set < Missing_data_map.size(); col_set++){
    Rcpp::List Missing_data_map_set = Rcpp::as<Rcpp::List>(Missing_data_map[col_set]);

    VectorXi Y_obs = as<Map<VectorXi>>(Missing_data_map_set["Y_obs"]);
    VectorXi Y_cols = as<Map<VectorXi>>(Missing_data_map_set["Y_cols"]);

    int n_obs = Y_obs.size();
    int n_cols = Y_cols.size();

    // pull out elements of Y corresponding to observed individuals
    MatrixXd Y_set = get_block_matrix(Y,Y_obs,Y_cols);

    VectorXi h2_index_set(n_cols);
    VectorXd tot_Eta_prec_set(n_cols);
    for(int j = 0; j < n_cols; j++){
      int col_j = Y_cols[j]-1;
      h2_index_set[j] = h2_index[col_j];
      tot_Eta_prec_set[j] = tot_Eta_prec[col_j];
    }

    Rcpp::List Sigma_Choleskys_set = Rcpp::as<Rcpp::List>(Sigma_Cholesky_sets[col_set]);
    std::vector<MSpMat> chol_R_list_set;
    VectorXd log_det_Sigmas_set(b);
    for(int i = 0; i < b; i++){
      Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys_set[i]);
      chol_R_list_set.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
      log_det_Sigmas_set[i] = Rcpp::as<double>(Sigma_Choleskys_i["log_det"]);
    }

    VectorXd r_draws_set = as<VectorXd>(runif(n_cols));
    VectorXd state_draws_set = as<VectorXd>(runif(n_cols));
    VectorXi new_index_set(n_cols);

    sample_h2s_discrete_MH_worker sampler(Y_set,h2s_matrix,chol_R_list_set,log_det_Sigmas_set,tot_Eta_prec_set,discrete_priors,r_draws_set,state_draws_set,h2_index_set,step_size,new_index_set);
    RcppParallel::parallelFor(0,n_cols,sampler,grainSize);

    // copy new_index_set values into new_index
    for(int j = 0; j < n_cols; j++){
      new_index(Y_cols[j]-1) = new_index_set(j);
    }
  }

  return(new_index);
}
