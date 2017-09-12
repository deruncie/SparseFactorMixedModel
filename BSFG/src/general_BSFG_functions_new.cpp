#include <math.h>
#include <iostream>
#include "BSFG_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace RcppParallel;

// -------------------------------------------- //
// ---------- sample_MME_fixedEffects --------- //
// -------------------------------------------- //

// draws a sample of the vector b from the model:
// y = X %*% beta + e
// with beta[i] ~ N(prior_mean[i],1/prior_prec[i]), i=1:b
// with e ~ N(0,t(chol_R) %*% chol(R) / tot_Eta_prec), i=1:n
// Uses sampling method from MCMCglmm, which requires draws b+n draws from N(0,1), which are passed as randn_theta and randn_e
// If b >= n, inverts the C matrix using Binomial Inverse Theorem
// [[Rcpp::export()]]
VectorXd sample_MME_single_diagK(  // returns b x 1 vector
    VectorXd y,           // nx1
    MatrixXd X,           // nxb
    VectorXd prior_mean,  // bx1
    VectorXd prior_prec,  // bx1
    MSpMat chol_R,        // nxn upper triangular Cholesky decomposition of R. Format: dgCMatrix
    double tot_Eta_prec,  // double
    VectorXd randn_theta, // bx1
    VectorXd randn_e      // nx1
){
  int n = y.size();
  int b = X.cols();

  chol_R *= 1/sqrt(tot_Eta_prec);
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
      coefs.col(j) = sample_MME_single_diagK(Y.col(j), X, prior_mean.col(j), prior_prec.col(j), chol_R_list[h2_index], tot_Eta_prec[j], randn_theta.col(j),randn_e.col(j));
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
    Map<MatrixXd> Y,              // nxp
    Map<MatrixXd> X,              // nxb
    Rcpp::List Sigma_Cholesky_sets, // List. Each element contains: List of Cholesky decompositions of residuals. nxn upper triangular, dgCMatrix, Y_obs, Y_cols
    VectorXi h2s_index,           // px1 index of Cholesky matrix for each column
    Map<VectorXd> tot_Eta_prec,   // px1
    Map<MatrixXd> prior_mean,     // bxp
    Map<MatrixXd> prior_prec,     // bxp
    int grainSize) {

  int n = X.rows();
  int b = X.cols();
  int p = Y.cols();

  MatrixXd coefs(b,p);

  // sample sets of columns with same pattern of missing data as a block
  for(int col_set = 0; col_set < Sigma_Cholesky_sets.length(); col_set++){
    Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Cholesky_sets[col_set]);

    VectorXi Y_obs = as<Map<VectorXi>>(Sigma_Choleskys_i["Y_obs"]);
    VectorXi Y_cols = as<Map<VectorXi>>(Sigma_Choleskys_i["Y_cols"]);

    int n_obs = Y_obs.size();
    int n_cols = Y_cols.size();

    // pull out rows of X corresponding to observed individuals
    MatrixXd X_set(n_obs,b);
    for(int i = 0; i < n_obs; i++) {
      X_set.row(i) = X.row(Y_obs[i]-1);
    }

    // pull out elements of Y corresponding to observed individuals
    MatrixXd Y_set(n_obs,n_cols);
    for(int i = 0; i < n_obs; i++){
      int row_i = Y_obs[i]-1;
      for(int j = 0; j < n_cols; j++){
        int col_j = Y_cols[j]-1;
        Y_set.coeffRef(i,j) = Y.coeffRef(row_i,col_j);
      }
    }

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
    Rcpp::List Sigma_Choleskys = Rcpp::as<Rcpp::List>(Sigma_Choleskys_i["Sigma_Choleskys"]);
    std::vector<MSpMat> chol_R_list_set;
    for(int i = 0; i < h2s_index_set.maxCoeff(); i++){
      Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);
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


