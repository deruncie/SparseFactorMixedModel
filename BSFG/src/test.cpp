#include <math.h>
#include <iostream>
#include "BSFG_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace RcppParallel;

// [[Rcpp::export()]]
MatrixXd XDXt_c(Map<MatrixXd> X, Map<VectorXd> d){
  return(X * d.asDiagonal() * X.transpose());
}

// [[Rcpp::export()]]
VectorXd forwardSolve_c(Map<MatrixXd> chol_M, Map<VectorXd> y){
  return(chol_M.transpose().triangularView<Lower>().solve(y));
}

// [[Rcpp::export()]]
MatrixXd chol_c(Map<MatrixXd> C){
    LLT<MatrixXd> C_llt;
    C_llt.compute(C);
    return(C_llt.matrixU());
}

//
//
// // [[Rcpp::export()]]
// VectorXd sample_MME_single_diagK2(  // returns b x 1 vector
//     VectorXd y,           // nx1
//     MatrixXd X,           // nxb
//     VectorXd prior_mean,  // bx1
//     VectorXd prior_prec,  // bx1
//     MSpMat chol_R,        // nxn upper triangular Cholesky decomposition of R. Format: dgCMatrix
//     double tot_Eta_prec, // double
//     VectorXd randn_theta, // bx1
//     VectorXd randn_e      // 0x1 or nx1. 0x1 if b<n
// ){
//   if(randn_e.size() == 0){
//     MatrixXd RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
//     VectorXd XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y * sqrt(tot_Eta_prec));
//     VectorXd XtRinvy_std_mu = XtRinvy + prior_prec.asDiagonal()*prior_mean;
//     MatrixXd C = RinvSqX.transpose() * RinvSqX;
//     C.diagonal() += prior_prec;
//     LLT<MatrixXd> C_llt;
//     C_llt.compute(C);
//     MatrixXd chol_C = C_llt.matrixU();
//
//     VectorXd b = chol_C.transpose().triangularView<Lower>().solve(XtRinvy_std_mu);
//     b += randn_theta;
//     b = chol_C.triangularView<Upper>().solve(b);
//     return(b);
//   } else {
//     // Using algorithm from Bhattacharya et al 2016 Biometrika. https://academic.oup.com/biomet/article/103/4/985/2447851
//     MatrixXd Phi = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
//     VectorXd alpha = chol_R.transpose().triangularView<Lower>().solve(y * sqrt(tot_Eta_prec));
//
//     VectorXd u = randn_theta.array() / prior_prec.cwiseSqrt().array();
//     u += prior_mean;
//     VectorXd v = Phi * u + randn_e;
//     VectorXd alpha_v = alpha-v;
//
//     MatrixXd D_PhiT = prior_prec.cwiseInverse().asDiagonal() * Phi.transpose();
//     MatrixXd cov = Phi * D_PhiT;
//     cov.diagonal().array() += 1.0;
//
//     VectorXd w = cov.ldlt().solve(alpha_v);
//
//     VectorXd theta = u + D_PhiT * w;
//
//     return(theta);
//   }
// }
//
//
// // [[Rcpp::export()]]
// VectorXd sample_MME_single_diagK2b(  // returns b x 1 vector
//     VectorXd y,           // nx1
//     MatrixXd X,           // nxb
//     VectorXd prior_mean,  // bx1
//     VectorXd prior_prec,  // bx1
//     MSpMat inv_chol_Rt,        // nxn upper triangular Cholesky decomposition of R. Format: dgCMatrix
//     double tot_Eta_prec, // double
//     VectorXd randn_theta, // bx1
//     VectorXd randn_e      // 0x1 or nx1. 0x1 if b<n
// ){
//   if(randn_e.size() == 0){
//     MatrixXd RinvSqX = inv_chol_Rt * X * sqrt(tot_Eta_prec);
//     VectorXd XtRinvy = RinvSqX.transpose() * (inv_chol_Rt * (y * sqrt(tot_Eta_prec)));
//     VectorXd XtRinvy_std_mu = XtRinvy + prior_prec.asDiagonal()*prior_mean;
//     MatrixXd C = RinvSqX.transpose() * RinvSqX;
//     C.diagonal() += prior_prec;
//     LLT<MatrixXd> C_llt;
//     C_llt.compute(C);
//     MatrixXd chol_C = C_llt.matrixU();
//
//     VectorXd b = chol_C.transpose().triangularView<Lower>().solve(XtRinvy_std_mu);
//     b += randn_theta;
//     b = chol_C.triangularView<Upper>().solve(b);
//     return(b);
//   } else {
//     // Using algorithm from Bhattacharya et al 2016 Biometrika. https://academic.oup.com/biomet/article/103/4/985/2447851
//     MatrixXd Phi = inv_chol_Rt * X * sqrt(tot_Eta_prec);
//     VectorXd alpha = inv_chol_Rt * y * sqrt(tot_Eta_prec);
//
//     VectorXd u = randn_theta.array() / prior_prec.cwiseSqrt().array();
//     u += prior_mean;
//     VectorXd v = Phi * u + randn_e;
//     VectorXd alpha_v = alpha-v;
//
//     MatrixXd D_PhiT = prior_prec.cwiseInverse().asDiagonal() * Phi.transpose();
//     MatrixXd cov = Phi * D_PhiT;
//     cov.diagonal().array() += 1.0;
//
//     VectorXd w = cov.ldlt().solve(alpha_v);
//
//     VectorXd theta = u + D_PhiT * w;
//
//     return(theta);
//   }
// }
//
//
// // [[Rcpp::export()]]
// VectorXd sample_MME_single_diagK3(  // returns b x 1 vector
//     VectorXd y,           // nx1
//     MatrixXd X,           // nxb
//     VectorXd prior_mean,  // bx1
//     VectorXd prior_prec,  // bx1
//     MSpMat chol_R,        // nxn upper triangular Cholesky decomposition of R. Format: dgCMatrix
//     double tot_Eta_prec, // double
//     VectorXd randn_theta, // bx1
//     VectorXd randn_e      // 0x1 or nx1. 0x1 if b<n
// ){
//   if(randn_e.size() == 0){
//     MatrixXd RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
//     VectorXd XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y * sqrt(tot_Eta_prec));
//     VectorXd XtRinvy_std_mu = XtRinvy + prior_prec.asDiagonal()*prior_mean;
//     MatrixXd C = RinvSqX.transpose() * RinvSqX;
//     C.diagonal() += prior_prec;
//     LLT<MatrixXd> C_llt;
//     C_llt.compute(C);
//     MatrixXd chol_C = C_llt.matrixU();
//
//     VectorXd b = chol_C.transpose().triangularView<Lower>().solve(XtRinvy_std_mu);
//     b += randn_theta;
//     b = chol_C.triangularView<Upper>().solve(b);
//     return(b);
//   } else {
//     // should check that randn_e.size() == n
//     VectorXd theta_star = randn_theta.array() / prior_prec.cwiseSqrt().array();
//     theta_star += prior_mean;
//     VectorXd e_star = chol_R.transpose() * (randn_e / sqrt(tot_Eta_prec));
//     MatrixXd X_theta_star = X * theta_star;
//     VectorXd y_resid = y - X_theta_star - e_star;
//
//     MatrixXd RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
//     VectorXd XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y_resid * sqrt(tot_Eta_prec));
//
//     VectorXd theta_tilda;
//     MatrixXd VAi = X * prior_prec.cwiseInverse().asDiagonal();
//     MatrixXd inner = VAi*X.transpose() + (chol_R.transpose() * chol_R) / tot_Eta_prec;
//     VectorXd VAiXtURinvy = VAi * XtRinvy;
//     VectorXd outerXtURinvy = VAi.transpose() * inner.ldlt().solve(VAiXtURinvy);
//     theta_tilda = XtRinvy.array() / prior_prec.array();
//     theta_tilda -= outerXtURinvy;
//     VectorXd theta = theta_star + theta_tilda;
//
//     return(theta);
//   }
// }


//
// // [[Rcpp::export()]]
// VectorXd cumprod_c(VectorXd x) {
//   int n = x.size();
//   VectorXd y(n);
//   y(0) = x[0];
//   if(n>1){
//     for(int i = 1; i < n; i++){
//       y[i] = y[i-1]*x[i];
//     }
//   }
//   return(y);
// }
//
// // [[Rcpp::export()]]
// VectorXd multiply_vec(VectorXd xx, int n, double y){
//   VectorXd x(xx);
//   x.tail(x.size()-n) *= y;
//   return(x);
// }
//
// // [[Rcpp::export()]]
// VectorXd sample_delta_c_Eigen2(
//     VectorXd delta,
//     VectorXd tauh,
//     Map<VectorXd> scores,
//     double delta_1_rate,
//     double delta_2_rate,
//     Map<MatrixXd> randg_draws  // all done with rate = 1;
// ) {
//   int times = randg_draws.rows();
//   int k = tauh.size();
//
//   double rate,delta_old;
//   for(int i = 0; i < times; i++){
//     delta_old = delta(0);
//     rate = delta_1_rate + (1/delta(0)) * tauh.dot(scores);
//     delta(0) = randg_draws(i,0) / rate;
//     // tauh = cumprod(delta);
//     tauh *= delta(0)/delta_old;   // replaces re-calculating cumprod
//
//     for(int h = 1; h < k; h++) {
//       delta_old = delta(h);
//       rate = delta_2_rate + (1/delta(h))*tauh.tail(k-h).dot(scores.tail(k-h));
//       delta(h) = randg_draws(i,h) / rate;
//       // tauh = cumprod(delta);
//       tauh.tail(k-h) *= delta(h)/delta_old; // replaces re-calculating cumprod
//       // Rcout << (tauh - cumprod(delta)).sum() << std::endl;
//     }
//     // Rcout << tauh.transpose() << std::endl;
//     // Rcout << delta.transpose() << std::endl;
//     // Rcout << cumprod_c(delta).transpose() << std::endl;
//     // Rcout << i << " " << (tauh - cumprod_c(delta)).sum() << std::endl;
//   }
//   return(delta);
// }


//
//
// // [[Rcpp::export()]]
// MatrixXd rstdnorm_mat2(int n,int p) {  // returns nxp matrix
//   VectorXd X_vec = as<VectorXd>(rnorm(n*p));
//   Map<MatrixXd> X_mat(X_vec.data(),n,p);
//   return(X_mat);
// }
//
// // [[Rcpp::export()]]
// VectorXd sample_MME_single_hierarchical_diagK(  // returns b x 1 vector
//     VectorXd y,           // nx1
//     MSpMat   Z,           // nxr   // use when r < n < b
//     MatrixXd X,           // rxb
//     VectorXd prior_mean,  // bx1
//     VectorXd prior_prec,  // bx1
//     MSpMat chol_R,        // nxn upper triangular Cholesky decomposition of R. Format: dgCMatrix
//     double tot_Eta_prec, // double
//     VectorXd randn_theta, // bx1
//     VectorXd randn_e      // 0x1 or nx1. 0x1 if b<n
// ){
//   VectorXd theta_star = randn_theta.array() / prior_prec.cwiseSqrt().array();
//   theta_star += prior_mean;
//   VectorXd e_star = chol_R.transpose() * (randn_e / sqrt(tot_Eta_prec));
//   MatrixXd ZX = Z*X;
//   MatrixXd ZX_theta_star = ZX * theta_star;
//   VectorXd y_resid = y - ZX_theta_star - e_star;
//
//   MatrixXd RinvSqZ = chol_R.transpose().triangularView<Lower>().solve(Z.toDense() * sqrt(tot_Eta_prec));
//   MatrixXd RinvSqZX = RinvSqZ*X;
//   VectorXd XtZtRinvy = RinvSqZX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y_resid * sqrt(tot_Eta_prec));
//
//   VectorXd theta_tilda;
//   MatrixXd VAi = X * prior_prec.cwiseInverse().asDiagonal();
//   MatrixXd I(Z.cols(),Z.cols());
//   I.setIdentity();
//   MatrixXd ZtRinvZ = RinvSqZ.transpose()*RinvSqZ;
//   MatrixXd inner = VAi*X.transpose() + ZtRinvZ.ldlt().solve(I);
//   VectorXd VAiXtURinvy = VAi * XtZtRinvy;
//   VectorXd outerXtURinvy = VAi.transpose() * inner.ldlt().solve(VAiXtURinvy);
//   theta_tilda = XtZtRinvy.array() / prior_prec.array();
//   theta_tilda -= outerXtURinvy;
//   VectorXd theta = theta_star + theta_tilda;
//
//   return(theta);
// }
//
// struct sample_MME_single_hierarchical_diagK_worker : public RcppParallel::Worker {
//   MatrixXd Y;
//   MSpMat Z;
//   MatrixXd X;
//   MatrixXd prior_mean, prior_prec, randn_theta, randn_e;
//   const std::vector<MSpMat> chol_R_list;
//   VectorXi h2s_index;
//   VectorXd tot_Eta_prec;
//   MatrixXd &coefs;
//
//   sample_MME_single_hierarchical_diagK_worker(
//     MatrixXd Y,           // nxp
//     MSpMat   Z,           // nxb
//     MatrixXd X,           // nxb
//     MatrixXd prior_mean,  // bxp
//     MatrixXd prior_prec,  // bxp
//     const std::vector<MSpMat> &chol_R_list, // each element a MSpMat nxn upper-triangular cholesky decomposition of R
//     VectorXi h2s_index,   // px1, 1-based index
//     VectorXd tot_Eta_prec,// px1
//     MatrixXd randn_theta, // bxp
//     MatrixXd randn_e,     // 0xp or nxp. The former if b<n, the latter if b >= n
//     MatrixXd &coefs       // bxp
//   ):
//     Y(Y), Z(Z), X(X), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
//     chol_R_list(chol_R_list), h2s_index(h2s_index), tot_Eta_prec(tot_Eta_prec), coefs(coefs) {}
//
//   void operator()(std::size_t begin, std::size_t end) {
//     for(std::size_t j = begin; j < end; j++){
//       int h2_index = h2s_index[j] - 1;
//       MSpMat chol_R = chol_R_list[h2_index];
//       coefs.col(j) = sample_MME_single_hierarchical_diagK(Y.col(j), Z, X, prior_mean.col(j), prior_prec.col(j), chol_R, tot_Eta_prec[j], randn_theta.col(j),randn_e.col(j));
//     }
//   }
// };
//
// // Samples from B in model:
// // Y = ZXB + E
// // [[Rcpp::export()]]
// MatrixXd sample_MME_fixedEffects_hierarchical_c(  // returns bxp matrix
//     Map<MatrixXd> Y,              // nxp
//     MSpMat        Z,              // nxr   // use when r < n < b
//     Map<MatrixXd> X,              // rxb
//     Rcpp::List Sigma_Choleskys,   // list of Cholesky decompositions of residuals. nxn upper triangular, dgCMatrix
//     VectorXi h2s_index,           // px1 index of Cholesky matrix for each column
//     Map<VectorXd> tot_Eta_prec,   // px1
//     Map<MatrixXd> prior_mean,     // bxp
//     Map<MatrixXd> prior_prec,     // bxp
//     int grainSize) {
//
//   int b = X.cols();
//   int p = Y.cols();
//   int n = Y.rows();
//
//   MatrixXd randn_theta = rstdnorm_mat2(b,p);
//   MatrixXd randn_e = rstdnorm_mat2(n,p);
//
//   std::vector<MSpMat> chol_R_list;
//   for(int i = 0; i < h2s_index.maxCoeff(); i++){
//     Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);
//     chol_R_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
//   }
//
//   MatrixXd coefs(b,p);
//
//   sample_MME_single_hierarchical_diagK_worker sampler(Y,Z,X,prior_mean,prior_prec,chol_R_list,h2s_index,tot_Eta_prec,randn_theta,randn_e, coefs);
//   RcppParallel::parallelFor(0,p,sampler,grainSize);
//
//   return(coefs);
// }


// //
// // [[Rcpp::export()]]
// MatrixXd get_fitted_set_c2(  // returns n_tot x p matrix in same order as data
//     Rcpp::List model_matrices,  // List. Each element contains: y (n_i x t), X (n_i x b), s (n_i x 1), position (n_i x 1)
//     Map<MatrixXd> coefs,  // b x n matrix
//     int grainSize){
//
//   std::vector<MatrixXd> X_list;
//   std::vector<ArrayXi> position_list;
//   std::vector<ArrayXi> nonZero_cols_X;
//   int total_obs = 0;
//   int n = model_matrices.size();
//   for(int i = 0; i < n; i++){
//     Rcpp::List model_matrix_i = Rcpp::as<Rcpp::List>(model_matrices[i]);
//     X_list.push_back(Rcpp::as<MatrixXd>(model_matrix_i["X"]));
//     position_list.push_back(Rcpp::as<ArrayXi>(model_matrix_i["position"]));
//     nonZero_cols_X.push_back(Rcpp::as<ArrayXi>(model_matrix_i["nonZero_cols_X"])); // list of which columns b_X correspond to in full X matrix
//
//     int n_obs = X_list[i].rows();
//     total_obs += n_obs;
//   }
//
//   Rcpp::List model_matrix_1 = Rcpp::as<Rcpp::List>(model_matrices[0]);
//   MatrixXd Y = Rcpp::as<MatrixXd>(model_matrix_1["y"]);
//   int n_traits = Y.cols();
//
//   MatrixXd Y_fitted(total_obs,n_traits);
//
//   struct sampleColumn : public RcppParallel::Worker {
//     std::vector<MatrixXd> X_list;
//     std::vector<ArrayXi> position_list;
//     std::vector<ArrayXi> nonZero_cols_X;
//     int n_traits;
//     MatrixXd coefs;
//     MatrixXd &Y_fitted;
//
//     sampleColumn(
//       std::vector<MatrixXd> &X_list,
//       std::vector<ArrayXi>  &position_list,
//       std::vector<ArrayXi>  &nonZero_cols_X,
//       int n_traits,
//       MatrixXd coefs,
//       MatrixXd &Y_fitted) :
//       X_list(X_list), position_list(position_list),nonZero_cols_X(nonZero_cols_X),
//       n_traits(n_traits),
//       coefs(coefs),Y_fitted(Y_fitted)
//     {}
//
//     void operator()(std::size_t begin, std::size_t end) {
//       for(std::size_t j = begin; j < end; j++){
//         int b_X = X_list[j].cols();
//         int b_tot = coefs.rows() / n_traits;
//         VectorXd coefs_j(b_X * n_traits);
//         for(int t = 0; t < n_traits; t++){
//           for(int k = 0; k < b_X; k++){
//             coefs_j[t*b_X+k] = coefs.coeffRef(t*b_tot+nonZero_cols_X[j][k]-1,j);
//           }
//         }
//         Map<MatrixXd> Eta_i(coefs_j.data(),b_X,n_traits);
//         MatrixXd Y_fitted_j = X_list[j] * Eta_i;
//         for(int i = 0; i < position_list[j].size(); i++) {
//           Y_fitted.row(position_list[j][i]-1) = Y_fitted_j.row(i);
//         }
//       }
//     }
//   };
//
//   sampleColumn sampler(X_list,position_list,nonZero_cols_X,n_traits,coefs,Y_fitted);
//   RcppParallel::parallelFor(0,n,sampler,grainSize);
//
//   return(Y_fitted);
// }

//
// // [[Rcpp::export()]]
// MatrixXd rstdnorm_mat2(int n,int p) {  // returns nxp matrix
//   VectorXd X_vec = as<VectorXd>(rnorm(n*p));
//   Map<MatrixXd> X_mat(X_vec.data(),n,p);
//   return(X_mat);
// }
//
// // [[Rcpp::export()]]
// MatrixXd slice(const MatrixXd &X, VectorXi rows, VectorXi cols){
//   int n_rows = rows.size();
//   int n_cols = cols.size();
//   MatrixXd X_sub(n_rows,n_cols);
//   for(int j = 0; j < n_cols; j++){
//     for(int i = 0; i < n_rows; i++){
//       X_sub.coeffRef(i,j) = X.coeffRef(rows[i]-1,cols[j]-1);
//     }
//   }
//   return X_sub;
// }
//
// // [[Rcpp::export()]]
// MatrixXd mat_slice(int m, int n, VectorXi rows, VectorXi cols){
//   MatrixXd X = rstdnorm_mat2(m,n);
//   return(slice(X,rows,cols));
// }
//
// // [[Rcpp::export()]]
// VectorXd sample_MME_single_diagK2(  // returns b x 1 vector
//     VectorXd y,           // nx1
//     MatrixXd X,           // nxb
//     VectorXd prior_mean,  // bx1
//     VectorXd prior_prec,  // bx1
//     MSpMat chol_R,        // nxn upper triangular Cholesky decomposition of R. Format: dgCMatrix
//     double tot_Eta_prec, // double
//     VectorXd randn_theta, // bx1
//     VectorXd randn_e      // 0x1 or nx1. 0x1 if b<n
// ){
//   if(randn_e.size() == 0){
//     MatrixXd RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
//     VectorXd XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y * sqrt(tot_Eta_prec));
//     VectorXd XtRinvy_std_mu = XtRinvy + prior_prec.asDiagonal()*prior_mean;
//     MatrixXd C = RinvSqX.transpose() * RinvSqX;
//     C.diagonal() += prior_prec;
//     LLT<MatrixXd> C_llt;
//     C_llt.compute(C);
//     MatrixXd chol_C = C_llt.matrixU();
//
//     VectorXd b = chol_C.transpose().triangularView<Lower>().solve(XtRinvy_std_mu);
//     b += randn_theta;
//     b = chol_C.triangularView<Upper>().solve(b);
//     return(b);
//   } else {
//     // should check that randn_e.size() == n
//     VectorXd theta_star = randn_theta.array() / prior_prec.cwiseSqrt().array();
//     theta_star += prior_mean;
//     VectorXd e_star = chol_R.transpose() * (randn_e / sqrt(tot_Eta_prec));
//     MatrixXd X_theta_star = X * theta_star;
//     VectorXd y_resid = y - X_theta_star - e_star;
//
//     MatrixXd RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
//     VectorXd XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y_resid * sqrt(tot_Eta_prec));
//
//     VectorXd theta_tilda;
//     MatrixXd VAi = X * prior_prec.cwiseInverse().asDiagonal();
//     MatrixXd inner = VAi*X.transpose() + (chol_R.transpose() * chol_R) / tot_Eta_prec;
//     VectorXd VAiXtURinvy = VAi * XtRinvy;
//     VectorXd outerXtURinvy = VAi.transpose() * inner.ldlt().solve(VAiXtURinvy);
//     theta_tilda = XtRinvy.array() / prior_prec.array();
//     theta_tilda -= outerXtURinvy;
//     VectorXd theta = theta_star + theta_tilda;
//
//     return(theta);
//   }
// }
//
//
// // [[Rcpp::export()]]
// MatrixXd sample_coefs_set_c2(    // return pxn matrix
//     Rcpp::List model_matrices,  // List. Each element contains: y (n_i x t), X (n_i x p), position (n_i x 1)
//     Map<VectorXd> tot_Y_prec,   // tx1
//     Map<MatrixXd> prior_mean,   // pxn
//     Map<MatrixXd> prior_prec,   // pxn
//     int grainSize){
//
//   int n = model_matrices.size();
//   int p = prior_mean.rows();
//
//   std::vector<MatrixXd> y_list;
//   std::vector<MatrixXd> X_list;
//   std::vector<ArrayXi> nonZero_cols_X;
//   std::vector<MatrixXd> randn_theta_list;
//   for(int i = 0; i < n; i++){
//     Rcpp::List model_matrix_i = Rcpp::as<Rcpp::List>(model_matrices[i]);
//     y_list.push_back(Rcpp::as<MatrixXd>(model_matrix_i["y"]));
//     X_list.push_back(Rcpp::as<MatrixXd>(model_matrix_i["X"]));
//     nonZero_cols_X.push_back(Rcpp::as<ArrayXi>(model_matrix_i["nonZero_cols_X"]));
//
//     int t = y_list[i].cols();
//     MatrixXd randn_theta = rstdnorm_mat2(p/t,t);
//     randn_theta_list.push_back(randn_theta);
//   }
//
//   int n_traits = y_list[0].cols();
//
//   MatrixXd coefs(p,n);
//
//   struct sampleColumn : public RcppParallel::Worker {
//     std::vector<MatrixXd> y_list;
//     std::vector<MatrixXd> X_list;
//     std::vector<ArrayXi> nonZero_cols_X;
//     std::vector<MatrixXd> randn_theta_list;
//     VectorXd tot_Y_prec;
//     MatrixXd prior_mean,prior_prec;
//     int n_traits;
//     MatrixXd &coefs;
//
//     sampleColumn(
//       std::vector<MatrixXd> &y_list,
//       std::vector<MatrixXd> &X_list,
//       std::vector<ArrayXi>  &nonZero_cols_X,
//       std::vector<MatrixXd> &randn_theta_list,
//       VectorXd tot_Y_prec,
//       MatrixXd prior_mean,
//       MatrixXd prior_prec,
//       int n_traits,
//       MatrixXd &coefs) :
//       y_list(y_list), X_list(X_list), nonZero_cols_X(nonZero_cols_X), randn_theta_list(randn_theta_list),
//       tot_Y_prec(tot_Y_prec),prior_mean(prior_mean), prior_prec(prior_prec),n_traits(n_traits),
//       coefs(coefs)
//     {}
//
//     void operator()(std::size_t begin, std::size_t end) {
//       for(std::size_t j = begin; j < end; j++){
//         MatrixXd Y = y_list[j];
//         MatrixXd X = X_list[j];
//         int b = randn_theta_list[j].rows();
//         int b_X = X.cols();
//         int n_obs = Y.rows();
//         SpMat I = MatrixXd::Identity(n_obs,n_obs).sparseView();
//         MSpMat chol_R(I.rows(),I.cols(), I.nonZeros(),I.outerIndexPtr(),I.innerIndexPtr(),I.valuePtr());
//         MatrixXd randn_e = MatrixXd::Zero(0,n_traits);
//         for(int t = 0; t < n_traits; t++) {
//           // first assign the result vector to prior_mean + randn/sqrt(prec)
//           // will then replace values with sampled values.
//           VectorXd prior_mean_tj = prior_mean.block(t*b,j,b,1);
//           VectorXd prior_prec_tj = prior_prec.block(t*b,j,b,1);
//           VectorXd randn_theta_tj = randn_theta_list[j].col(t);
//           coefs.block(t*b,j,b,1) = prior_mean_tj.array() + randn_theta_tj.array() / prior_prec_tj.array().sqrt();
//
//           // now, pull out parameters for the coefficients corresponding to the columns of X
//           VectorXd prior_mean_tj_X(b_X);
//           VectorXd prior_prec_tj_X(b_X);
//           VectorXd randn_theta_tj_X(b_X);
//           for(int k = 0; k < b_X; k++){
//             int element = nonZero_cols_X[j][k]-1;
//             prior_mean_tj_X.coeffRef(k) = prior_mean_tj.coeffRef(element);
//             prior_prec_tj_X.coeffRef(k) = prior_prec_tj.coeffRef(element);
//             randn_theta_tj_X.coeffRef(k) = randn_theta_tj.coeffRef(element);
//           }
//           VectorXd coefs_X = sample_MME_single_diagK2(Y.col(t), X,
//                                                       prior_mean_tj_X, prior_prec_tj_X,
//                                                       chol_R,tot_Y_prec[t],
//                                                       randn_theta_tj_X,randn_e.col(t));
//
//           // now replace the values in coef with the corresponding ones in coefs_X
//           for(int k = 0; k < b_X; k++){
//             int element = nonZero_cols_X[j][k]-1;
//             coefs.coeffRef(t*b + element,j) = coefs_X.coeffRef(k);
//           }
//         }
//       }
//     }
//   };
//
//
//   sampleColumn sampler(y_list, X_list,nonZero_cols_X,randn_theta_list,tot_Y_prec,prior_mean,prior_prec,n_traits,coefs);
//   RcppParallel::parallelFor(0,n,sampler,grainSize);
//
//   return(coefs);
// }



//
//
// // [[Rcpp::export()]]
// VectorXd sample_MME_single_diagK2(  // returns b x 1 vector
//     VectorXd y,           // nx1
//     MatrixXd X,           // nxb
//     VectorXd prior_mean,  // bx1
//     VectorXd prior_prec,  // bx1
//     MSpMat chol_R,        // nxn upper triangular Cholesky decomposition of R. Format: dgCMatrix
//     double tot_Eta_prec, // double
//     VectorXd randn_theta, // bx1
//     VectorXd randn_e,      // 0x1 or nx1. 0x1 if b<n
//     bool chol_R_triangular = true
// ){
//   if(randn_e.size() == 0){
//     MatrixXd C;
//     VectorXd XtRinvy;
//     if(chol_R_triangular == true){
//       MatrixXd RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
//       XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y * sqrt(tot_Eta_prec));
//       C = RinvSqX.transpose() * RinvSqX;
//     } else{
//       MatrixXd R = (chol_R.transpose() * chol_R) / tot_Eta_prec;
//       MatrixXd RinvX = R.colPivHouseholderQr().solve(X);
//       C = X.transpose() * RinvX;
//       XtRinvy = RinvX.transpose() * y;
//     }
//     // Rcout << XtRinvy.transpose() << std::endl;
//     // Rcout << C.row(0) << std::endl;
//     VectorXd XtRinvy_std_mu = XtRinvy + prior_prec.asDiagonal()*prior_mean;
//     C.diagonal() += prior_prec;
//     LLT<MatrixXd> C_llt;
//     C_llt.compute(C);
//     MatrixXd chol_C = C_llt.matrixU();
//
//     VectorXd b = chol_C.transpose().triangularView<Lower>().solve(XtRinvy_std_mu);
//     b += randn_theta;
//     b = chol_C.triangularView<Upper>().solve(b);
//     return(b);
//   } else {
//     // should check that randn_e.size() == n
//     VectorXd theta_star = randn_theta.array() / prior_prec.cwiseSqrt().array();
//     theta_star += prior_mean;
//     VectorXd e_star = chol_R.transpose() * (randn_e / sqrt(tot_Eta_prec));
//     MatrixXd X_theta_star = X * theta_star;
//     VectorXd y_resid = y - X_theta_star - e_star;
//
//     MatrixXd R = (chol_R.transpose() * chol_R) / tot_Eta_prec;
//
//     VectorXd XtRinvy;
//     if(chol_R_triangular == true) {
//       MatrixXd RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
//       XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y_resid * sqrt(tot_Eta_prec));
//     } else{
//       XtRinvy = X.transpose() * R.ldlt().solve(y_resid);
//     }
//
//     VectorXd theta_tilda;
//     MatrixXd VAi = X * prior_prec.cwiseInverse().asDiagonal();
//     MatrixXd inner = VAi*X.transpose() + R;
//     VectorXd VAiXtURinvy = VAi * XtRinvy;
//     VectorXd outerXtURinvy = VAi.transpose() * inner.ldlt().solve(VAiXtURinvy);
//     theta_tilda = XtRinvy.array() / prior_prec.array();
//     theta_tilda -= outerXtURinvy;
//     VectorXd theta = theta_star + theta_tilda;
//
//     return(theta);
//   }
// }
//
// // [[Rcpp::export()]]
// VectorXd sample_coefs_single_hierarchical(
//     VectorXd y,         // nx1
//     MatrixXd ZX,
//     MatrixXd Z,              // nxr
//     MatrixXd X,             // nxp
//     VectorXd prior_mean,    // bx1
//     VectorXd prior_prec,    // bx1
//     MSpMat chol_R,
//     double tot_Eta_prec,
//     VectorXd randn_theta,   // bx1
//     VectorXd randn_e       // nx1
// ) {
//   if(randn_e.size() == 0){
//     return(sample_MME_single_diagK2(y,ZX,prior_mean,prior_prec,chol_R,tot_Eta_prec,randn_theta,randn_e));
//   } else {
//     // should check that randn_e.size() == n
//     VectorXd theta_star = randn_theta.array() / prior_prec.cwiseSqrt().array();
//     theta_star += prior_mean;
//     VectorXd e_star = chol_R.transpose() * (randn_e / sqrt(tot_Eta_prec));
//     MatrixXd ZX_theta_star = ZX * theta_star;
//     VectorXd y_resid = y - ZX_theta_star - e_star;
//
//     MatrixXd RinvSqZX = chol_R.transpose().triangularView<Lower>().solve(ZX * sqrt(tot_Eta_prec));
//     VectorXd XtZtRinvy = RinvSqZX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y_resid * sqrt(tot_Eta_prec));
//
//     int r = X.rows();
//     MatrixXd R = (chol_R.transpose() * chol_R) / tot_Eta_prec;
//
//     MatrixXd B = Z.transpose() * R.llt().solve(Z);
//     SimplicialLLT<SparseMatrix<double> >solver;
//     // MatrixXd I(r,r);
//     // I.setIdentity();
//     // MatrixXd Bi = solver.compute(B).solve(MatrixXd::Identity(r,r));
//     MatrixXd Bi = B.ldlt().solve(MatrixXd::Identity(r,r));
//     MatrixXd VAi = X * prior_prec.cwiseInverse().asDiagonal();
//     MatrixXd inner = VAi * X.transpose() + Bi;
//     VectorXd VAiXtZtRinvy = VAi * XtZtRinvy;
//     VectorXd outerXtWtQRinvy = VAi.transpose() * inner.ldlt().solve(VAiXtZtRinvy);
//     VectorXd theta_tilda = XtZtRinvy.array() / prior_prec.array();
//     theta_tilda -= outerXtWtQRinvy;
//     VectorXd coefs = theta_tilda + theta_star;
//     return(coefs);
//   }
// }



// // [[Rcpp::export()]]
// void sample_MME_fixedEffects_cis_c2(
//     Map<MatrixXd> Y,
//     Map<MatrixXd> X,
//     Rcpp::List cis_genotypes,
//     Rcpp::List Sigma_Choleskys,
//     VectorXi h2s_index,
//     Map<VectorXd> tot_Eta_prec,
//     Map<MatrixXd> prior_mean,
//     Map<MatrixXd> prior_prec,
//     Map<VectorXd> cis_effect_index,
//     int total_cis_effects,
//     int grainSize) {
//
//   int p = Y.cols();
//   int b = X.cols();
//
//   std::vector<MSpMat> chol_R_list;
//   for(int i = 0; i < h2s_index.maxCoeff(); i++){
//     Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);
//     chol_R_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
//   }
//
//   std::vector<MatrixXd> cis_X;
//   for(int i = 0; i < p; i++){
//     cis_X.push_back(Rcpp::as<MatrixXd>(cis_genotypes[i]));
//   }
//   //
//   // MatrixXd coefs(b,p);
//   // VectorXd cis_effects(total_cis_effects);
//   // MatrixXd randn_theta = rstdnorm_mat2(b,p);
//   // VectorXd randn_cis   = rstdnorm_mat2(total_cis_effects,1).col(0);
// }
//
// // // [[Rcpp::export]]
// // SEXP funx() {
// //   /* creating a pointer to a vector<int> */
// //   // std::vector<int> *v = new std::vector<int>;
// //   // v->push_back( 1 ) ;
// //   // v->push_back( 2 ) ;
// //   std::vector<int> v;
// //   v.push_back( 1 ) ;
// //   v.push_back( 2 ) ;
// //
// //   /* wrap the pointer as an external pointer */
// //   /* this automatically protected the external pointer from R garbage
// //   collection until p goes out of scope. */
// //   Rcpp::XPtr< std::vector<int> > p(&v, true) ;
// //
// //   /* return it back to R, since p goes out of scope after the return
// //   the external pointer is no more protected by p, but it gets
// //   protected by being on the R side */
// //   return( p ) ;
// // }
// //
// // // [[Rcpp::export]]
// // VectorXd funx2(SEXP p1){
// //   Rcpp::XPtr< std::vector<int> > p(p1) ;
// //   std::vector<int> v = *p;
// //   VectorXd res(v.size());
// //   for(int i = 0; i < res.size(); i++){
// //     res[i] = v[i];
// //   }
// //   return(res);
// // }
//
// // [[Rcpp::export]]
// SEXP to_vector(List x){
//   std::vector<MSpMat> *y = new std::vector<MSpMat>;
//   for(int i = 0; i < x.size(); i++){
//     List xi = as<List>(x[i]);
//     // MSpMat chol_Sigma = as<MSpMat>(xi["chol_Sigma"]);
//     // y->push_back(&chol_Sigma);
//     y->push_back(as<MSpMat>(xi["chol_Sigma"]));
//   }
//   Rcpp::XPtr< std::vector<MSpMat> > p(y, true) ;
//
//   /* return it back to R, since p goes out of scope after the return
//   the external pointer is no more protected by p, but it gets
//   protected by being on the R side */
//   return( p ) ;
// }
//
//
// double get_score(const MSpMat &chol_R, const VectorXd y){
//   VectorXd score = chol_R * y;
//   return(score.dot(score));
// }
//
// struct tot_prec_scores_worker2 : public RcppParallel::Worker {
//   const MatrixXd Y;
//   const std::vector<MSpMat> chol_R_list;
//   VectorXi h2s_index;
//   VectorXd &scores;
//
//   tot_prec_scores_worker2(const MatrixXd Y,                                // nxp
//                          const std::vector<MSpMat> &chol_R_list,     // std::vector of nxn SpMat upper-triangule
//                          VectorXi h2s_index,                        // px1, 1-based index
//                          VectorXd &scores                           // px1
//   ):
//     Y(Y), chol_R_list(chol_R_list), h2s_index(h2s_index),
//     scores(scores) {}
//
//   void operator()(std::size_t begin, std::size_t end) {
//     for(std::size_t j = begin; j < end; j++){
//       int h2_index = h2s_index[j] - 1;
//       // MSpMat *chol_Rp = chol_R_list[h2_index];
//       MSpMat chol_R = chol_R_list[h2_index];
//       // const MSpMat chol_R = chol_R_list->operator[](h2_index);
//       // MSpMat * chol_Rp = &chol_R;
//       // chol_R += h2_index;
//         // &chol_R_list->operator[](h2_index);
//       VectorXd score = chol_R.transpose().triangularView<Lower>().solve(Y.col(j));
//       // VectorXd score = chol_R * Y.col(j);
//       scores[j] = score.dot(score);
//       // scores[j] = get_score(chol_R,Y.col(j));
//     }
//   }
// };
//
// // calculates normal scores for:
// // Y_resid %*% Sigma_inv %*% t(Y_resid)
// // assumes complete data Y
// // [[Rcpp::export()]]
// VectorXd tot_prec_scores2(
//     Map<MatrixXd> Y,              // nxp
//     SEXP Sigma_Choleskys_ptr,   // List of nxn dgCMatrix upper-triangule
//     VectorXi h2s_index,           // px1
//     int grainSize)
// {
//
//
//   int p = Y.cols();
//   VectorXd scores(p);
//
//   Rcpp::XPtr< std::vector<MSpMat> > Sigma_Choleskys2(Sigma_Choleskys_ptr) ;
//   std::vector<MSpMat> Sigma_Choleskys = *Sigma_Choleskys2;
//   tot_prec_scores_worker2 sampler(Y,Sigma_Choleskys,h2s_index,scores);
//
//
//   RcppParallel::parallelFor(0,p,sampler,grainSize);
//   return scores;
// }
//
//
// double get_score3(const MSpMat &chol_R, const VectorXd y){
//   VectorXd score = chol_R * y;
//   return(score.dot(score));
// }
//
// struct tot_prec_scores_worker3 : public RcppParallel::Worker {
//   const MatrixXd Y;
//   const std::vector<MSpMat> chol_R_list;
//   VectorXi h2s_index;
//   VectorXd &scores;
//
//   tot_prec_scores_worker3(const MatrixXd &Y,                                // nxp
//                           const std::vector<MSpMat> &chol_R_list,     // std::vector of nxn SpMat upper-triangule
//                           VectorXi h2s_index,                        // px1, 1-based index
//                           VectorXd &scores                           // px1
//   ):
//     Y(Y), chol_R_list(chol_R_list), h2s_index(h2s_index),
//     scores(scores) {}
//
//   void operator()(std::size_t begin, std::size_t end) {
//     for(std::size_t j = begin; j < end; j++){
//       int h2_index = h2s_index[j] - 1;
//       // MSpMat *chol_Rp = chol_R_list[h2_index];
//       MSpMat chol_R = chol_R_list[h2_index];
//       // const MSpMat chol_R = chol_R_list->operator[](h2_index);
//       // MSpMat * chol_Rp = &chol_R;
//       // chol_R += h2_index;
//       // &chol_R_list->operator[](h2_index);
//       VectorXd score = chol_R.transpose().triangularView<Lower>().solve(Y.col(j));
//       // VectorXd score = chol_R * Y.col(j);
//       scores[j] = score.dot(score);
//       // scores[j] = get_score3(chol_R,Y.col(j));
//     }
//   }
// };
//
// // [[Rcpp::export()]]
// VectorXd tot_prec_scores3(
//     Map<MatrixXd> Y,              // nxp
//     Rcpp::List Sigma_Choleskys,   // List of nxn dgCMatrix upper-triangule
//     VectorXi h2s_index,           // px1
//     int grainSize)
// {
//
//
//   int p = Y.cols();
//   VectorXd scores(p);
//
//   std::vector<MSpMat> chol_R_list;
//   for(int i = 0; i < h2s_index.maxCoeff(); i++){
//     Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);
//     chol_R_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
//   }
//   tot_prec_scores_worker3 sampler(Y,chol_R_list,h2s_index,scores);
//   RcppParallel::parallelFor(0,p,sampler,grainSize);
//   return scores;
// }
//
// // [[Rcpp::export()]]
// MatrixXd rstdnorm_mat2(int n,int p) {  // returns nxp matrix
//   VectorXd X_vec = as<VectorXd>(rnorm(n*p));
//   Map<MatrixXd> X_mat(X_vec.data(),n,p);
//   return(X_mat);
// }
//
//
// // [[Rcpp::export()]]
// VectorXd sample_MME_single_diagK2(  // returns b x 1 vector
//     const VectorXd& y,           // nx1
//     const MatrixXd& X,           // nxb
//     const VectorXd& prior_mean,  // bx1
//     const VectorXd& prior_prec,  // bx1
//     const MSpMat&   chol_Ri,        // nxn upper triangular Cholesky decomposition of R. Format: dgCMatrix
//     double tot_Eta_prec,
//     const VectorXd& randn_theta, // bx1
//     const VectorXd& randn_e      // 0x1 or nx1. 0x1 if b<n
// ){
//   if(randn_e.size() == 0){
//     MatrixXd RinvSqX = chol_Ri * (X * sqrt(tot_Eta_prec)); //
//     VectorXd XtRinvy = RinvSqX.transpose() * (chol_Ri * y * sqrt(tot_Eta_prec));
//     VectorXd XtRinvy_std_mu = XtRinvy + prior_prec.asDiagonal()*prior_mean;
//     MatrixXd C = RinvSqX.transpose() * RinvSqX;
//     C.diagonal() += prior_prec;
//     LLT<MatrixXd> C_llt;
//     C_llt.compute(C);
//     MatrixXd chol_C = C_llt.matrixU();
//
//     VectorXd b = chol_C.transpose().triangularView<Lower>().solve(XtRinvy_std_mu);
//     b += randn_theta;
//     b = chol_C.triangularView<Upper>().solve(b);
//     return(b);
//   } else {
//     // should check that randn_e.size() == n
//     VectorXd theta_star = randn_theta.array() / prior_prec.cwiseSqrt().array();
//     theta_star += prior_mean;
//     VectorXd e_star = chol_Ri.transpose().triangularView<Lower>().solve(randn_e / sqrt(tot_Eta_prec));
//     MatrixXd X_theta_star = X * theta_star;
//     VectorXd y_resid = y - X_theta_star - e_star;
//
//     MatrixXd RinvSqX = chol_Ri * (X * sqrt(tot_Eta_prec));
//     VectorXd XtRinvy = RinvSqX.transpose() * chol_Ri * (y_resid * sqrt(tot_Eta_prec));
//
//     VectorXd theta_tilda;
//     MatrixXd VAi = X * prior_prec.cwiseInverse().asDiagonal();
//     MatrixXd inner = VAi*X.transpose() + (chol_Ri.transpose() * chol_Ri) / tot_Eta_prec;  // wrong!
//     VectorXd VAiXtURinvy = VAi * XtRinvy;
//     VectorXd outerXtURinvy = VAi.transpose() * inner.ldlt().solve(VAiXtURinvy);
//     theta_tilda = XtRinvy.array() / prior_prec.array();
//     theta_tilda -= outerXtURinvy;
//     VectorXd theta = theta_star + theta_tilda;
//
//     return(theta);
//   }
// }
//
//
// // RcppParallel struct for sampling coefficients from a set of parallel regression problems:
// // Y = X %*% B + E, where
// // b[i,j] ~ N(prior_mean[i,j],1/prior_prec[i,j])
// // E[,j] ~ N(0,t(chol_R) %*% chol(R) / tot_Eta_prec[j])
// // where chol_R is selected from a list based on h2s_index[j]
// // Y is complete, so everythign has the same dimensions
// struct sample_MME_single_diagK_worker2 : public RcppParallel::Worker {
//   const MatrixXd Y;
//   const MatrixXd X;
//   const MatrixXd prior_mean;
//   const MatrixXd prior_prec;
//   const MatrixXd randn_theta;
//   const MatrixXd randn_e;
//   const std::vector<MSpMat> chol_R_list;
//   const VectorXi h2s_index;
//   const VectorXd tot_Eta_prec;
//   MatrixXd &coefs;
//
//   sample_MME_single_diagK_worker2(
//     const MatrixXd &Y,           // nxp
//     const MatrixXd &X,           // nxb
//     const MatrixXd &prior_mean,  // bxp
//     const MatrixXd &prior_prec,  // bxp
//     const std::vector<MSpMat> &chol_R_list, // each element a MSpMat nxn upper-triangular cholesky decomposition of R
//     const VectorXi &h2s_index,   // px1, 1-based index
//     const VectorXd &tot_Eta_prec,// px1
//     const MatrixXd &randn_theta, // bxp
//     const MatrixXd &randn_e,     // 0xp or nxp. The former if b<n, the latter if b >= n
//     MatrixXd &coefs       // bxp
//   ):
//     Y(Y), X(X), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
//     chol_R_list(chol_R_list), h2s_index(h2s_index), tot_Eta_prec(tot_Eta_prec), coefs(coefs) {}
//
//   void operator()(std::size_t begin, std::size_t end) {
//     for(std::size_t j = begin; j < end; j++){
//       int h2_index = h2s_index[j] - 1;
//       // const MSpMat chol_R = chol_R_list->operator[](h2_index);
//       MSpMat chol_R = chol_R_list[h2_index];
//       // chol_R *= 1/sqrt(tot_Eta_prec[j]);
//       coefs.col(j) = sample_MME_single_diagK2(Y.col(j), X, prior_mean.col(j), prior_prec.col(j), chol_R, tot_Eta_prec[j],randn_theta.col(j),randn_e.col(j));
//     }
//   }
// };
//
//
// // -------------------------------------------------- //
// // -- Versions of the independent residuals regression model --- //
// // -------------------------------------------------- //
//
//
// // basic version - all regression have same dimensions
// // [[Rcpp::export()]]
// MatrixXd sample_MME_fixedEffects_c2(  // returns bxp matrix
//     Map<MatrixXd> Y,              // nxp
//     Map<MatrixXd> X,              // nxb
//     SEXP Sigma_Choleskys_ptr,   // List of nxn dgCMatrix upper-triangule
//     VectorXi h2s_index,           // px1 index of Cholesky matrix for each column
//     Map<VectorXd> tot_Eta_prec,   // px1
//     Map<MatrixXd> prior_mean,     // bxp
//     Map<MatrixXd> prior_prec,     // bxp
//     int grainSize) {
//
//   int b = X.cols();
//   int p = Y.cols();
//   int n = Y.rows();
//
//   MatrixXd randn_theta = rstdnorm_mat2(b,p);
//   MatrixXd randn_e;
//   if(b < n) {
//     randn_e = rstdnorm_mat2(0,p);
//   } else{
//     randn_e = rstdnorm_mat2(n,p);
//   }
//
//   Rcpp::XPtr< std::vector<MSpMat> > Sigma_Choleskys2(Sigma_Choleskys_ptr) ;
//   std::vector<MSpMat> Sigma_Choleskys = *Sigma_Choleskys2;
//
//   MatrixXd coefs(b,p);
//
//   sample_MME_single_diagK_worker2 sampler(Y,X,prior_mean,prior_prec,Sigma_Choleskys,h2s_index,tot_Eta_prec,randn_theta,randn_e, coefs);
//   RcppParallel::parallelFor(0,p,sampler,grainSize);
//   return(coefs);
// }
//
// // basic version - all regression have same dimensions
// // [[Rcpp::export()]]
// MatrixXd sample_MME_fixedEffects_c3(  // returns bxp matrix
//     Map<MatrixXd> Y,              // nxp
//     Map<MatrixXd> X,              // nxb
//     Rcpp::List Sigma_Choleskys,   // list of Cholesky decompositions of residuals. nxn upper triangular, dgCMatrix
//     VectorXi h2s_index,           // px1 index of Cholesky matrix for each column
//     Map<VectorXd> tot_Eta_prec,   // px1
//     Map<MatrixXd> prior_mean,     // bxp
//     Map<MatrixXd> prior_prec,     // bxp
//     int grainSize) {
//
//   int b = X.cols();
//   int p = Y.cols();
//   int n = Y.rows();
//
//   MatrixXd randn_theta = rstdnorm_mat2(b,p);
//   MatrixXd randn_e;
//   if(b < n) {
//     randn_e = rstdnorm_mat2(0,p);
//   } else{
//     randn_e = rstdnorm_mat2(n,p);
//   }
//
//   std::vector<MSpMat> chol_R_list;
//   for(int i = 0; i < h2s_index.maxCoeff(); i++){
//     Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);
//     chol_R_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
//   }
//
//   MatrixXd coefs(b,p);
//
//   sample_MME_single_diagK_worker2 sampler(Y,X,prior_mean,prior_prec,chol_R_list,h2s_index,tot_Eta_prec,randn_theta,randn_e, coefs);
//   RcppParallel::parallelFor(0,p,sampler,grainSize);
//   return(coefs);
// }



// // [[Rcpp::export()]]
// MatrixXd diagU_solve(MSpMat chol_R, MatrixXd Y){
//   return(chol_R.triangularView<Upper>().solve(Y));
// }
// // [[Rcpp::export()]]
// MatrixXd diagU_multiply(MSpMat chol_R, MatrixXd Y){
//   return(chol_R * Y);
// }
//
//
// // [[Rcpp::export()]]
// VectorXd sample_MME_single_diagK2(  // returns b x 1 vector
//     VectorXd y,           // nx1
//     MatrixXd X,           // nxb
//     VectorXd prior_mean,  // bx1
//     VectorXd prior_prec,  // bx1
//     SpMat chol_Sigma_inv,        // nxn upper triangular Cholesky decomposition of R. Format: dgCMatrix
//     VectorXd randn_theta, // bx1
//     VectorXd randn_e      // 0x1 or nx1. 0x1 if b<n
// ){
//     MatrixXd Sigma_invSq_X = chol_Sigma_inv * X;
//     VectorXd Xt_Sigma_inv_Y = Sigma_invSq_X.transpose() * (chol_Sigma_inv * y);
//     VectorXd Xt_Sigma_inv_Y_std_mu = Xt_Sigma_inv_Y + prior_prec.asDiagonal()*prior_mean;
//     MatrixXd C = Sigma_invSq_X.transpose() * Sigma_invSq_X;
//     C.diagonal() += prior_prec;
//     LLT<MatrixXd> C_llt;
//     C_llt.compute(C);
//     MatrixXd chol_C = C_llt.matrixU();
//
//     VectorXd b = chol_C.transpose().triangularView<Lower>().solve(Xt_Sigma_inv_Y_std_mu);
//     b += randn_theta;
//     b = chol_C.triangularView<Upper>().solve(b);
//     return(b);
//
// }
//
// // [[Rcpp::export()]]
// VectorXd sample_MME_single_diagK3(  // returns b x 1 vector
//     VectorXd y,           // nx1
//     MatrixXd X,           // nxb
//     VectorXd prior_mean,  // bx1
//     VectorXd prior_prec,  // bx1
//     SpMat chol_R,        // nxn upper triangular Cholesky decomposition of R. Format: dgCMatrix
//     VectorXd randn_theta, // bx1
//     VectorXd randn_e      // 0x1 or nx1. 0x1 if b<n
// ){
//     MatrixXd RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X);
//     VectorXd XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y);
//     VectorXd XtRinvy_std_mu = XtRinvy + prior_prec.asDiagonal()*prior_mean;
//
//     MatrixXd C = RinvSqX.transpose() * RinvSqX;
//     C.diagonal() += prior_prec;
//     LLT<MatrixXd> C_llt;
//     C_llt.compute(C);
//     MatrixXd chol_C = C_llt.matrixU();
//
//     VectorXd b = chol_C.transpose().triangularView<Lower>().solve(XtRinvy_std_mu);
//     b += randn_theta;
//     b = chol_C.triangularView<Upper>().solve(b);
//     return(b);
// }
//
//
//
// //
// //
// // [[Rcpp::export()]]
// MatrixXd rstdnorm_mat2(int n,int p) {  // returns nxp matrix
//   VectorXd X_vec = as<VectorXd>(rnorm(n*p));
//   Map<MatrixXd> X_mat(X_vec.data(),n,p);
//   return(X_mat);
// }
// //
// // // [[Rcpp::export()]]
// // VectorXd sample_MME_single_diagR2(
// //     VectorXd y,           // nx1
// //     SpMat ZQt,              // nxr dgCMatrix
// //     SpMat Q,              // nxr dgCMatrix
// //     SpMat chol_S,         // rxr dgCMatrix upper triangular: chol(ZtRinvZ + diag(Kinv))
// //     double pe,            // double
// //     VectorXd randn_theta // rx1
// // ){
// //   VectorXd b = ZQt * y * pe;
// //   b = chol_S.transpose().triangularView<Lower>().solve(b);
// //   b += randn_theta;
// //   b = Q * chol_S.triangularView<Upper>().solve(b);
// //   return(b);
// // }
// //
// // // [[Rcpp::export()]]
// // VectorXd sample_MME_single_diagK(  // returns b x 1 vector
// //     VectorXd y,           // nx1
// //     MatrixXd X,           // nxb
// //     VectorXd prior_mean,  // bx1
// //     VectorXd prior_prec,  // bx1
// //     SpMat chol_R,        // nxn upper triangular Cholesky decomposition of R. Format: dgCMatrix
// //     VectorXd randn_theta, // bx1
// //     VectorXd randn_e      // 0x1 or nx1. 0x1 if b<n
// // ){
// //   if(randn_e.size() == 0){
// //     MatrixXd RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X);
// //     VectorXd XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y);
// //     MatrixXd C = RinvSqX.transpose() * RinvSqX;
// //     C.diagonal() += prior_prec;
// //     LLT<MatrixXd> C_llt;
// //     C_llt.compute(C);
// //     MatrixXd chol_C = C_llt.matrixU();
// //
// //     VectorXd b = chol_C.transpose().triangularView<Lower>().solve(XtRinvy);
// //     b += randn_theta;
// //     b = chol_C.triangularView<Upper>().solve(b);
// //     return(b);
// //   } else {
// //     // should check that randn_e.size() == n
// //     VectorXd theta_star = randn_theta.array() / prior_prec.cwiseSqrt().array();
// //     theta_star += prior_mean;
// //     VectorXd e_star = chol_R * randn_e;
// //     MatrixXd X_theta_star = X * theta_star;
// //     VectorXd y_resid = y - X_theta_star - e_star;
// //
// //     MatrixXd RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X);
// //     VectorXd XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y_resid);
// //
// //     VectorXd theta_tilda;
// //     MatrixXd VAi = X * prior_prec.cwiseInverse().asDiagonal();
// //     MatrixXd inner = VAi*X.transpose() + chol_R.transpose() * chol_R;
// //     VectorXd VAiXtURinvy = VAi * XtRinvy;
// //     VectorXd outerXtURinvy = VAi.transpose() * inner.ldlt().solve(VAiXtURinvy);
// //     theta_tilda = XtRinvy.array() / prior_prec.array();
// //     theta_tilda -= outerXtURinvy;
// //     VectorXd theta = theta_star + theta_tilda;
// //
// //     return(theta);
// //   }
// // }
// //
// //
// // // RcppParallel struct for sampling coefficients from a set of parallel regression problems:
// // // Y = X %*% B + E, where
// // // b[i,j] ~ N(prior_mean[i,j],1/prior_prec[i,j])
// // // E[,j] ~ N(0,t(chol_R) %*% chol(R) / tot_Eta_prec[j])
// // // where chol_R is selected from a list based on h2s_index[j]
// // // Y is complete, so everythign has the same dimensions
// // struct sample_MME_single_diagK_worker : public RcppParallel::Worker {
// //   MatrixXd Y;
// //   MatrixXd X;
// //   MatrixXd prior_mean, prior_prec, randn_theta, randn_e;
// //   const std::vector<MSpMat> chol_R_list;
// //   VectorXi h2s_index;
// //   VectorXd tot_Eta_prec;
// //   MatrixXd &coefs;
// //
// //   sample_MME_single_diagK_worker(
// //     MatrixXd Y,           // nxp
// //     MatrixXd X,           // nxb
// //     MatrixXd prior_mean,  // bxp
// //     MatrixXd prior_prec,  // bxp
// //     const std::vector<MSpMat> chol_R_list, // each element a MSpMat nxn upper-triangular cholesky decomposition of R
// //     VectorXi h2s_index,   // px1, 1-based index
// //     VectorXd tot_Eta_prec,// px1
// //     MatrixXd randn_theta, // bxp
// //     MatrixXd randn_e,     // 0xp or nxp. The former if b<n, the latter if b >= n
// //     MatrixXd &coefs       // bxp
// //   ):
// //     Y(Y), X(X), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
// //     chol_R_list(chol_R_list), h2s_index(h2s_index), tot_Eta_prec(tot_Eta_prec), coefs(coefs) {}
// //
// //   void operator()(std::size_t begin, std::size_t end) {
// //     for(std::size_t j = begin; j < end; j++){
// //       int h2_index = h2s_index[j] - 1;
// //       SpMat chol_R = chol_R_list[h2_index];
// //       chol_R *= 1/sqrt(tot_Eta_prec[j]);
// //       coefs.col(j) = sample_MME_single_diagK(Y.col(j), X, prior_mean.col(j), prior_prec.col(j), chol_R, randn_theta.col(j),randn_e.col(j));
// //     }
// //   }
// // };
// //
// //
// // // -------------------------------------------------- //
// // // -- Versions of the independent residuals regression model --- //
// // // -------------------------------------------------- //
// //
// //
// // // basic version - all regression have same dimensions
// // // [[Rcpp::export()]]
// // MatrixXd sample_MME_fixedEffects_c(  // returns bxp matrix
// //     Map<MatrixXd> Y,              // nxp
// //     Map<MatrixXd> X,              // nxb
// //     Rcpp::List Sigma_Choleskys,   // list of Cholesky decompositions of residuals. nxn upper triangular, dgCMatrix
// //     VectorXi h2s_index,           // px1 index of Cholesky matrix for each column
// //     Map<VectorXd> tot_Eta_prec,   // px1
// //     Map<MatrixXd> prior_mean,     // bxp
// //     Map<MatrixXd> prior_prec,     // bxp
// //     int grainSize) {
// //
// //   int b = X.cols();
// //   int p = Y.cols();
// //   int n = Y.rows();
// //
// //   MatrixXd randn_theta = rstdnorm_mat2(b,p);
// //   MatrixXd randn_e;
// //   if(b < n) {
// //     randn_e = rstdnorm_mat2(0,p);
// //   } else{
// //     randn_e = rstdnorm_mat2(n,p);
// //   }
// //
// //   std::vector<MSpMat> chol_R_list;
// //   for(int i = 0; i < h2s_index.maxCoeff(); i++){
// //     Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);
// //     chol_R_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
// //   }
// //
// //   MatrixXd coefs(b,p);
// //
// //   sample_MME_single_diagK_worker sampler(Y,X,prior_mean,prior_prec,chol_R_list,h2s_index,tot_Eta_prec,randn_theta,randn_e, coefs);
// //   RcppParallel::parallelFor(0,p,sampler,grainSize);
// //   return(coefs);
// // }
