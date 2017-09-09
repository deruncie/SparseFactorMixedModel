// #include <math.h>
// #include <iostream>
// #include "BSFG_types.h"
//
// // [[Rcpp::depends(RcppEigen)]]
// using namespace Eigen;
// using namespace RcppParallel;
//
//
// // -------------------------- //
// // -------------------------- //
// // -------------------------- //
//
// // functions to speed up sparse multiplication and conversion to dense matrices
//
// // [[Rcpp::export()]]
// MatrixXd SxD(MSpMat X, Map<MatrixXd> Y){
//   return(X * Y);
// }
// // [[Rcpp::export()]]
// MatrixXd SxS(MSpMat X, MSpMat Y){
//   return(X * Y);
// }
//
//
// /*** R
// `%**%` = function(X1,X2){
//   if(inherits(X1,'dgCMatrix') && inherits(X2,'matrix')) return(SxD(X1,X2))
//   if(inherits(X1,'dgCMatrix') && inherits(X2,'dgCMatrix')) return(SxS(X1,X2))
//   if(inherits(X1,'matrix') & inherits(X2,'matrix')) return(X1 %*% X2)
//   return(as.matrix(X1 %*% X2))
// }
// */
//
//
//
// // -------------------------- //
// // ----- Helper functions --- //
// // -------------------------- //
// // [[Rcpp::export()]]
// MatrixXd uncorrelated_prec_mat(
//     VectorXd h2,
//     VectorXd tot_prec,
//     VectorXd s
// ){
//   int n = s.size();
//   int p = h2.size();
//
//   MatrixXd Prec(n,p);
//
//   struct perColumn : public RcppParallel::Worker {
//     VectorXd h2;
//     VectorXd tot_prec;
//     VectorXd s;
//     MatrixXd &Prec;
//
//     perColumn(VectorXd h2, VectorXd tot_prec, VectorXd s,
//               MatrixXd &Prec):
//       h2(h2), tot_prec(tot_prec),s(s), Prec(Prec){}
//
//     void operator()(std::size_t begin, std::size_t end) {
//       for(std::size_t j = begin; j < end; j++){
//         Prec.col(j) = tot_prec(j) * (h2(j) * s.array() + (1.0-h2(j))).inverse();
//       }
//     }
//   };
//
//   perColumn job(h2,tot_prec,s,Prec);
//   RcppParallel::parallelFor(0,p,job,1);
//   return(Prec);
// }
//
//
// // -------------------------- //
// // ----- sample coeffients ----- //
// // -------------------------- //
//
// // basic functions for sampling from regression models with a-priori independent coefficients and residuals
//
// // draws a sample of the vector b from the model:
// // y = X %*% beta + e
// // with beta[i] ~ N(prior_mean[i],1/prior_prec[i]), i=1:b
// // with e[j] ~ N(0,1/resid_prec[i]), i=1:n
// // Uses sampling method from MCMCglmm, which requires draws b+n draws from N(0,1), which are passed as randn_theta and randn_e
// // If b >= n, inverts the C matrix using Binomial Inverse Theorem
// VectorXd sample_coefs_uncorrelated(
//     VectorXd y,
//     MatrixXd X,
//     VectorXd prior_mean,
//     VectorXd prior_prec,
//     ArrayXd  resid_prec,
//     VectorXd randn_theta,
//     VectorXd randn_e,
//     int b,
//     int n
// ) {
//   VectorXd R_sq_diag = resid_prec.sqrt().inverse();
//   VectorXd theta_star = randn_theta.array()/prior_prec.array().sqrt();
//   theta_star += prior_mean;
//   VectorXd e_star = randn_e.array() * R_sq_diag.array();
//   MatrixXd UtW_theta_star = X * theta_star;
//   VectorXd eta_resid = y - UtW_theta_star - e_star;
//   MatrixXd RinvSqUtW = R_sq_diag.cwiseInverse().asDiagonal() * X;
//   VectorXd eta_std = eta_resid.array()/R_sq_diag.array();
//   VectorXd WtURinvy = RinvSqUtW.transpose() * eta_std;
//
//   VectorXd theta_tilda;
//   if(b < n) {
//     MatrixXd C = RinvSqUtW.transpose() * RinvSqUtW;
//     C.diagonal() += prior_prec;
//     theta_tilda = C.llt().solve(WtURinvy);
//   } else{
//     MatrixXd VAi = X * prior_prec.cwiseInverse().asDiagonal();
//     MatrixXd inner = VAi*X.transpose();
//     inner.diagonal() += R_sq_diag.cwiseProduct(R_sq_diag);
//     VectorXd VAiWtURinvy = VAi * WtURinvy;
//     VectorXd outerWtURinvy = VAi.transpose() * inner.ldlt().solve(VAiWtURinvy);
//     theta_tilda = WtURinvy.array() / prior_prec.array();
//     theta_tilda -= outerWtURinvy;
//   }
//
//   VectorXd coefs = theta_tilda + theta_star;
//   return coefs;
// }
//
// // RcppParallel struct for sampling coefficients from a set of parallel regression problems:
// // Y = X %*% B + E, where
// // b[i,j] ~ N(prior_mean[i,j],1/prior_prec[i,j])
// // E[i,j] ~ N(0,1/resid_prec[i,j])
// // Where Y is complete, so everythign has the same dimensions
// struct sample_coefs_uncorrelated_worker : public Worker {
//   MatrixXd Y, X;
//   MatrixXd resid_prec, prior_mean, prior_prec, randn_theta, randn_e;
//   int b,n;
//   MatrixXd &coefs;
//
//   sample_coefs_uncorrelated_worker( MatrixXd Y, MatrixXd X,MatrixXd resid_prec,MatrixXd prior_mean, MatrixXd prior_prec,
//                                     MatrixXd randn_theta, MatrixXd randn_e,
//                                     int b, int n,
//                                     MatrixXd &coefs) :
//     Y(Y), X(X), resid_prec(resid_prec), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
//     b(b), n(n),
//     coefs(coefs) {}
//
//   void operator()(std::size_t begin, std::size_t end) {
//     for(std::size_t j = begin; j < end; j++){
//       coefs.col(j) = sample_coefs_uncorrelated(Y.col(j), X, prior_mean.col(j), prior_prec.col(j), resid_prec.col(j), randn_theta.col(j),randn_e.col(j),b,n);
//     }
//   }
// };
//
//
// // -------------------------------------------------- //
// // -- Versions of the independent regression model --- //
// // -------------------------------------------------- //
//
// // basic version - all regression have same dimensions
// // [[Rcpp::export()]]
// MatrixXd sample_coefs_parallel_sparse_c_Eigen(
//     Map<MatrixXd> Y,
//     Map<MatrixXd> X,
//     Map<MatrixXd> resid_prec,
//     Map<MatrixXd> prior_mean,
//     Map<MatrixXd> prior_prec,
//     Map<MatrixXd> randn_theta,
//     Map<MatrixXd> randn_e,
//     int grainSize){
//
//   // Sample regression coefficients
//   // columns of matrices are independent
//   // each column conditional posterior is a MVN due to conjugacy
//
//   int p = Y.cols();
//   int b = X.cols();
//   int n = X.rows();
//
//   MatrixXd coefs(b,p);
//
//   sample_coefs_uncorrelated_worker sampler(Y, X, resid_prec, prior_mean, prior_prec, randn_theta,randn_e,b, n, coefs);
//   RcppParallel::parallelFor(0,p,sampler,grainSize);
//
//   return(coefs);
// }
//
//
// // missing data version - sets of columns of Eta have different patterns of missing data
// // each set can be sampled at once using sample_coefs_uncorrelated_worker
// // [[Rcpp::export()]]
// MatrixXd sample_coefs_parallel_sparse_missing_c_Eigen(
//     Map<MatrixXd> Eta,
//     Map<MatrixXd> W,
//     Map<VectorXd> h2,
//     Map<VectorXd> tot_Eta_prec,
//     Map<MatrixXd> prior_mean,
//     Map<MatrixXd> prior_prec,
//     Map<MatrixXd> randn_theta,
//     Map<VectorXd> randn_e,
//     Rcpp::List invert_aI_bZKZ,
//     int grainSize){
//
//   // Sample regression coefficients
//   // columns of matrices are independent
//   // each column conditional posterior is a MVN due to conjugacy
//
//   int p = Eta.cols();
//   int b = W.cols();
//   MatrixXd coefs(b,p);
//
//   // sample sets of columns with same pattern of missing data as a block
//   int randn_e_index = 0;
//   for(int col_set = 0; col_set < invert_aI_bZKZ.length(); col_set++){
//     Rcpp::List invert_aI_bZKZ_i = Rcpp::as<Rcpp::List>(invert_aI_bZKZ[col_set]);
//     SpMat Ut = as<MSpMat>(invert_aI_bZKZ_i["Ut"]);
//     VectorXd s = as<Map<VectorXd>>(invert_aI_bZKZ_i["s"]);
//     VectorXi Y_obs = as<Map<VectorXi>>(invert_aI_bZKZ_i["Y_obs"]);
//     VectorXi Y_cols = as<Map<VectorXi>>(invert_aI_bZKZ_i["Y_cols"]);
//
//     int n_obs = Y_obs.size();
//     int n_cols = Y_cols.size();
//
//     // pull out rows of X corresponding to observed individuals, pre-multiply by Ut
//     MatrixXd X_set(n_obs,b);
//     for(int i = 0; i < n_obs; i++) {
//       X_set.row(i) = W.row(Y_obs[i]-1);
//     }
//     X_set = Ut * X_set;
//
//     // pull out elements of Eta corresponding to observed individuals, pre-multiply by Ut
//     MatrixXd Y_set(n_obs,n_cols);
//     for(int i = 0; i < n_obs; i++){
//       int row_i = Y_obs[i]-1;
//       for(int j = 0; j < n_cols; j++){
//         int col_j = Y_cols[j]-1;
//         Y_set.coeffRef(i,j) = Eta.coeffRef(row_i,col_j);
//       }
//     }
//     Y_set = Ut * Y_set;
//
//     // pull out h2s, tot_precs for col_set, calculate resid_prec
//     VectorXd h2_set(n_cols);
//     VectorXd tot_Eta_prec_set(n_cols);
//     for(int j = 0; j < n_cols; j++){
//       int col_j = Y_cols[j]-1;
//       h2_set[j] = h2[col_j];
//       tot_Eta_prec_set[j] = tot_Eta_prec[col_j];
//     }
//     MatrixXd resid_prec_set = uncorrelated_prec_mat(h2_set,tot_Eta_prec_set,s);
//
//     // pull out prior_mean, prior_prec, randn_theta, randn_e for col_set
//     MatrixXd prior_mean_set(b,n_cols);
//     MatrixXd prior_prec_set(b,n_cols);
//     MatrixXd randn_theta_set(b,n_cols);
//     MatrixXd randn_e_set(n_obs,n_cols);
//     for(int j = 0; j < n_cols; j++){
//       int col_j = Y_cols[j]-1;
//       prior_mean_set.col(j) = prior_mean.col(col_j);
//       prior_prec_set.col(j) = prior_prec.col(col_j);
//       randn_theta_set.col(j) = randn_theta.col(col_j);
//       randn_e_set.col(j) = randn_e.segment(randn_e_index,n_obs);
//       randn_e_index += n_obs;
//     }
//
//     MatrixXd coefs_set(b,n_cols);
//
//     // sample coefs_set as a block
//     sample_coefs_uncorrelated_worker sampler(Y_set, X_set, resid_prec_set, prior_mean_set, prior_prec_set, randn_theta_set, randn_e_set, b, n_obs, coefs_set);
//     RcppParallel::parallelFor(0,n_cols,sampler,grainSize);
//
//     // copy new coefs into full matrix
//     for(int j = 0; j < n_cols; j++){
//       coefs.col(Y_cols[j]-1) = coefs_set.col(j);
//     }
//   }
//   return(coefs);
// }
//
//
//
// // [[Rcpp::export()]]
// MatrixXd sample_coefs_set_c(
//     Rcpp::List model_matrices,
//     Map<VectorXd> randn_theta_vec,
//     Map<VectorXd> randn_e_vec,
//     Map<MatrixXd> h2s,
//     Map<MatrixXd> tot_Eta_prec,
//     Map<MatrixXd> prior_mean,
//     Map<MatrixXd> prior_prec,
//     int n,
//     int grainSize){
//
//   std::vector<MatrixXd> y_list;
//   std::vector<MatrixXd> X_list;
//   std::vector<VectorXd> s_list;
//   std::vector<MatrixXd> randn_theta_list;
//   std::vector<MatrixXd> randn_e_list;
//   int randn_theta_index = 0;
//   int randn_e_index = 0;
//   int total_obs = 0;
//   for(int i = 0; i < n; i++){
//     Rcpp::List model_matrix_i = Rcpp::as<Rcpp::List>(model_matrices[i]);
//     y_list.push_back(Rcpp::as<MatrixXd>(model_matrix_i["y"]));
//     X_list.push_back(Rcpp::as<MatrixXd>(model_matrix_i["X"]));
//     s_list.push_back(Rcpp::as<VectorXd>(model_matrix_i["s"]));
//
//     int p = y_list[i].cols();
//     int n_obs = y_list[i].rows();
//     int b = X_list[i].cols();
//
//     total_obs += n_obs;
//
//     MatrixXd r_theta = randn_theta_vec.segment(randn_theta_index,b*p);
//     Map<MatrixXd> randn_theta(r_theta.data(),b,p);
//     randn_theta_list.push_back(randn_theta);
//     randn_theta_index += b*p;
//     MatrixXd r_e = randn_e_vec.segment(randn_e_index,n_obs*p);
//     Map<MatrixXd> randn_e(r_e.data(),n_obs,p);
//     randn_e_list.push_back(randn_e);
//     randn_e_index += n_obs*p;
//   }
//
//
//   int n_traits = y_list[0].cols();
//   int b = randn_theta_list[0].rows()*n_traits;
//
//   MatrixXd coefs(b,n);
//
//   struct sampleColumn : public RcppParallel::Worker {
//     std::vector<MatrixXd> y_list;
//     std::vector<MatrixXd> X_list;
//     std::vector<MatrixXd> randn_theta_list;
//     std::vector<MatrixXd> randn_e_list;
//     std::vector<VectorXd> s_list;
//     MatrixXd h2s,tot_Eta_prec;
//     MatrixXd prior_mean,prior_prec;
//     int n_traits;
//     MatrixXd &coefs;
//
//     sampleColumn(
//       std::vector<MatrixXd> y_list,
//       std::vector<MatrixXd> X_list,
//       std::vector<MatrixXd> randn_theta_list,
//       std::vector<MatrixXd> randn_e_list,
//       std::vector<VectorXd> s_list,
//       MatrixXd h2s,
//       MatrixXd tot_Eta_prec,
//       MatrixXd prior_mean,
//       MatrixXd prior_prec,
//       int n_traits,
//       MatrixXd &coefs) :
//       y_list(y_list), X_list(X_list), randn_theta_list(randn_theta_list), randn_e_list(randn_e_list),s_list(s_list),
//       h2s(h2s), tot_Eta_prec(tot_Eta_prec),prior_mean(prior_mean), prior_prec(prior_prec),n_traits(n_traits),
//       coefs(coefs)
//     {}
//
//     void operator()(std::size_t begin, std::size_t end) {
//       for(std::size_t j = begin; j < end; j++){
//         int b = randn_theta_list[j].rows();
//         int n_obs = randn_e_list[j].rows();
//         VectorXd s = s_list[j];
//         for(int t = 0; t < n_traits; t++) {
//           VectorXd resid_prec = tot_Eta_prec(j,t) * (h2s(j,t) * s.array() + (1.0-h2s(j,t))).inverse();
//           coefs.block(t*b,j,b,1) = sample_coefs_uncorrelated(y_list[j].col(t), X_list[j], prior_mean.block(t*b,j,b,1),
//                       prior_prec.block(t*b,j,b,1), resid_prec, randn_theta_list[j].col(t),randn_e_list[j].col(t),b,n_obs);
//         }
//       }
//     }
//   };
//
//
//   sampleColumn sampler(y_list, X_list,randn_theta_list,randn_e_list,s_list,h2s,tot_Eta_prec,prior_mean,prior_prec,n_traits,coefs);
//   RcppParallel::parallelFor(0,n,sampler,grainSize);
//
//   return(coefs);
// }
//
//
// // [[Rcpp::export()]]
// MatrixXd get_fitted_set_c(
//     Rcpp::List model_matrices,
//     Map<MatrixXd> coefs,
//     int grainSize){
//
//   std::vector<MatrixXd> X_list;
//   std::vector<VectorXd> position_list;
//   int total_obs = 0;
//   int n = model_matrices.size();
//   for(int i = 0; i < n; i++){
//     Rcpp::List model_matrix_i = Rcpp::as<Rcpp::List>(model_matrices[i]);
//     X_list.push_back(Rcpp::as<MatrixXd>(model_matrix_i["X"]));
//     position_list.push_back(Rcpp::as<VectorXd>(model_matrix_i["position"]));
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
//     std::vector<VectorXd> position_list;
//     std::vector<VectorXd> s_list;
//     int n_traits;
//     MatrixXd coefs;
//     MatrixXd &Y_fitted;
//
//     sampleColumn(
//       std::vector<MatrixXd> X_list,
//       std::vector<VectorXd> position_list,
//       int n_traits,
//       MatrixXd coefs,
//       MatrixXd &Y_fitted) :
//       X_list(X_list), position_list(position_list),
//       n_traits(n_traits),
//       coefs(coefs),Y_fitted(Y_fitted)
//     {}
//
//     void operator()(std::size_t begin, std::size_t end) {
//       for(std::size_t j = begin; j < end; j++){
//         int b = X_list[j].cols();
//         Map<MatrixXd> Eta_i(coefs.col(j).data(),b,n_traits);
//         MatrixXd Y_fitted_j = X_list[j] * Eta_i;
//         for(int i = 0; i < position_list[j].size(); i++) {
//           Y_fitted.row(position_list[j][i]-1) = Y_fitted_j.row(i);
//         }
//       }
//     }
//   };
//
//   sampleColumn sampler(X_list,position_list,n_traits,coefs,Y_fitted);
//   RcppParallel::parallelFor(0,n,sampler,grainSize);
//
//   return(Y_fitted);
// }
//
//
// // needs to be tested!
// // [[Rcpp::export()]]
// Rcpp::List sample_cis_coefs_parallel_sparse_c_Eigen(
//     Map<MatrixXd> Y,
//     Map<MatrixXd> X,
//     Rcpp::List cis_genotypes,  // note: must be Ut*cis_X
//     Map<MatrixXd> resid_prec,
//     Map<MatrixXd> prior_mean,
//     Map<MatrixXd> prior_prec,
//     Map<MatrixXd> randn_theta,
//     Map<MatrixXd> randn_e,
//     Map<VectorXd> randn_cis,
//     Map<VectorXd> cis_effect_index,
//     int grainSize){
//
//   // Sample regression coefficients
//   // columns of matrices are independent
//   // each column conditional posterior is a MVN due to conjugacy
//
//   int p = Y.cols();
//   int b = X.cols();
//   int n = X.rows();
//
//   std::vector<MatrixXd> cis_X;
//   int length_cis = 0;
//   for(int i = 0; i < p; i++){
//     MatrixXd cXi = Rcpp::as<MatrixXd>(cis_genotypes[i]);
//     cis_X.push_back(cXi);
//     length_cis += cXi.cols();
//   }
//
//   MatrixXd coefs(b,p);
//   VectorXd cis_effects(length_cis);
//
//   struct sampleColumn : public Worker {
//     MatrixXd Y,X;
//     std::vector<MatrixXd> cis_X;
//     MatrixXd resid_prec, prior_mean, prior_prec, randn_theta, randn_e;
//     VectorXd randn_cis;
//     VectorXd cis_effect_index;
//     int b,n;
//     MatrixXd &coefs;
//     VectorXd &cis_effects;
//
//     sampleColumn(MatrixXd Y, MatrixXd X, std::vector<MatrixXd> cis_X,
//                  MatrixXd resid_prec, MatrixXd prior_mean, MatrixXd prior_prec,
//                  VectorXd cis_effect_index,
//                  MatrixXd randn_theta, MatrixXd randn_e,VectorXd randn_cis,
//                  int b, int n,
//                  MatrixXd &coefs, VectorXd &cis_effects) :
//       Y(Y),X(X),
//       cis_X(cis_X),
//       resid_prec(resid_prec), prior_mean(prior_mean), prior_prec(prior_prec),
//       randn_theta(randn_theta), randn_e(randn_e), randn_cis(randn_cis),
//       cis_effect_index(cis_effect_index),
//       b(b), n(n),
//       coefs(coefs), cis_effects(cis_effects) {}
//
//     void operator()(std::size_t begin, std::size_t end) {
//       for(std::size_t j = begin; j < end; j++){
//         int b_cis = cis_X[j].cols();
//         MatrixXd X_cisj(n,b+b_cis);
//         X_cisj << X, cis_X[j];
//
//         VectorXd prior_mean_j = VectorXd::Zero(b+b_cis);
//         prior_mean_j.head(b) = prior_mean.col(j);
//
//         VectorXd prior_prec_j = VectorXd::Constant(b+b_cis,1e-10);
//         prior_prec_j.head(b) = prior_prec.col(j);
//
//         VectorXd randn_theta_j(b+b_cis);
//         randn_theta_j.head(b) = randn_theta.col(j);
//         randn_theta_j.tail(b_cis) = randn_cis.segment(cis_effect_index[j],b_cis);
//
//         VectorXd result = sample_coefs_uncorrelated(Y.col(j), X_cisj, prior_mean_j, prior_prec_j, resid_prec.col(j), randn_theta_j,randn_e.col(j),b,n);
//         coefs.col(j) = result.head(b);
//         cis_effects.segment(cis_effect_index[j],b_cis) = result.tail(b_cis);
//       }
//     }
//   };
//
//   sampleColumn sampler(Y, X, cis_X, resid_prec, prior_mean, prior_prec, cis_effect_index, randn_theta,randn_e, randn_cis, b, n, coefs,cis_effects);
//   RcppParallel::parallelFor(0,p,sampler,grainSize);
//
//   return(Rcpp::List::create(coefs,cis_effects));
// }
//
//
// // -------------------------------------------------- //
// // -- Sample random effects --- //
// // -------------------------------------------------- //
//
// // Parallel worker for sampling random effects from model
// // Y = ZU + E
// // U[,j] ~ N(0,1/a_prec[j] * K)
// // E[,j] ~ N(0,1/e_prec[j] * I_n)
// // Where matrix Q and vectors s1 and s2 diagonalize ZZt + K
// struct sample_random_effects_worker : public Worker {
//   MatrixXd Eta;
//   SpMat Q, ZQt;
//   ArrayXd s1, s2;
//   ArrayXd a_prec, e_prec;
//   ArrayXXd randn_draws;
//   MatrixXd &effects;
//
//   sample_random_effects_worker(
//     MatrixXd Eta,SpMat Q, SpMat ZQt,
//     ArrayXd s1, ArrayXd s2,
//     ArrayXd a_prec, ArrayXd e_prec,
//     ArrayXXd randn_draws,
//     MatrixXd &effects) :
//     Eta(Eta), Q(Q), ZQt(ZQt),
//     s1(s1), s2(s2),
//     a_prec(a_prec), e_prec(e_prec),
//     randn_draws(randn_draws),
//     effects(effects) {}
//
//   void operator()(std::size_t begin, std::size_t end) {
//     for(std::size_t j = begin; j < end; j++){
//       ArrayXd b = ZQt * Eta.col(j) * e_prec(j);
//       ArrayXd d = s2*a_prec(j) + s1*e_prec(j);
//       ArrayXd mlam = b / d;
//       effects.col(j) = Q * (mlam + randn_draws.col(j) / sqrt(d)).matrix();
//     }
//   }
// };
//
//
// // samples random effects from model:
// // Y = ZU + E
// // U[,j] ~ N(0,1/a_prec[j] * K)
// // E[,j] ~ N(0,1/e_prec[j] * I_n)
// // Where matrix Q and vectors s1 and s2 diagonalize ZZt + K
// // invert_aZZt_Kinv has parameters to diagonalize a*Z*Z' + b*K for fast inversion
// // no missing observations
// // [[Rcpp::export()]]
// MatrixXd sample_randomEffects_parallel_sparse_c_Eigen (
//     Map<MatrixXd> Eta,
//     MSpMat Z,
//     Map<ArrayXd> tot_prec,
//     Map<ArrayXd> h2,
//     List invert_aZZt_Kinv,
//     Map<ArrayXXd> randn_draws,
//     int grainSize) {
//
//   ArrayXd a_prec = tot_prec / h2;
//   ArrayXd e_prec = tot_prec / (1.0-h2);
//
//   SpMat Q = as<MSpMat>(invert_aZZt_Kinv["Q"]);
//   VectorXd s1 = as<VectorXd>(invert_aZZt_Kinv["s1"]);
//   VectorXd s2 = as<VectorXd>(invert_aZZt_Kinv["s2"]);
//
//   int p = Eta.cols();
//   int r = Z.cols();
//   SpMat ZQt = (Z*Q).transpose();
//
//   MatrixXd effects(r,p);
//
//   sample_random_effects_worker sampler(Eta,Q,ZQt, s1, s2, a_prec, e_prec, randn_draws, effects);
//   RcppParallel::parallelFor(0,p,sampler,grainSize);
//
//   return(effects);
// }
//
// // samples random effects from model:
// // Y = ZU + E
// // U[,j] ~ N(0,1/a_prec[j] * K)
// // E[,j] ~ N(0,1/e_prec[j] * I_n)
// // Where matrix Q and vectors s1 and s2 diagonalize ZZt + K
// // invert_aZZt_Kinv has parameters to diagonalize a*Z*Z' + b*K for fast inversion
// // with missing observations; Samples sets of columns with same patterns of missingness
// // [[Rcpp::export()]]
// MatrixXd sample_randomEffects_parallel_sparse_missing_c_Eigen (
//     Map<MatrixXd> Eta,
//     Map<ArrayXd> tot_prec,
//     Map<ArrayXd> h2,
//     List invert_aZZt_Kinv,
//     Map<ArrayXXd> randn_draws,
//     int grainSize) {
//
//   ArrayXd a_prec = tot_prec / h2;
//   ArrayXd e_prec = tot_prec / (1.0-h2);
//
//
//   int p = Eta.cols();
//   int r = randn_draws.rows();
//   MatrixXd U(r,p);
//
//   // process sets of columns with identical patterns of missing observations
//   for(int col_set = 0; col_set < invert_aZZt_Kinv.length(); col_set++){
//     // Extract pre-calculated matrices for this set of columns
//     Rcpp::List invert_aZZt_Kinv_i = Rcpp::as<Rcpp::List>(invert_aZZt_Kinv[col_set]);
//     SpMat Q = as<MSpMat>(invert_aZZt_Kinv_i["Q"]);
//     SpMat ZQt = as<MSpMat>(invert_aZZt_Kinv_i["ZQt"]);
//     VectorXd s1 = as<VectorXd>(invert_aZZt_Kinv_i["s1"]);
//     VectorXd s2 = as<VectorXd>(invert_aZZt_Kinv_i["s2"]);
//     VectorXi Y_obs = as<VectorXi>(invert_aZZt_Kinv_i["Y_obs"]);
//     VectorXi Y_cols = as<VectorXi>(invert_aZZt_Kinv_i["Y_cols"]);
//
//     int n_obs = Y_obs.size();
//     int n_cols = Y_cols.size();
//
//     // extract matrix of observations
//     MatrixXd Eta_set(n_obs,n_cols);
//     for(int j = 0; j < n_cols; j++){
//       int col_j = Y_cols[j]-1;
//       for(int i = 0; i < n_obs; i++){
//         int row_i = Y_obs[i]-1;
//         Eta_set.coeffRef(i,j) = Eta.coeffRef(row_i,col_j);
//       }
//     }
//
//     // extract parameters and random draws for this set
//     VectorXd a_prec_set(n_cols);
//     VectorXd e_prec_set(n_cols);
//     MatrixXd randn_draws_set(r,n_cols);
//     for(int j = 0; j < n_cols; j++){
//       int col_j = Y_cols[j]-1;
//       a_prec_set[j] = a_prec[col_j];
//       e_prec_set[j] = e_prec[col_j];
//       randn_draws_set.col(j) = randn_draws.col(col_j);
//     }
//
//     // sample U values for this set
//     MatrixXd U_set(r,n_cols);
//     sample_random_effects_worker sampler(Eta_set,Q,ZQt, s1, s2, a_prec_set, e_prec_set, randn_draws_set, U_set);
//     RcppParallel::parallelFor(0,n_cols,sampler,grainSize);
//
//     // copy new U values into full matrix
//     for(int j = 0; j < n_cols; j++){
//       U.col(Y_cols[j]-1) = U_set.col(j);
//     }
//   }
//   return(U);
// }
//
//
//
//
//
// // -------------------------------------------------- //
// // -- Sample factor scores --- //
// // -------------------------------------------------- //
//
// // Sample factor scores given factor loadings (F_a), factor residual variances (F_e_prec) and
// // phenotype residuals
// // Y - ZU = F * Lambda + E
// // F[,k] ~ N(XFBF[,k] + ZUF[,k],1/F_e_prec[,j])
// // E[,j] ~ N(0,1/resid_Eta_prec[,j])
// // Sampling is done separately for each block of rows with the same pattern of missing observations
// // [[Rcpp::export()]]
// MatrixXd sample_factors_scores_sparse_c_Eigen(
//     Map<MatrixXd> Eta_tilde,
//     Map<MatrixXd> prior_mean,
//     Map<MatrixXd> Lambda,
//     Map<VectorXd> resid_Eta_prec,
//     Map<VectorXd> F_e_prec,
//     Map<MatrixXd> randn_draws
// ) {
//   MatrixXd Lmsg = resid_Eta_prec.asDiagonal() * Lambda;
//   MatrixXd Sigma = Lambda.transpose() * Lmsg;
//   Sigma.diagonal() += F_e_prec;
//   Eigen::LLT<MatrixXd> chol_Sigma;
//   chol_Sigma.compute(Sigma);
//   MatrixXd R = chol_Sigma.matrixU();
//
//   MatrixXd Meta = R.transpose().triangularView<Lower>().solve((Eta_tilde * Lmsg + prior_mean * F_e_prec.asDiagonal()).transpose());
//
//   MatrixXd Ft = R.triangularView<Upper>().solve(Meta + randn_draws.transpose());
//
//   return Ft.transpose();
// }
//
//
// // Sample factor scores given factor loadings (F_a), factor residual variances (F_e_prec) and
// // phenotype residuals
// // Y - ZU = F * Lambda + E
// // F[,k] ~ N(XFBF[,k] + ZUF[,k],1/F_e_prec[,j])
// // E[,j] ~ N(0,1/resid_Eta_prec[,j])
// // Sampling is done separately for each block of rows with the same pattern of missing observations
// // [[Rcpp::export()]]
// MatrixXd sample_factors_scores_sparse_mising_c_Eigen(
//     Map<MatrixXd> Eta_tilde,
//     Map<MatrixXd> prior_mean,
//     Map<MatrixXd> Lambda,
//     Map<VectorXd> resid_Eta_prec,
//     Map<VectorXd> F_e_prec,
//     Map<MatrixXd> randn_draws,
//     Rcpp::List Y_row_obs_sets
// ) {
//
//   int k = randn_draws.rows();
//   int n = randn_draws.cols();
//
//   MatrixXd Ft(k,n);
//   MatrixXd Lmsg = resid_Eta_prec.asDiagonal() * Lambda;
//
//   std::vector<VectorXi> obs_sets;
//   std::vector<MatrixXd> R_s;
//   std::vector<MatrixXd> Lmsg_s;
//   for(int row_set = 0; row_set < Y_row_obs_sets.length(); row_set++){
//     Rcpp::List Y_row_obs_sets_i = Rcpp::as<Rcpp::List>(Y_row_obs_sets[row_set]);
//     VectorXi obs_set_i = as<VectorXi>(Y_row_obs_sets_i["columns"]);
//     VectorXi rows = as<VectorXi>(Y_row_obs_sets_i["rows"]);
//
//     int n_traits = obs_set_i.size();
//     int n_rows = rows.size();
//     MatrixXd Lmsg_set(n_traits,k);
//     MatrixXd Lambda_set(n_traits,k);
//     for(int j = 0; j < n_traits; j++){
//       int index = obs_set_i[j]-1;
//       Lmsg_set.row(j) = Lmsg.row(index);
//       Lambda_set.row(j) = Lambda.row(index);
//     }
//     MatrixXd Sigma_set = Lambda_set.transpose() * Lmsg_set;
//     Sigma_set.diagonal() += F_e_prec;
//     Eigen::LLT<MatrixXd> chol_Sigma_set;
//     chol_Sigma_set.compute(Sigma_set);
//     MatrixXd R_set = chol_Sigma_set.matrixU();
//
//     for(int i = 0; i < n_rows; i++){
//       int index = rows[i]-1;
//       Eigen::RowVectorXd Eta_i(n_traits);
//       for(int j = 0; j < n_traits; j++){
//         int trait_index = obs_set_i[j]-1;
//         Eta_i[j] = Eta_tilde.coeffRef(index,trait_index);
//       }
//       VectorXd Meta = R_set.transpose().triangularView<Lower>().solve((Eta_i * Lmsg_set + prior_mean.row(index) * F_e_prec.asDiagonal()).transpose());
//       Ft.col(index) = R_set.triangularView<Upper>().solve(Meta + randn_draws.col(index));
//     }
//   }
//   return Ft.transpose();
// }
//
//
//
//
// // -------------------------------------------------- //
// // -- calculate total prec scores --- //
// // -------------------------------------------------- //
//
// // [[Rcpp::export()]]
// VectorXd tot_prec_scores_c (
//     MatrixXd Y,
//     ArrayXXd resid_prec
// ) {
//
//   int p = Y.cols();
//
//   VectorXd scores(p);
//
//   for(int j = 0; j < p; j++){
//     ArrayXd Sigma_sqrt = resid_prec.col(j).sqrt().inverse();
//     VectorXd SiY_j = Y.col(j).array() / Sigma_sqrt;
//     scores(j) = SiY_j.dot(SiY_j);
//   }
//   return scores;
// }
//
// // [[Rcpp::export()]]
// VectorXd tot_prec_scores_missing_c (
//     Map<MatrixXd> Eta,
//     Map<VectorXd> h2,
//     Rcpp::List invert_aI_bZKZ
// ) {
//
//   int p = Eta.cols();
//   VectorXd scores(p);
//
//
//   for(int set = 0; set < invert_aI_bZKZ.size(); set++){
//     Rcpp::List invert_aI_bZKZ_i = Rcpp::as<Rcpp::List>(invert_aI_bZKZ[set]);
//     SpMat Qt = as<MSpMat>(invert_aI_bZKZ_i["Qt"]);
//     ArrayXd s = as<ArrayXd>(invert_aI_bZKZ_i["s"]);
//     VectorXi Y_obs = as<VectorXi>(invert_aI_bZKZ_i["Y_obs"]);
//     VectorXi Y_cols = as<VectorXi>(invert_aI_bZKZ_i["Y_cols"]);
//
//     int n_cols = Y_cols.size();
//     int n_obs = Y_obs.size();
//     for(int j = 0; j < n_cols; j++){
//       int col = Y_cols[j]-1;
//       VectorXd UtEta_j(n_obs);
//       for(int i = 0; i < n_obs; i++){
//         UtEta_j[i] = Eta.coeffRef(Y_obs[i]-1,col);
//       }
//       UtEta_j = Qt * UtEta_j;
//       ArrayXd Sigma_sqrt = sqrt(h2(col) * s + (1.0 - h2(col)));
//       VectorXd SiUtEta_j = UtEta_j.array() / Sigma_sqrt;
//       scores(col) = SiUtEta_j.dot(SiUtEta_j);
//     }
//   }
//   return scores;
// }
//
//
//
// // -------------------------------------------------- //
// // -- Sample h2s effects --- //
// // -------------------------------------------------- //
//
// // -- full scan of h2s --- //
//
// struct calc_log_ps_fast_worker : public Worker {
//   MatrixXd std_scores_b2;
//   VectorXd discrete_priors;
//   VectorXd tot_Eta_prec;
//   VectorXd s;
//   MatrixXd &log_ps;
//
//   calc_log_ps_fast_worker(MatrixXd std_scores_b2,
//                           VectorXd discrete_priors,
//                           VectorXd tot_Eta_prec,
//                           VectorXd s,
//                           MatrixXd &log_ps):
//     std_scores_b2(std_scores_b2), discrete_priors(discrete_priors), tot_Eta_prec(tot_Eta_prec),s(s), log_ps(log_ps) {}
//
//   void operator()(std::size_t begin, std::size_t end) {
//     int n = std_scores_b2.cols();
//     int b = discrete_priors.size();
//     for(std::size_t i = begin; i < end; i++){
//       double h2 = double(i)/b;
//       VectorXd s2s = h2*s.array() + (1.0-h2);
//       MatrixXd std_scores2 = -0.5 * std_scores_b2 * s2s.cwiseInverse().asDiagonal();
//       double det = -n/2 * log(2.0*M_PI) - 0.5*s2s.array().log().sum();
//       log_ps.row(i) = std_scores2.rowwise().sum().array() + det + log(discrete_priors[i]) + n/2*tot_Eta_prec.array().log();
//     }
//   }
// };
//
// // [[Rcpp::export()]]
// MatrixXd log_p_h2s_fast(
//     Map<MatrixXd> UtEta,
//     Map<VectorXd> tot_Eta_prec,
//     Map<VectorXd> discrete_priors,
//     Map<VectorXd> s,
//     int grainSize)
// {
//
//   int b = discrete_priors.size();
//   int p = UtEta.cols();
//
//   MatrixXd log_ps(b,p);
//
//   MatrixXd std_scores_b = tot_Eta_prec.cwiseSqrt().asDiagonal() * UtEta.transpose();
//   MatrixXd std_scores_b2 = std_scores_b.cwiseProduct(std_scores_b);
//
//   calc_log_ps_fast_worker sampler(std_scores_b2,discrete_priors,tot_Eta_prec,s,log_ps);
//   RcppParallel::parallelFor(0,b,sampler,grainSize);
//   return(log_ps);
// }
//
// // [[Rcpp::export()]]
// MatrixXd log_p_h2s_fast_missing(
//     Map<MatrixXd> Eta,
//     Map<VectorXd> tot_Eta_prec,
//     Map<VectorXd> discrete_priors,
//     Rcpp::List invert_aI_bZKZ,
//     int grainSize)
// {
//
//   int b = discrete_priors.size();
//   int p = Eta.cols();
//
//   MatrixXd log_ps(b,p);
//
//   // do by set of columns with same pattern of missingness
//   for(int col_set = 0; col_set < invert_aI_bZKZ.length(); col_set++){
//     Rcpp::List invert_aI_bZKZ_i = Rcpp::as<Rcpp::List>(invert_aI_bZKZ[col_set]);
//     SpMat Qt_set = as<MSpMat>(invert_aI_bZKZ_i["Qt"]);
//     VectorXd s_set = as<VectorXd>(invert_aI_bZKZ_i["s"]);
//     VectorXi Y_obs = as<VectorXi>(invert_aI_bZKZ_i["Y_obs"]);
//     VectorXi Y_cols = as<VectorXi>(invert_aI_bZKZ_i["Y_cols"]);
//
//     int n_obs = Y_obs.size();
//     int n_cols = Y_cols.size();
//     VectorXd Eta_std_j(n_obs);
//     MatrixXd std_scores_b2_set(n_cols,n_obs);
//     VectorXd tot_Eta_prec_set(n_cols);
//     for(int j = 0; j < n_cols; j++) {
//       int col = Y_cols[j]-1;
//       tot_Eta_prec_set(j) = tot_Eta_prec(col);
//       for(int i = 0; i < n_obs; i++){
//         Eta_std_j[i] = Eta.coeffRef(Y_obs[i]-1,col);
//       }
//       VectorXd QtEta_std_j = Qt_set * Eta_std_j * sqrt(tot_Eta_prec[col]);
//
//       std_scores_b2_set.row(j) = QtEta_std_j.cwiseProduct(QtEta_std_j);
//     }
//
//     MatrixXd log_ps_set(b,n_cols);
//
//     calc_log_ps_fast_worker sampler(std_scores_b2_set,discrete_priors,tot_Eta_prec_set,s_set,log_ps_set);
//     RcppParallel::parallelFor(0,b,sampler,grainSize);
//
//     // copy new log_ps_set values into full matrix
//     for(int j = 0; j < n_cols; j++){
//       log_ps.col(Y_cols[j]-1) = log_ps_set.col(j);
//     }
//   }
//   return(log_ps);
// }
//
//
// // -- using MH proposal
//
// double log_prob_h2_fast_c(
//     ArrayXd Qty,
//     ArrayXd s,
//     double h2,
//     int n,
//     double tot_Eta_prec,
//     double discrete_prior
// ){
//   ArrayXd s2s = h2*s + (1.0-h2);
//   double score2 = (Qty.pow(2)/s2s).sum() * tot_Eta_prec;
//   double log_det_Sigma = s2s.log().sum();
//
//   double log_p = -n/2.0 * log(2*M_PI) - 0.5*(log_det_Sigma - n*log(tot_Eta_prec)) - 0.5 * score2 + log(discrete_prior);
//   return log_p;
// }
//
// struct sample_h2_MH_fast_worker : public RcppParallel::Worker {
//   const MatrixXd QtEta;
//   const MatrixXd h2s_matrix;
//   const VectorXd s;
//   const VectorXd tot_Eta_prec;
//   const VectorXd discrete_priors;
//   const VectorXd r_draws;
//   const VectorXd state_draws;
//   const VectorXi h2_index;
//   const double step_size;
//   VectorXi &new_index;
//
//   sample_h2_MH_fast_worker(const MatrixXd QtEta,
//                            const MatrixXd h2s_matrix,
//                            const VectorXd s,
//                            const VectorXd tot_Eta_prec,
//                            const VectorXd discrete_priors,
//                            const VectorXd r_draws,
//                            const VectorXd state_draws,
//                            const VectorXi h2_index,
//                            const double step_size,
//                            VectorXi &new_index):
//
//     QtEta(QtEta), h2s_matrix(h2s_matrix),s(s),
//     tot_Eta_prec(tot_Eta_prec),  discrete_priors(discrete_priors),
//     r_draws(r_draws),state_draws(state_draws),h2_index(h2_index),step_size(step_size),new_index(new_index) {}
//
//   void operator()(std::size_t begin, std::size_t end) {
//     int n = QtEta.rows();
//     for(std::size_t j = begin; j < end; j++){
//       int old_state = h2_index[j] - 1;
//       double old_h2 = h2s_matrix.coeffRef(0,old_state);
//       double old_log_p = log_prob_h2_fast_c(QtEta.col(j),s,old_h2,n,tot_Eta_prec[j],discrete_priors[old_state]);
//
//       VectorXd candidate_new_states = find_candidate_states(h2s_matrix,step_size,old_state);
//       int r = state_draws[j] * (candidate_new_states.size());
//       int proposed_state = candidate_new_states[r];
//
//       double new_h2 = h2s_matrix.coeffRef(0,proposed_state);
//       double new_log_p = log_prob_h2_fast_c(QtEta.col(j),s,new_h2,n,tot_Eta_prec[j],discrete_priors[proposed_state]);
//
//       VectorXd candidate_states_from_new_state = find_candidate_states(h2s_matrix,step_size,proposed_state);
//
//       double forward_prob = 1.0 / candidate_new_states.size();
//       double back_prob = 1.0 / candidate_states_from_new_state.size();
//
//       double log_MH_ratio = new_log_p - old_log_p + log(forward_prob) - log(back_prob);
//
//       if(log(r_draws[j]) < log_MH_ratio) {
//         new_index[j] = proposed_state;
//       } else {
//         new_index[j] = old_state;
//       }
//     }
//   }
// };
//
//
// // [[Rcpp::export()]]
// VectorXi sample_h2s_discrete_MH_fast_c(
//     Map<MatrixXd> UtEta,
//     Map<VectorXd> tot_Eta_prec,
//     Map<VectorXd> discrete_priors,
//     VectorXi h2_index,
//     Map<MatrixXd> h2s_matrix,
//     Map<VectorXd> s,
//     Map<VectorXd> r_draws,
//     Map<VectorXd> state_draws,
//     double step_size,
//     int grainSize
// ){
//
//   int p = UtEta.cols();
//
//   VectorXi new_index(p);
//
//   sample_h2_MH_fast_worker sampler(UtEta,h2s_matrix,s,tot_Eta_prec,discrete_priors,r_draws,state_draws,h2_index,step_size,new_index);
//   RcppParallel::parallelFor(0,p,sampler,grainSize);
//   return new_index;
// }
//
//
// // [[Rcpp::export()]]
// VectorXi sample_h2s_discrete_MH_fast_missing_c(
//     Map<MatrixXd> Eta,
//     Map<VectorXd> tot_Eta_prec,
//     Map<VectorXd> discrete_priors,
//     VectorXi h2_index,
//     Map<MatrixXd> h2s_matrix,
//     Map<VectorXd> r_draws,
//     Map<VectorXd> state_draws,
//     double step_size,
//     Rcpp::List invert_aI_bZKZ,
//     int grainSize
// ){
//
//   int p = Eta.cols();
//
//   VectorXi new_index(p);
//
//   for(int col_set = 0; col_set < invert_aI_bZKZ.length(); col_set++){
//     Rcpp::List invert_aI_bZKZ_i = Rcpp::as<Rcpp::List>(invert_aI_bZKZ[col_set]);
//     SpMat Qt_set = as<MSpMat>(invert_aI_bZKZ_i["Qt"]);
//     VectorXd s_set = as<VectorXd>(invert_aI_bZKZ_i["s"]);
//     VectorXi Y_obs = as<VectorXi>(invert_aI_bZKZ_i["Y_obs"]);
//     VectorXi Y_cols = as<VectorXi>(invert_aI_bZKZ_i["Y_cols"]);
//
//     int n_obs = Y_obs.size();
//     int n_cols = Y_cols.size();
//
//     MatrixXd Eta_set(n_obs,n_cols);
//     VectorXd tot_Eta_prec_set(n_cols);
//     VectorXd r_draws_set(n_cols);
//     VectorXd state_draws_set(n_cols);
//     VectorXi h2_index_set(n_cols);
//     for(int j = 0; j < n_cols; j++) {
//       int col = Y_cols[j]-1;
//       tot_Eta_prec_set(j) = tot_Eta_prec(col);
//       r_draws_set(j)      = r_draws(col);
//       state_draws_set(j)  = state_draws(col);
//       h2_index_set(j)     = h2_index(col);
//       for(int i = 0; i < n_obs; i++){
//         Eta_set.coeffRef(i,j) = Eta.coeffRef(Y_obs[i]-1,col);
//       }
//     }
//
//     VectorXi new_index_set(n_cols);
//     MatrixXd QtEta_set = Qt_set * Eta_set;
//
//     sample_h2_MH_fast_worker sampler(QtEta_set,h2s_matrix,s_set,tot_Eta_prec_set,discrete_priors,r_draws_set,state_draws_set,h2_index_set,step_size,new_index_set);
//     RcppParallel::parallelFor(0,n_cols,sampler,grainSize);
//
//     for(int j = 0; j < n_cols; j++){
//       new_index(Y_cols[j]-1) = new_index_set(j);
//     }
//   }
//   return new_index;
// }
