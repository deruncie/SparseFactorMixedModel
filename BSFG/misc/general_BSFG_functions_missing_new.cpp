// #include <math.h>
// #include <iostream>
// #include "BSFG_types.h"
//
// // [[Rcpp::depends(RcppEigen)]]
// using namespace Eigen;
// using namespace RcppParallel;
//
// MatrixXd subset_MatrixXd_block(MatrixXd Y, VectorXi Y_obs, VectorXi Y_cols){
//   int n_obs = Y_obs.size();
//   int n_cols = Y_cols.size();
//   MatrixXd Y_set(n_obs,n_cols);
//   for(int i = 0; i < n_obs; i++){
//     int row_i = Y_obs[i]-1;
//     for(int j = 0; j < n_cols; j++){
//       int col_j = Y_cols[j]-1;
//       Y_set.coeffRef(i,j) = Y.coeffRef(row_i,col_j);
//     }
//   }
//   return(Y_set);
// }
//
// VectorXd subset_VectorXd(VectorXd x, VectorXi indices){
//   int n = indices.size();
//   VectorXd out(n);
//   for(int i = 0; i < n; i++){
//     out(i) = x[indices[i]-1];
//   }
//   return(out);
// }
//
// VectorXi subset_VectorXi(VectorXi x, VectorXi indices){
//   int n = indices.size();
//   VectorXi out(n);
//   for(int i = 0; i < n; i++){
//     out(i) = x[indices[i]-1];
//   }
//   return(out);
// }
//
// MatrixXd subset_MatrixXd_cols(MatrixXd Y, VectorXi Y_cols) {
//   int n = Y.rows();
//   int n_cols = Y_cols.size();
//   MatrixXd Y_sub(n,n_cols);
//   for(int j = 0; j < n_cols; j++){
//     Y_sub.col(j) = Y.col(Y_cols[j]-1);
//   }
//   return(Y_sub);
// }
//
//
// MatrixXd subset_MatrixXd_rows(MatrixXd Y, VectorXi Y_rows) {
//   MatrixXd Y_sub = subset_MatrixXd_cols(Y.transpose(),Y_rows).transpose();
//   return(Y_sub);
//   // int p = Y.cols();
//   // int n_rows = Y_rows.size();
//   // MatrixXd Y_sub(n_rows,p);
//   // for(int i = 0; i < n_rows; i++){
//   //   Y_sub.row(i) = Y.row(Y_rows[i]-1);
//   // }
//   // return(Y_sub);
// }
// // missing data version - sets of columns of Y have different patterns of missing data
// // each set can be sampled at once using sample_MME_single_diagK_worker
// // [[Rcpp::export()]]
// MatrixXd sample_MME_fixedEffects_missing_c(  // returns bxp matrix
//     Map<MatrixXd> Y,                // nxp
//     Rcpp::List X_sets,              // List. Each element contains a n_obs x b matrix, a subset of the full X matrix
//     Rcpp::List Sigma_Cholesky_sets, // List. Each element contains: List of Cholesky decompositions of residuals. nxn upper triangular, dgCMatrix.
//     Rcpp::List Missing_data_map,    // List. Each element containts Y_obs, Y_cols. The length should correspond to the length of Sigma_Choleksy_sets
//     VectorXi h2s_index,             // px1 index of Cholesky matrix for each column
//     Map<VectorXd> tot_Eta_prec,     // px1
//     Map<MatrixXd> prior_mean,       // bxp
//     Map<MatrixXd> prior_prec,       // bxp
//     int grainSize) {
//
//   int b = prior_mean.cols();
//   int p = Y.cols();
//
//   MatrixXd coefs(b,p);
//
//   // sample sets of columns with same pattern of missing data as a block
//   for(int col_set = 0; col_set < Missing_data_map.length(); col_set++){
//     Rcpp::List Missing_data_map_set = Rcpp::as<Rcpp::List>(Missing_data_map[col_set]);
//
//     VectorXi Y_obs = as<Map<VectorXi>>(Missing_data_map_set["Y_obs"]);
//     VectorXi Y_cols = as<Map<VectorXi>>(Missing_data_map_set["Y_cols"]);
//
//     int n_obs = Y_obs.size();
//     int n_cols = Y_cols.size();
//
//     // pull out rows of X corresponding to observed individuals
//     MatrixXd X_set = as<MatrixXd>(X_sets[col_set]);
//
//     // pull out elements of Y corresponding to observed individuals
//     MatrixXd Y_set = subset_MatrixXd_block(Y,Y_obs,Y_cols);
//
//     // pull out prior_mean, prior_prec, randn_theta, randn_e for col_set
//     MatrixXd prior_mean_set = subset_MatrixXd_cols(prior_mean,Y_cols);
//     MatrixXd prior_prec_set = subset_MatrixXd_cols(prior_prec,Y_cols);
//
//     // pull out h2s_index_set and tot_Eta_prec_set
//     VectorXi h2s_index_set = subset_VectorXi(h2s_index, Y_cols);
//     VectorXd tot_Eta_prec_set = subset_VectorXd(tot_Eta_prec, Y_cols);
//
//     // pull out necessary chol_R matrices for set
//     Rcpp::List Sigma_Choleskys_set = Rcpp::as<Rcpp::List>(Sigma_Cholesky_sets[col_set]);
//     std::vector<MSpMat> chol_R_list_set;
//     for(int i = 0; i < h2s_index_set.maxCoeff(); i++){
//       Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys_set[i]);
//       chol_R_list_set.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
//     }
//
//     // draw standard normals
//     MatrixXd randn_theta_set = rstdnorm_mat(b,n_cols);
//     MatrixXd randn_e_set = rstdnorm_mat(n_obs,n_cols);
//
//     MatrixXd coefs_set(b,n_cols);
//
//     // sample coefs_set as a block
//     sample_MME_single_diagK_worker sampler(Y_set,X_set,prior_mean_set,prior_prec_set,chol_R_list_set,h2s_index_set,tot_Eta_prec_set,randn_theta_set,randn_e_set,coefs_set);
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
// // samples random effects from model:
// // Y = ZU + E
// // U[,j] ~ N(0,1/a_prec[j] * K)
// // E[,j] ~ N(0,1/e_prec[j] * I_n)
// // Where matrix Q and vectors s1 and s2 diagonalize ZZt + K
// // invert_aZZt_Kinv has parameters to diagonalize a*Z*Z' + b*K for fast inversion
// // with missing observations; Samples sets of columns with same patterns of missingness
// // [[Rcpp::export()]]
// MatrixXd sample_MME_ZKZts_missing_c (  // returns rxp matrix
//     Map<MatrixXd> Y,                    // nxp
//     Rcpp::List Z_sets,                  // List. Each element contains a n_obs x r dgCMatrix, a subset of rows of the full Z
//     Map<VectorXd> tot_Eta_prec,         // px1
//     Rcpp::List randomEffect_C_Cholesky_sets, // List. Each element contains: List of pairs chol_C and chol_K_inv (both rxr dgCMatrix upper-triangle)
//     Rcpp::List Missing_data_map,        // List. Each element containts Y_obs, Y_cols. The length should correspond to the length of randomEffect_C_Cholesky_sets
//     Map<MatrixXd> h2s,                  // n_RE x p
//     Map<VectorXi> h2s_index,            // px1
//     int grainSize) {
//
//   int p = Y.cols();
//   MSpMat Z_set1 = as<MSpMat>(Z_sets[0]);
//   int r = Z_set1.cols();
//
//   MatrixXd U(r,p);
//
//   // process sets of columns with identical patterns of missing observations
//   for(int col_set = 0; col_set < Missing_data_map.length(); col_set++){
//     Rcpp::List Missing_data_map_set = Rcpp::as<Rcpp::List>(Missing_data_map[col_set]);
//     Rcpp::List randomEffect_C_Choleskys_set = Rcpp::as<Rcpp::List>(randomEffect_C_Cholesky_sets[col_set]);
//
//     VectorXi Y_obs = as<Map<VectorXi>>(Missing_data_map_set["Y_obs"]);
//     VectorXi Y_cols = as<Map<VectorXi>>(Missing_data_map_set["Y_cols"]);
//
//     int n_obs = Y_obs.size();
//     int n_cols = Y_cols.size();
//
//     // pull out rows of Z corresponding to observed individuals
//     MSpMat Z_set = as<MSpMat>(Z_sets[col_set]);
//
//     // pull out elements of Y corresponding to observed individuals
//     MatrixXd Y_set = subset_MatrixXd_block(Y,Y_obs,Y_cols);
//
//     // draw random numbers
//     MatrixXd randn_theta_set = rstdnorm_mat(r,n_cols);
//     MatrixXd randn_e_set = rstdnorm_mat(n_obs,n_cols);
//
//     VectorXi h2s_index_set = subset_VectorXi(h2s_index, Y_cols);
//     VectorXd tot_Eta_prec_set = subset_VectorXd(tot_Eta_prec, Y_cols);
//     MatrixXd h2s_set = subset_MatrixXd_cols(h2s,Y_cols);
//     VectorXd h2_e_set = 1.0 - h2s_set.colwise().sum().array();
//     VectorXd pes_set = tot_Eta_prec_set.array() / h2_e_set.array();
//
//     std::vector<MSpMat> chol_C_list_set,chol_K_inv_list_set;
//     for(int i = 0; i < h2s_index_set.maxCoeff(); i++){
//       Rcpp::List randomEffect_C_Cholesky_i = Rcpp::as<Rcpp::List>(randomEffect_C_Choleskys_set[i]);
//       chol_C_list_set.push_back(Rcpp::as<MSpMat>(randomEffect_C_Cholesky_i["chol_C"]));
//       chol_K_inv_list_set.push_back(Rcpp::as<MSpMat>(randomEffect_C_Cholesky_i["chol_K_inv"]));
//     }
//
//     // sample U values for this set
//     MatrixXd U_set(r,n_cols);
//     sample_MME_single_diagR_worker sampler(Y_set,Z_set,chol_C_list_set,chol_K_inv_list_set,pes_set,tot_Eta_prec_set,h2s_index_set,randn_theta_set,randn_e_set,U_set);
//     RcppParallel::parallelFor(0,p,sampler,grainSize);
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
// // -------------------------------------------- //
// // -------------- tot_prec_scores ------------- //
// // -------------------------------------------- //
//
// // calculates normal scores for:
// // Y_resid %*% Sigma_inv %*% t(Y_resid)
// // works on sets of columns of Y with same patterns of missingness
// // [[Rcpp::export()]]
// VectorXd tot_prec_scores_missing_c (
//     Map<MatrixXd> Y,                // nxp
//     Map<VectorXi> h2s_index,        // px1
//     Rcpp::List Sigma_Cholesky_sets, // List. Each element contains: List of Cholesky decompositions of residuals. nxn upper triangular, dgCMatrix.
//     Rcpp::List Missing_data_map,    // List. Each element containts Y_obs, Y_cols. The length should correspond to the length of Sigma_Choleksy_sets
//     int grainSize
// ) {
//
//   int p = Y.cols();
//   VectorXd scores(p);
//
//
//   for(int col_set = 0; col_set < Missing_data_map.size(); col_set++){
//     Rcpp::List Missing_data_map_set = Rcpp::as<Rcpp::List>(Missing_data_map[col_set]);
//
//     VectorXi Y_obs = as<Map<VectorXi>>(Missing_data_map_set["Y_obs"]);
//     VectorXi Y_cols = as<Map<VectorXi>>(Missing_data_map_set["Y_cols"]);
//
//     int n_cols = Y_cols.size();
//
//     // pull out elements of Y corresponding to observed individuals
//     MatrixXd Y_set = subset_MatrixXd_block(Y,Y_obs,Y_cols);
//
//     // pull out h2s_index_set
//     VectorXi h2s_index_set = subset_VectorXi(h2s_index, Y_cols);
//
//     // pull out necessary chol_R matrices for set
//     Rcpp::List Sigma_Choleskys_set = Rcpp::as<Rcpp::List>(Sigma_Cholesky_sets[col_set]);
//     std::vector<MSpMat> chol_R_list_set;
//     for(int i = 0; i < h2s_index_set.maxCoeff(); i++){
//       Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys_set[i]);
//       chol_R_list_set.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
//     }
//
//     VectorXd scores_set(n_cols);
//
//     tot_prec_scores_worker sampler(Y_set,chol_R_list_set,h2s_index_set,scores_set);
//     RcppParallel::parallelFor(0,n_cols,sampler,grainSize);
//
//     // copy scores_set into scores
//     for(int j = 0; j < n_cols; j++){
//       scores(Y_cols[j]-1) = scores_set(j);
//     }
//
//   }
//   return scores;
// }
//
//
// // -------------------------------------------- //
// // ---------------- sample h2s ---------------- //
// // -------------------------------------------- //
// // [[Rcpp::export()]]
// MatrixXd log_p_h2s_missing(
//     Map<MatrixXd> Y,              // nxp
//     Map<VectorXd> tot_Eta_prec,   // px1
//     Rcpp::List Sigma_Cholesky_sets, // List. Each element contains: List. Each element contains: chol_Sigma:. nxn upper triangular, dgCMatrix, log_det_Sigma
//     Rcpp::List Missing_data_map,    // List. Each element containts Y_obs, Y_cols. The length should correspond to the length of Sigma_Choleksy_sets
//     Map<VectorXd> discrete_priors,// n_h2 x 1
//     int grainSize)
// {
//
//   int b = discrete_priors.size();
//   int p = Y.cols();
//   MatrixXd log_ps(b,p);
//
//   for(int col_set = 0; col_set < Missing_data_map.size(); col_set++){
//     Rcpp::List Missing_data_map_set = Rcpp::as<Rcpp::List>(Missing_data_map[col_set]);
//
//     VectorXi Y_obs = as<Map<VectorXi>>(Missing_data_map_set["Y_obs"]);
//     VectorXi Y_cols = as<Map<VectorXi>>(Missing_data_map_set["Y_cols"]);
//
//     int n_cols = Y_cols.size();
//
//     // pull out elements of Y corresponding to observed individuals
//     MatrixXd Y_set = subset_MatrixXd_block(Y,Y_obs,Y_cols);
//
//     VectorXd tot_Eta_prec_set = subset_VectorXd(tot_Eta_prec, Y_cols);
//
//     Rcpp::List Sigma_Choleskys_set = Rcpp::as<Rcpp::List>(Sigma_Cholesky_sets[col_set]);
//     std::vector<MSpMat> chol_R_list_set;
//     VectorXd log_det_Sigmas_set(b);
//     for(int i = 0; i < b; i++){
//       Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys_set[i]);
//       chol_R_list_set.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
//       log_det_Sigmas_set[i] = Rcpp::as<double>(Sigma_Choleskys_i["log_det"]);
//     }
//
//     MatrixXd log_ps_set(b,n_cols);
//
//     log_ps_worker sampler(Y_set,tot_Eta_prec_set,chol_R_list_set,log_det_Sigmas_set,discrete_priors,log_ps_set);
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
// // -------------------------------------------- //
// // --------------- sample h2s MH -------------- //
// // -------------------------------------------- //
// // [[Rcpp::export()]]
// VectorXi sample_h2s_discrete_MH_missing_c(
//     Map<MatrixXd> Y,                // nxp
//     Map<VectorXd> tot_Eta_prec,     // px1
//     Map<VectorXd> discrete_priors,  // n_h2 x 1
//     VectorXi h2s_index,             // px1
//     Map<MatrixXd> h2s_matrix,       // n_RE x n_h2
//     Rcpp::List Sigma_Cholesky_sets, // List. Each element contains: List. Each element contains: chol_Sigma:. nxn upper triangular, dgCMatrix, log_det_Sigma
//     Rcpp::List Missing_data_map,    // List. Each element containts Y_obs, Y_cols. The length should correspond to the length of Sigma_Choleksy_sets
//     double step_size,               // double
//     int grainSize
// ){
//
//   int p = Y.cols();
//   int b = discrete_priors.size();
//   VectorXi new_index(p);
//
//   for(int col_set = 0; col_set < Missing_data_map.size(); col_set++){
//     Rcpp::List Missing_data_map_set = Rcpp::as<Rcpp::List>(Missing_data_map[col_set]);
//
//     VectorXi Y_obs = as<Map<VectorXi>>(Missing_data_map_set["Y_obs"]);
//     VectorXi Y_cols = as<Map<VectorXi>>(Missing_data_map_set["Y_cols"]);
//
//     int n_cols = Y_cols.size();
//
//     // pull out elements of Y corresponding to observed individuals
//     MatrixXd Y_set = subset_MatrixXd_block(Y,Y_obs,Y_cols);
//
//     VectorXi h2s_index_set = subset_VectorXi(h2s_index, Y_cols);
//     VectorXd tot_Eta_prec_set = subset_VectorXd(tot_Eta_prec, Y_cols);
//
//     Rcpp::List Sigma_Choleskys_set = Rcpp::as<Rcpp::List>(Sigma_Cholesky_sets[col_set]);
//     std::vector<MSpMat> chol_R_list_set;
//     VectorXd log_det_Sigmas_set(b);
//     for(int i = 0; i < b; i++){
//       Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys_set[i]);
//       chol_R_list_set.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
//       log_det_Sigmas_set[i] = Rcpp::as<double>(Sigma_Choleskys_i["log_det"]);
//     }
//
//     VectorXd r_draws_set = as<VectorXd>(runif(n_cols));
//     VectorXd state_draws_set = as<VectorXd>(runif(n_cols));
//     VectorXi new_index_set(n_cols);
//
//     sample_h2s_discrete_MH_worker sampler(Y_set,h2s_matrix,chol_R_list_set,log_det_Sigmas_set,tot_Eta_prec_set,discrete_priors,r_draws_set,state_draws_set,h2s_index_set,step_size,new_index_set);
//     RcppParallel::parallelFor(0,n_cols,sampler,grainSize);
//
//     // copy new_index_set values into new_index
//     for(int j = 0; j < n_cols; j++){
//       new_index(Y_cols[j]-1) = new_index_set(j);
//     }
//   }
//
//   return(new_index);
// }
