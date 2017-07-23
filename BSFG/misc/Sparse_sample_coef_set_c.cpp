// #include <math.h>
// #include <iostream>
// #include "BSFG_types.h"
// #include<Eigen/SparseCholesky>
//
// VectorXd sample_coefs_single2(
//     VectorXd UtEta,
//     SpMat UtW,
//     VectorXd prior_mean,
//     VectorXd prior_prec,
//     double h2,
//     double tot_Eta_prec,
//     VectorXd randn_theta,
//     VectorXd randn_e,
//     VectorXd s,
//     int b,
//     int n
// ) {
//
//   VectorXd R_sq_diag = ((h2 * s.array() + (1.0-h2))/tot_Eta_prec).sqrt();
//   VectorXd theta_star = randn_theta.array()/prior_prec.array().sqrt();
//   theta_star += prior_mean;
//   VectorXd e_star = randn_e.array() * R_sq_diag.array();
//   MatrixXd UtW_theta_star = UtW * theta_star;
//   VectorXd eta_resid = UtEta - UtW_theta_star - e_star;
//   SpMat RinvSqUtW = R_sq_diag.cwiseInverse().asDiagonal() * UtW;
//   VectorXd eta_std = eta_resid.array()/R_sq_diag.array();
//   VectorXd WtURinvy = RinvSqUtW.transpose() * eta_std;
//
//   VectorXd theta_tilda;
//   if(b < n) {
//     MatrixXd C = RinvSqUtW.transpose() * RinvSqUtW;
//     C.diagonal() = C.diagonal() + prior_prec;
//     theta_tilda = C.householderQr().solve(WtURinvy);
//   } else{
//     SpMat VAi = UtW * prior_prec.cwiseInverse().asDiagonal();
//     SpMat inner = VAi*UtW.transpose();
//     for(int i = 0; i < n; i++) {
//       inner.coeffRef(i,i) += (h2 * s(i) + (1.0-h2)) / tot_Eta_prec;
//     }
//     VectorXd VAiWtURinvy = VAi * WtURinvy;
//     SimplicialLDLT <SparseMatrix<double>, Eigen::Upper > solver;
//     VectorXd outerWtURinvy = VAi.transpose() * solver.compute(inner).solve(VAiWtURinvy);
//     theta_tilda = WtURinvy.array() / prior_prec.array();
//     theta_tilda -= outerWtURinvy;
//   }
//   VectorXd coefs = theta_tilda + theta_star;
//   return coefs;
// }
//
// // [[Rcpp::export()]]
// List sample_coefs_set_c2(
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
//   std::vector<MatrixXd> UtEta_list;
//   std::vector<SpMat> UtW_list;
//   std::vector<VectorXd> s_list;
//   std::vector<VectorXd> position_list;
//   std::vector<MatrixXd> randn_theta_list;
//   std::vector<MatrixXd> randn_e_list;
//   int randn_theta_index = 0;
//   int randn_e_index = 0;
//   int total_obs = 0;
//   for(int i = 0; i < n; i++){
//     Rcpp::List model_matrix_i = Rcpp::as<Rcpp::List>(model_matrices[i]);
//     UtEta_list.push_back(Rcpp::as<MatrixXd>(model_matrix_i["y"]));
//     UtW_list.push_back(Rcpp::as<SpMat>(model_matrix_i["X"]));
//     position_list.push_back(Rcpp::as<VectorXd>(model_matrix_i["position"]));
//     s_list.push_back(Rcpp::as<VectorXd>(model_matrix_i["s"]));
//
//     int p = UtEta_list[i].cols();
//     int n_obs = UtEta_list[i].rows();
//     int b = UtW_list[i].cols();
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
//   int n_traits = UtEta_list[0].cols();
//   int b = randn_theta_list[0].rows()*n_traits;
//
//   MatrixXd coefs(b,n);
//   MatrixXd Y_fitted(total_obs,n_traits);
//
//   struct sampleColumn : public RcppParallel::Worker {
//     std::vector<MatrixXd> UtEta_list;
//     std::vector<SpMat> UtW_list;
//     std::vector<MatrixXd> randn_theta_list;
//     std::vector<MatrixXd> randn_e_list;
//     std::vector<VectorXd> position_list;
//     std::vector<VectorXd> s_list;
//     MatrixXd h2s,tot_Eta_prec;
//     MatrixXd prior_mean,prior_prec;
//     int n_traits;
//     MatrixXd &coefs;
//     MatrixXd &Y_fitted;
//
//     sampleColumn(
//       std::vector<MatrixXd> UtEta_list,
//       std::vector<SpMat> UtW_list,
//       std::vector<MatrixXd> randn_theta_list,
//       std::vector<MatrixXd> randn_e_list,
//       std::vector<VectorXd> position_list,
//       std::vector<VectorXd> s_list,
//       MatrixXd h2s,
//       MatrixXd tot_Eta_prec,
//       MatrixXd prior_mean,
//       MatrixXd prior_prec,
//       int n_traits,
//       MatrixXd &coefs,
//       MatrixXd &Y_fitted) :
//       UtEta_list(UtEta_list), UtW_list(UtW_list), randn_theta_list(randn_theta_list), randn_e_list(randn_e_list),position_list(position_list),s_list(s_list),
//       h2s(h2s), tot_Eta_prec(tot_Eta_prec),prior_mean(prior_mean), prior_prec(prior_prec),n_traits(n_traits),
//       coefs(coefs),Y_fitted(Y_fitted)
//     {}
//
//     void operator()(std::size_t begin, std::size_t end) {
//       for(std::size_t j = begin; j < end; j++){
//         int b = randn_theta_list[j].rows();
//         int n_obs = randn_e_list[j].rows();
//         for(int t = 0; t < n_traits; t++) {
//           coefs.block(t*b,j,b,1) = sample_coefs_single2(UtEta_list[j].col(t), UtW_list[j], prior_mean.block(t*b,j,b,1),
//                       prior_prec.block(t*b,j,b,1), h2s(j,t), tot_Eta_prec(j,t), randn_theta_list[j].col(t),randn_e_list[j].col(t),s_list[j],b,n_obs);
//         }
//         Map<MatrixXd> Eta_i(coefs.col(j).data(),b,n_traits);
//         MatrixXd Y_fitted_j = UtW_list[j] * Eta_i;
//         for(int i = 0; i < position_list[j].size(); i++) {
//           Y_fitted.row(position_list[j][i]-1) = Y_fitted_j.row(i);
//         }
//       }
//     }
//   };
//
//
//   sampleColumn sampler(UtEta_list, UtW_list,randn_theta_list,randn_e_list,position_list,s_list,h2s,tot_Eta_prec,prior_mean,prior_prec,n_traits,coefs,Y_fitted);
//   RcppParallel::parallelFor(0,n,sampler,grainSize);
//
//   return(Rcpp::List::create(
//       Rcpp::Named("coefs") = coefs,
//       Rcpp::Named("Y_fitted") = Y_fitted));
// }
