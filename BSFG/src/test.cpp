// #include <math.h>
// #include <iostream>
// #include "BSFG_types.h"
// // #include<Eigen/SparseCholesky>
//
// // [[Rcpp::export()]]
// MatrixXd add_Sp_Mat(SpMat X1,Map<MatrixXd> X2){
//   return(X1+X2);
// }
// // [[Rcpp::export()]]
// SpMat add_Sp_Sp(SpMat X1,SpMat X2){
//   return(X1+X2);
// }
// // [[Rcpp::export()]]
// MatrixXd to_Mat(SpMat X1){
//   return(X1.toDense());
// }

//
//
//
//
//
//
// VectorXd sample_coefs_single_hierarchical(
//     VectorXd UtEta,
//     SpMat UtW,
//     SpMat UtWX,
//     SpMat X,
//     VectorXd prior_mean,
//     VectorXd prior_prec,
//     double h2,
//     double tot_Eta_prec,
//     VectorXd randn_theta,
//     VectorXd randn_e,
//     VectorXd s,
//     int b,
//     int n,
//     int r
// ) {
//
//   VectorXd R_sq_diag = ((h2 * s.array() + (1.0-h2))/tot_Eta_prec).sqrt();
//   VectorXd theta_star = randn_theta.array()/prior_prec.array().sqrt();
//   theta_star += prior_mean;
//   VectorXd e_star = randn_e.array() * R_sq_diag.array();
//   MatrixXd UtWX_theta_star = UtWX*theta_star;
//   VectorXd eta_resid = UtEta - UtWX_theta_star - e_star;
//   SpMat RinvSqUtWX = R_sq_diag.cwiseInverse().asDiagonal() * UtWX;
//   VectorXd eta_std = eta_resid.array()/R_sq_diag.array();
//   VectorXd XtWtURinvy = RinvSqUtWX.transpose() * eta_std;
//
//   VectorXd theta_tilda;
//   if(b < r) {
//     MatrixXd C = RinvSqUtWX.transpose() * RinvSqUtWX;
//     C.diagonal() = C.diagonal() + prior_prec;
//     theta_tilda = C.householderQr().solve(XtWtURinvy);
//   } else{
//     SpMat B = UtW.transpose() * R_sq_diag.array().pow(2).inverse().matrix().asDiagonal() * UtW;
//     SimplicialLLT<SparseMatrix<double> >solver;
//     SpMat I(r,r);
//     I.setIdentity();
//     SpMat Bi = solver.compute(B).solve(I);
//     MatrixXd VAi = X * prior_prec.cwiseInverse().asDiagonal();
//     MatrixXd inner = VAi * X.transpose() + Bi;
//     VectorXd VAiXtWtURinvy = VAi * XtWtURinvy;
//     VectorXd outerXtWtURinvy = VAi.transpose() * inner.householderQr().solve(VAiXtWtURinvy);
//     theta_tilda = XtWtURinvy.array() / prior_prec.array();
//     theta_tilda -= outerXtWtURinvy;
//   }
//   VectorXd coefs = theta_tilda + theta_star;
//   return coefs;
// }
//
// // [[Rcpp::export()]]
// MatrixXd sample_coefs_hierarchical_parallel_sparse_c_Eigen(
//     MSpMat Ut,
//     Map<MatrixXd> Eta,
//     MSpMat W,
//     MSpMat X,
//     Map<VectorXd> h2,
//     Map<VectorXd> tot_Eta_prec,
//     Map<VectorXd> s,
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
//
//   MatrixXd UtEta = Ut*Eta;
//   SpMat UtW = Ut*W;
//   SpMat UtWX = UtW*X;
//
//   int p = UtEta.cols();
//   int b = X.cols();
//   int n = UtW.rows();
//   int r = X.rows();
//
//   MatrixXd coefs(b,p);
//
//   struct sampleColumn : public RcppParallel::Worker {
//     SpMat UtW,UtWX,X;
//     MatrixXd UtEta;
//     MatrixXd prior_mean, prior_prec, randn_theta, randn_e;
//     VectorXd h2, tot_Eta_prec, s;
//     int b,n,r;
//     MatrixXd &coefs;
//
//     sampleColumn(SpMat UtW, SpMat UtWX, SpMat X,MatrixXd UtEta, MatrixXd prior_mean, MatrixXd prior_prec, VectorXd h2, VectorXd tot_Eta_prec,
//                  MatrixXd randn_theta, MatrixXd randn_e,
//                  VectorXd s, int b, int n,int r,
//                  MatrixXd &coefs) :
//       UtW(UtW), UtWX(UtWX), X(X),UtEta(UtEta), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
//       h2(h2), tot_Eta_prec(tot_Eta_prec),
//       s(s), b(b), n(n),r(r),
//       coefs(coefs) {}
//
//     void operator()(std::size_t begin, std::size_t end) {
//       for(std::size_t j = begin; j < end; j++){
//       // for(int j = 0; j < 1; j++){
//         coefs.col(j) = sample_coefs_single_hierarchical(UtEta.col(j), UtW, UtWX, X,prior_mean.col(j), prior_prec.col(j), h2(j), tot_Eta_prec(j), randn_theta.col(j),randn_e.col(j),s,b,n,r);
//       }
//     }
//   };
//
//   sampleColumn sampler(UtW, UtWX, X,UtEta, prior_mean, prior_prec, h2, tot_Eta_prec, randn_theta,randn_e, s, b, n, r,coefs);
//   RcppParallel::parallelFor(0,p,sampler,grainSize);
//
//   return(coefs);
// }
//
//
//
// //
// //
// // // -------------------------------------------- //
// // // --------------- sample h2s MH fast --------- //
// // // -------------------------------------------- //
// //
// // // [[Rcpp::export()]]
// // double log_prob_h2_fast_c(
// //     ArrayXd Uty,
// //     ArrayXd s,
// //     double h2,
// //     int n,
// //     double tot_Eta_prec,
// //     double discrete_prior
// // ){
// //   ArrayXd s2s = h2*s + (1.0-h2);
// //   double score2 = (Uty.pow(2)/s2s).sum() * tot_Eta_prec;
// //   double log_det_Sigma = s2s.log().sum();
// //
// //   double log_p = -n/2.0 * log(2*M_PI) - 0.5*(log_det_Sigma - n*log(tot_Eta_prec)) - 0.5 * score2 + log(discrete_prior);
// //   return log_p;
// // }
// //
// //
// // // [[Rcpp::export()]]
// // VectorXd find_candidate_states(
// //     MatrixXd h2s_matrix,
// //     double step_size,
// //     int old_state
// // ) {
// //   VectorXd dists = (h2s_matrix.colwise() - h2s_matrix.col(old_state)).cwiseAbs().colwise().sum();
// //   VectorXd indices(dists.size());
// //   int count = 0;
// //   for(int i = 0; i < dists.size(); i++){
// //     if(dists[i] < step_size & dists[i] > 0) {
// //       indices[count] = i;
// //       count++;
// //     }
// //   }
// //   if(count == 0) {  // return all indices as candidates
// //     for(int i = 0; i < dists.size(); i++){
// //       indices[count] = i;
// //       count++;
// //     }
// //   }
// //   return indices.head(count);
// // }
// //
// // // [[Rcpp::export()]]
// // VectorXi sample_h2s_discrete_MH_fast_c(
// //     Map<MatrixXd> UtEta,
// //     Map<VectorXd> tot_Eta_prec,
// //     Map<VectorXd> discrete_priors,
// //     VectorXi h2_index,
// //     Map<MatrixXd> h2s_matrix,
// //     Map<VectorXd> s,
// //     Map<VectorXd> r_draws,
// //     Map<VectorXd> state_draws,
// //     double step_size,
// //     int grainSize
// // ){
// //
// //   int p = UtEta.cols();
// //
// //   struct sampleColumn : public RcppParallel::Worker {
// //     const MatrixXd UtEta;
// //     const MatrixXd h2s_matrix;
// //     const VectorXd s;
// //     const VectorXd tot_Eta_prec;
// //     const VectorXd discrete_priors;
// //     const VectorXd r_draws;
// //     const VectorXd state_draws;
// //     const VectorXi h2_index;
// //     const double step_size;
// //     VectorXi &new_index;
// //
// //     sampleColumn(const MatrixXd UtEta,
// //                  const MatrixXd h2s_matrix,
// //                  const VectorXd s,
// //                  const VectorXd tot_Eta_prec,
// //                  const VectorXd discrete_priors,
// //                  const VectorXd r_draws,
// //                  const VectorXd state_draws,
// //                  const VectorXi h2_index,
// //                  const double step_size,
// //                  VectorXi &new_index):
// //
// //       UtEta(UtEta), h2s_matrix(h2s_matrix),s(s),
// //       tot_Eta_prec(tot_Eta_prec),  discrete_priors(discrete_priors),
// //       r_draws(r_draws),state_draws(state_draws),h2_index(h2_index),step_size(step_size),new_index(new_index) {}
// //
// //     void operator()(std::size_t begin, std::size_t end) {
// //       int n = UtEta.rows();
// //       for(std::size_t j = begin; j < end; j++){
// //         int old_state = h2_index[j] - 1;
// //         double old_h2 = h2s_matrix.coeffRef(0,old_state-1);
// //         double old_log_p = log_prob_h2_fast_c(UtEta.col(j),s,old_h2,n,tot_Eta_prec[j],discrete_priors[old_state]);
// //
// //         VectorXd candidate_new_states = find_candidate_states(h2s_matrix,step_size,old_state);
// //         int r = state_draws[j] * (candidate_new_states.size());
// //         int proposed_state = candidate_new_states[r];
// //
// //         double new_h2 = h2s_matrix.coeffRef(0,proposed_state-1);
// //         double new_log_p = log_prob_h2_fast_c(UtEta.col(j),s,new_h2,n,tot_Eta_prec[j],discrete_priors[proposed_state]);
// //
// //         VectorXd candidate_states_from_new_state = find_candidate_states(h2s_matrix,step_size,proposed_state);
// //
// //         double forward_prob = 1.0 / candidate_new_states.size();
// //         double back_prob = 1.0 / candidate_states_from_new_state.size();
// //
// //         double log_MH_ratio = new_log_p - old_log_p + log(forward_prob) - log(back_prob);
// //
// //         if(log(r_draws[j]) < log_MH_ratio) {
// //           new_index[j] = proposed_state;
// //         } else {
// //           new_index[j] = old_state;
// //         }
// //       }
// //     }
// //   };
// //   VectorXi new_index(p);
// //
// //   sampleColumn sampler(UtEta,h2s_matrix,s,tot_Eta_prec,discrete_priors,r_draws,state_draws,h2_index,step_size,new_index);
// //   RcppParallel::parallelFor(0,p,sampler,grainSize);
// //   return new_index;
// // }
