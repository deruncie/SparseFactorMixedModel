// #include <math.h>
// #include <iostream>
// #include "BSFG_types.h"
// // // #include<Eigen/SparseCholesky>
//
// VectorXd sample_coefs_single(
//     VectorXd UtEta,
//     MatrixXd UtW,
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
//   MatrixXd RinvSqUtW = R_sq_diag.cwiseInverse().asDiagonal() * UtW;
//   VectorXd eta_std = eta_resid.array()/R_sq_diag.array();
//   VectorXd WtURinvy = RinvSqUtW.transpose() * eta_std;
//
//   VectorXd theta_tilda;
//   if(b < n) {
//     MatrixXd C = RinvSqUtW.transpose() * RinvSqUtW;
//     C.diagonal() += prior_prec;
//     theta_tilda = C.llt().solve(WtURinvy);
//   } else{
//     MatrixXd VAi = UtW * prior_prec.cwiseInverse().asDiagonal();
//     MatrixXd inner = VAi*UtW.transpose();
//     for(int i = 0; i < n; i++) {
//       inner(i,i) += (h2 * s(i) + (1.0-h2)) / tot_Eta_prec;
//     }
//     VectorXd VAiWtURinvy = VAi * WtURinvy;
//     VectorXd outerWtURinvy = VAi.transpose() * inner.llt().solve(VAiWtURinvy);
//     theta_tilda = WtURinvy.array() / prior_prec.array();
//     theta_tilda -= outerWtURinvy;
//   }
//
//   VectorXd coefs = theta_tilda + theta_star;
//   return coefs;
// }
//
// // [[Rcpp::export()]]
// MatrixXd sample_randomEffects_parallel_sparse_missing_c_Eigen2 (
//     Map<MatrixXd> Eta,
//     Map<ArrayXd> tot_prec,
//     Map<ArrayXd> h2,
//     List invert_aZZt_Kinv,
//     VectorXi Y_obs_index,
//     Map<ArrayXXd> randn_draws,
//     int grainSize) {
//   //samples genetic effects on factors (F_a) conditional on the factor scores F:
//   // F_i = F_{a_i} + E_i, E_i~N(0,s2*(1-h2)*I) for each latent trait i
//   // U_i = zeros(r,1) if h2_i = 0
//   // it is assumed that s2 = 1 because this scaling factor is absorbed in
//   // Lambda
//   // invert_aZZt_Kinv has parameters to diagonalize a*Z*Z' + b*K for fast
//   // inversion:
//
//   ArrayXd a_prec = tot_prec / h2;
//   ArrayXd e_prec = tot_prec / (1.0-h2);
//
//   std::vector<MSpMat> U_s;
//   std::vector<MSpMat> ZUt_s;
//   std::vector<VectorXd> s1_s;
//   std::vector<VectorXd> s2_s;
//   std::vector<VectorXi> Y_obs;
//   for(int i = 0; i < invert_aZZt_Kinv.length(); i++){
//     Rcpp::List invert_aZZt_Kinv_i = Rcpp::as<Rcpp::List>(invert_aZZt_Kinv[i]);
//     U_s.push_back(as<MSpMat>(invert_aZZt_Kinv_i["U"]));
//     ZUt_s.push_back(as<MSpMat>(invert_aZZt_Kinv_i["ZUt"]));
//     s1_s.push_back(as<VectorXd>(invert_aZZt_Kinv_i["s1"]));
//     s2_s.push_back(as<VectorXd>(invert_aZZt_Kinv_i["s2"]));
//     Y_obs.push_back(as<VectorXi>(invert_aZZt_Kinv_i["Y_obs"]));
//   }
//
//   int p = Eta.cols();
//   int r = randn_draws.rows();
//   MatrixXd b(r,p);
//   for(int j = 0; j < p; j++){
//     VectorXi Y_obs_j = Y_obs[Y_obs_index[j]-1];
//     int n_obs = Y_obs_j.size();
//     VectorXd Eta_j(n_obs);
//     for(int i = 0; i < n_obs; i++){
//       Eta_j[i] = Eta.coeffRef(Y_obs_j[i]-1,j);
//     }
//     b.col(j) = ZUt_s[Y_obs_index[j]-1] * Eta_j * e_prec[j];
//   }
//
//   MatrixXd effects(r,p);
//
//   struct sampleColumn : public Worker {
//     std::vector<MSpMat> U_s;
//     std::vector<VectorXd> s1_s;
//     std::vector<VectorXd> s2_s;
//     VectorXi Y_obs_index;
//     ArrayXd a_prec, e_prec;
//     MatrixXd b;
//     ArrayXXd randn_draws;
//     MatrixXd &effects;
//
//     sampleColumn(std::vector<MSpMat> U_s, std::vector<VectorXd> s1_s, std::vector<VectorXd> s2_s, VectorXi Y_obs_index,
//                  ArrayXd a_prec, ArrayXd e_prec, MatrixXd b, ArrayXXd randn_draws, MatrixXd &effects):
//       U_s(U_s), s1_s(s1_s), s2_s(s2_s), Y_obs_index(Y_obs_index),
//       a_prec(a_prec), e_prec(e_prec), b(b), randn_draws(randn_draws), effects(effects) {}
//
//     void operator()(std::size_t begin, std::size_t end) {
//       for(std::size_t j = begin; j < end; j++){
//         int index = Y_obs_index[j]-1;
//         ArrayXd d = s2_s[index]*a_prec(j) + s1_s[index]*e_prec(j);
//         ArrayXd mlam = b.col(j).array() / d;
//         effects.col(j) = U_s[index] * (mlam + randn_draws.col(j) / sqrt(d)).matrix();
//       }
//     }
//   };
//
//   sampleColumn sampler(U_s,s1_s,s2_s,Y_obs_index, a_prec, e_prec, b, randn_draws, effects);
//   RcppParallel::parallelFor(0,p,sampler,grainSize);
//
//   return(effects);
// }

//
//
// // [[Rcpp::export()]]
// MatrixXd sample_coefs_parallel_sparse_missing_c_Eigen2(
//     Map<MatrixXd> Eta,
//     Map<MatrixXd> W,
//     Map<VectorXd> h2,
//     Map<VectorXd> tot_Eta_prec,
//     Map<MatrixXd> prior_mean,
//     Map<MatrixXd> prior_prec,
//     Map<MatrixXd> randn_theta,
//     Map<VectorXd> randn_e,
//     Rcpp::List invert_aI_bZKZ,
//     VectorXi Y_obs_index,
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
//   std::vector<MSpMat> Ut_s;
//   std::vector<VectorXd> s_s;
//   std::vector<VectorXi> Y_obs;
//   std::vector<MatrixXd> UtW_s;
//   for(int i = 0; i < invert_aI_bZKZ.length(); i++){
//     Rcpp::List invert_aI_bZKZ_i = Rcpp::as<Rcpp::List>(invert_aI_bZKZ[i]);
//     Ut_s.push_back(as<MSpMat>(invert_aI_bZKZ_i["Ut"]));
//     s_s.push_back(as<VectorXd>(invert_aI_bZKZ_i["s"]));
//     Y_obs.push_back(as<VectorXi>(invert_aI_bZKZ_i["Y_obs"]));
//     MatrixXd UtW(Y_obs[i].size(),b);
//     for(int j = 0; j < Y_obs[i].size(); j++) {
//       UtW.row(j) = W.row(Y_obs[i][j]-1);
//     }
//     UtW = Ut_s[i] * UtW;
//     UtW_s.push_back(UtW);
//   }
//
//   std::vector<VectorXd> randn_e_list;
//   int index = 0;
//   for(int i = 0; i < Eta.cols(); i++){
//     int n_obs = Y_obs[Y_obs_index[i]-1].size();
//     randn_e_list.push_back(randn_e.segment(index,n_obs));
//     index += n_obs;
//   }
//
//   struct sampleColumn : public RcppParallel::Worker {
//     MatrixXd Eta;
//     MatrixXd prior_mean, prior_prec, randn_theta;
//     VectorXd h2, tot_Eta_prec;
//     std::vector<MSpMat> Ut_s;
//     std::vector<VectorXd> s_s;
//     std::vector<VectorXi> Y_obs;
//     std::vector<MatrixXd> UtW_s;
//     std::vector<VectorXd> randn_e_list;
//     VectorXi Y_obs_index;
//     int b;
//     MatrixXd &coefs;
//
//     sampleColumn(MatrixXd Eta, MatrixXd prior_mean, MatrixXd prior_prec, VectorXd h2, VectorXd tot_Eta_prec,
//                  MatrixXd randn_theta,
//                  std::vector<MSpMat> Ut_s, std::vector<VectorXd> s_s, std::vector<VectorXi> Y_obs, std::vector<MatrixXd> UtW_s,
//                  std::vector<VectorXd> randn_e_list,
//                  VectorXi Y_obs_index, int b,
//                  MatrixXd &coefs) :
//       Eta(Eta), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta),
//       h2(h2), tot_Eta_prec(tot_Eta_prec),
//       Ut_s(Ut_s), s_s(s_s), Y_obs(Y_obs), UtW_s(UtW_s),
//       randn_e_list(randn_e_list),
//       Y_obs_index(Y_obs_index), b(b),
//       coefs(coefs) {}
//
//     void operator()(std::size_t begin, std::size_t end) {
//       for(std::size_t j = begin; j < end; j++){
//         VectorXi Y_obs_j = Y_obs[Y_obs_index[j]-1];
//         int n_obs = Y_obs_j.size();
//         VectorXd UtEta_j(n_obs);
//         for(int i = 0; i < n_obs; i++){
//           UtEta_j[i] = Eta.coeffRef(Y_obs_j[i]-1,j);
//         }
//         UtEta_j = Ut_s[Y_obs_index[j]-1] * UtEta_j;
//
//         coefs.col(j) = sample_coefs_single(UtEta_j, UtW_s[Y_obs_index[j]-1], prior_mean.col(j), prior_prec.col(j), h2(j),
//                   tot_Eta_prec(j), randn_theta.col(j),randn_e_list[j],s_s[Y_obs_index[j]-1],b,n_obs);
//       }
//     }
//   };
//
//   sampleColumn sampler(Eta, prior_mean, prior_prec, h2, tot_Eta_prec, randn_theta,Ut_s, s_s, Y_obs, UtW_s,randn_e_list,
//                        Y_obs_index, b, coefs);
//   RcppParallel::parallelFor(0,p,sampler,grainSize);
//
//   return(coefs);
// }
//
//
// // [[Rcpp::export()]]
// VectorXd tot_prec_scores_missing_c (
//     Map<MatrixXd> Eta,
//     Map<VectorXd> h2,
//     Rcpp::List invert_aI_bZKZ,
//     VectorXi Y_obs_index
// ) {
//
//   int p = Eta.cols();
//
//   std::vector<MSpMat> Ut_s;
//   std::vector<VectorXd> s_s;
//   std::vector<VectorXi> Y_obs;
//   for(int i = 0; i < invert_aI_bZKZ.length(); i++){
//     Rcpp::List invert_aI_bZKZ_i = Rcpp::as<Rcpp::List>(invert_aI_bZKZ[i]);
//     Ut_s.push_back(as<MSpMat>(invert_aI_bZKZ_i["Ut"]));
//     s_s.push_back(as<VectorXd>(invert_aI_bZKZ_i["s"]));
//     Y_obs.push_back(as<VectorXi>(invert_aI_bZKZ_i["Y_obs"]));
//   }
//
//   VectorXd scores(p);
//
//   for(int j = 0; j < p; j++){
//     ArrayXd s = s_s[Y_obs_index[j]-1];
//     VectorXi Y_obs_j = Y_obs[Y_obs_index[j]-1];
//     int n_obs = Y_obs_j.size();
//     VectorXd UtEta_j(n_obs);
//     for(int i = 0; i < n_obs; i++){
//       UtEta_j[i] = Eta.coeffRef(Y_obs_j[i]-1,j);
//     }
//     UtEta_j = Ut_s[Y_obs_index[j]-1] * UtEta_j;
//     ArrayXd Sigma_sqrt = sqrt(h2(j) * s + (1.0 - h2(j)));
//     VectorXd SiUtEta_j = UtEta_j.array() / Sigma_sqrt;
//     scores(j) = SiUtEta_j.dot(SiUtEta_j);
//   }
//   return scores;
// }
//
// // [[Rcpp::export()]]
// MatrixXd log_p_h2s_fast_missing(
//     Map<MatrixXd> Eta,
//     Map<VectorXd> tot_Eta_prec,
//     Map<VectorXd> discrete_priors,
//     Rcpp::List invert_aI_bZKZ,
//     VectorXi Y_obs_index,
//     int grainSize)
// {
//
//   int b = discrete_priors.size();
//   int p = Eta.cols();
//
//   std::vector<MSpMat> Ut_s;
//   std::vector<VectorXd> s_s;
//   std::vector<VectorXi> Y_obs;
//   for(int i = 0; i < invert_aI_bZKZ.length(); i++){
//     Rcpp::List invert_aI_bZKZ_i = Rcpp::as<Rcpp::List>(invert_aI_bZKZ[i]);
//     Ut_s.push_back(as<MSpMat>(invert_aI_bZKZ_i["Ut"]));
//     s_s.push_back(as<VectorXd>(invert_aI_bZKZ_i["s"]));
//     Y_obs.push_back(as<VectorXi>(invert_aI_bZKZ_i["Y_obs"]));
//   }
//
//   MatrixXd log_ps(b,p);
//
//   std::vector<VectorXd> std_scores_b2;
//   for(int j = 0; j < p; j++) {
//     VectorXi Y_obs_j = Y_obs[Y_obs_index[j]-1];
//     int n_obs = Y_obs_j.size();
//     VectorXd Eta_std_j(n_obs);
//     for(int i = 0; i < n_obs; i++){
//       Eta_std_j[i] = Eta.coeffRef(Y_obs_j[i]-1,j);
//     }
//     VectorXd UtEta_std_j = Ut_s[Y_obs_index[j]-1] * Eta_std_j * sqrt(tot_Eta_prec[j]);
//     std_scores_b2.push_back(UtEta_std_j.cwiseProduct(UtEta_std_j));
//   }
//
//
//   struct sampleColumn : public Worker {
//     std::vector<VectorXd> std_scores_b2;
//     VectorXd discrete_priors;
//     VectorXd tot_Eta_prec;
//     std::vector<MSpMat> Ut_s;
//     std::vector<VectorXd> s_s;
//     std::vector<VectorXi> Y_obs;
//     VectorXi Y_obs_index;
//     MatrixXd &log_ps;
//
//     sampleColumn(std::vector<VectorXd> std_scores_b2,
//                  VectorXd discrete_priors,
//                  VectorXd tot_Eta_prec,
//                  std::vector<MSpMat> Ut_s,
//                  std::vector<VectorXd> s_s,
//                  std::vector<VectorXi> Y_obs,
//                  VectorXi Y_obs_index,
//                  MatrixXd &log_ps):
//       std_scores_b2(std_scores_b2), discrete_priors(discrete_priors), tot_Eta_prec(tot_Eta_prec),Ut_s(Ut_s), s_s(s_s), Y_obs(Y_obs),Y_obs_index(Y_obs_index), log_ps(log_ps) {}
//
//     void operator()(std::size_t begin, std::size_t end) {
//       int b = discrete_priors.size();
//       int p = std_scores_b2.size();
//       for(std::size_t i = begin; i < end; i++){
//         double h2 = double(i)/b;
//         for(int j = 0; j < p; j++) {
//           int index = Y_obs_index[j] - 1;
//           int n_obs = Y_obs[index].size();
//           VectorXd s2s = h2*s_s[index].array() + (1.0-h2);
//           VectorXd std_scores_b2_j = -0.5 * std_scores_b2[j].cwiseProduct(s2s.cwiseInverse());
//           double det = -n_obs/2 * log(2.0*M_PI) - 0.5*s2s.array().log().sum();
//           log_ps.coeffRef(i,j) = std_scores_b2_j.sum() + det + log(discrete_priors[i]) + n_obs/2*log(tot_Eta_prec[j]);
//         }
//       }
//     }
//   };
//
//   sampleColumn sampler(std_scores_b2,discrete_priors,tot_Eta_prec,Ut_s,s_s,Y_obs,Y_obs_index,log_ps);
//   RcppParallel::parallelFor(0,b,sampler,grainSize);
//   return(log_ps);
// }
//
//

//
//
//
// // [[Rcpp::export()]]
// MatrixXd sample_factors_scores_sparse_mising_c_Eigen(
//     Map<MatrixXd> Eta_tilde,
//     Map<MatrixXd> prior_mean,
//     Map<MatrixXd> Lambda,
//     Map<VectorXd> resid_Eta_prec,
//     Map<VectorXd> F_e_prec,
//     Map<MatrixXd> randn_draws,
//     VectorXi Y_row_obs_index,
//     Rcpp::List Y_row_obs_sets
// ) {
//   //Sample factor scores given factor loadings (F_a), factor residual variances (F_e_prec) and
//   //phenotype residuals
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
//   for(int i = 0; i < Y_row_obs_sets.length(); i++){
//     VectorXi obs_set_i = as<VectorXi>(Y_row_obs_sets[i]);
//     obs_sets.push_back(obs_set_i);
//     int n_traits = obs_set_i.size();
//     MatrixXd Lmsg_i(n_traits,k);
//     MatrixXd Lambda_i(n_traits,k);
//     for(int j = 0; j < n_traits; j++){
//       int index = obs_set_i[j]-1;
//       Lmsg_i.row(j) = Lmsg.row(index);
//       Lambda_i.row(j) = Lambda.row(index);
//     }
//     Lmsg_s.push_back(Lmsg_i);
//     MatrixXd Sigma = Lambda_i.transpose() * Lmsg_i;
//     Sigma.diagonal() += F_e_prec;
//     Eigen::LLT<MatrixXd> chol_Sigma;
//     chol_Sigma.compute(Sigma);
//     MatrixXd R = chol_Sigma.matrixU();
//     R_s.push_back(R);
//   }
//
//   for(int i = 0; i < n; i++){
//     int index = Y_row_obs_index[i]-1;
//     VectorXi obs_set_i = obs_sets[index];
//     int n_traits = obs_set_i.size();
//     MatrixXd R_i = R_s[index];
//     MatrixXd Lmsg_i = Lmsg_s[index];
//     Eigen::RowVectorXd Eta_i(n_traits);
//     for(int j = 0; j < n_traits; j++){
//       int trait_index = obs_set_i[j]-1;
//       Eta_i[j] = Eta_tilde.coeffRef(i,trait_index);
//     }
//     VectorXd Meta = R_i.transpose().triangularView<Lower>().solve((Eta_i * Lmsg_i + prior_mean.row(i) * F_e_prec.asDiagonal()).transpose());
//     Ft.col(i) = R_i.triangularView<Upper>().solve(Meta + randn_draws.col(i));
//   }
//   return Ft.transpose();
// }
//
//
// //
// // // [[Rcpp::export()]]
// // MatrixXd add_Sp_Mat(SpMat X1,Map<MatrixXd> X2){
// //   return(X1+X2);
// // }
// // // [[Rcpp::export()]]
// // SpMat add_Sp_Sp(SpMat X1,SpMat X2){
// //   return(X1+X2);
// // }
// // // [[Rcpp::export()]]
// // MatrixXd to_Mat(SpMat X1){
// //   return(X1.toDense());
// // }
//
// //
// //
// //
// //
// //
// //
// // VectorXd sample_coefs_single_hierarchical(
// //     VectorXd UtEta,
// //     SpMat UtW,
// //     SpMat UtWX,
// //     SpMat X,
// //     VectorXd prior_mean,
// //     VectorXd prior_prec,
// //     double h2,
// //     double tot_Eta_prec,
// //     VectorXd randn_theta,
// //     VectorXd randn_e,
// //     VectorXd s,
// //     int b,
// //     int n,
// //     int r
// // ) {
// //
// //   VectorXd R_sq_diag = ((h2 * s.array() + (1.0-h2))/tot_Eta_prec).sqrt();
// //   VectorXd theta_star = randn_theta.array()/prior_prec.array().sqrt();
// //   theta_star += prior_mean;
// //   VectorXd e_star = randn_e.array() * R_sq_diag.array();
// //   MatrixXd UtWX_theta_star = UtWX*theta_star;
// //   VectorXd eta_resid = UtEta - UtWX_theta_star - e_star;
// //   SpMat RinvSqUtWX = R_sq_diag.cwiseInverse().asDiagonal() * UtWX;
// //   VectorXd eta_std = eta_resid.array()/R_sq_diag.array();
// //   VectorXd XtWtURinvy = RinvSqUtWX.transpose() * eta_std;
// //
// //   VectorXd theta_tilda;
// //   if(b < r) {
// //     MatrixXd C = RinvSqUtWX.transpose() * RinvSqUtWX;
// //     C.diagonal() = C.diagonal() + prior_prec;
// //     theta_tilda = C.householderQr().solve(XtWtURinvy);
// //   } else{
// //     SpMat B = UtW.transpose() * R_sq_diag.array().pow(2).inverse().matrix().asDiagonal() * UtW;
// //     SimplicialLLT<SparseMatrix<double> >solver;
// //     SpMat I(r,r);
// //     I.setIdentity();
// //     SpMat Bi = solver.compute(B).solve(I);
// //     MatrixXd VAi = X * prior_prec.cwiseInverse().asDiagonal();
// //     MatrixXd inner = VAi * X.transpose() + Bi;
// //     VectorXd VAiXtWtURinvy = VAi * XtWtURinvy;
// //     VectorXd outerXtWtURinvy = VAi.transpose() * inner.householderQr().solve(VAiXtWtURinvy);
// //     theta_tilda = XtWtURinvy.array() / prior_prec.array();
// //     theta_tilda -= outerXtWtURinvy;
// //   }
// //   VectorXd coefs = theta_tilda + theta_star;
// //   return coefs;
// // }
// //
// // // [[Rcpp::export()]]
// // MatrixXd sample_coefs_hierarchical_parallel_sparse_c_Eigen(
// //     MSpMat Ut,
// //     Map<MatrixXd> Eta,
// //     MSpMat W,
// //     MSpMat X,
// //     Map<VectorXd> h2,
// //     Map<VectorXd> tot_Eta_prec,
// //     Map<VectorXd> s,
// //     Map<MatrixXd> prior_mean,
// //     Map<MatrixXd> prior_prec,
// //     Map<MatrixXd> randn_theta,
// //     Map<MatrixXd> randn_e,
// //     int grainSize){
// //
// //   // Sample regression coefficients
// //   // columns of matrices are independent
// //   // each column conditional posterior is a MVN due to conjugacy
// //
// //
// //   MatrixXd UtEta = Ut*Eta;
// //   SpMat UtW = Ut*W;
// //   SpMat UtWX = UtW*X;
// //
// //   int p = UtEta.cols();
// //   int b = X.cols();
// //   int n = UtW.rows();
// //   int r = X.rows();
// //
// //   MatrixXd coefs(b,p);
// //
// //   struct sampleColumn : public RcppParallel::Worker {
// //     SpMat UtW,UtWX,X;
// //     MatrixXd UtEta;
// //     MatrixXd prior_mean, prior_prec, randn_theta, randn_e;
// //     VectorXd h2, tot_Eta_prec, s;
// //     int b,n,r;
// //     MatrixXd &coefs;
// //
// //     sampleColumn(SpMat UtW, SpMat UtWX, SpMat X,MatrixXd UtEta, MatrixXd prior_mean, MatrixXd prior_prec, VectorXd h2, VectorXd tot_Eta_prec,
// //                  MatrixXd randn_theta, MatrixXd randn_e,
// //                  VectorXd s, int b, int n,int r,
// //                  MatrixXd &coefs) :
// //       UtW(UtW), UtWX(UtWX), X(X),UtEta(UtEta), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
// //       h2(h2), tot_Eta_prec(tot_Eta_prec),
// //       s(s), b(b), n(n),r(r),
// //       coefs(coefs) {}
// //
// //     void operator()(std::size_t begin, std::size_t end) {
// //       for(std::size_t j = begin; j < end; j++){
// //       // for(int j = 0; j < 1; j++){
// //         coefs.col(j) = sample_coefs_single_hierarchical(UtEta.col(j), UtW, UtWX, X,prior_mean.col(j), prior_prec.col(j), h2(j), tot_Eta_prec(j), randn_theta.col(j),randn_e.col(j),s,b,n,r);
// //       }
// //     }
// //   };
// //
// //   sampleColumn sampler(UtW, UtWX, X,UtEta, prior_mean, prior_prec, h2, tot_Eta_prec, randn_theta,randn_e, s, b, n, r,coefs);
// //   RcppParallel::parallelFor(0,p,sampler,grainSize);
// //
// //   return(coefs);
// // }
// //
// //
// //
// // //
// // //
// // // // -------------------------------------------- //
// // // // --------------- sample h2s MH fast --------- //
// // // // -------------------------------------------- //
// // //
// // // // [[Rcpp::export()]]
// // // double log_prob_h2_fast_c(
// // //     ArrayXd Uty,
// // //     ArrayXd s,
// // //     double h2,
// // //     int n,
// // //     double tot_Eta_prec,
// // //     double discrete_prior
// // // ){
// // //   ArrayXd s2s = h2*s + (1.0-h2);
// // //   double score2 = (Uty.pow(2)/s2s).sum() * tot_Eta_prec;
// // //   double log_det_Sigma = s2s.log().sum();
// // //
// // //   double log_p = -n/2.0 * log(2*M_PI) - 0.5*(log_det_Sigma - n*log(tot_Eta_prec)) - 0.5 * score2 + log(discrete_prior);
// // //   return log_p;
// // // }
// // //
// // //
// // // // [[Rcpp::export()]]
// // // VectorXd find_candidate_states(
// // //     MatrixXd h2s_matrix,
// // //     double step_size,
// // //     int old_state
// // // ) {
// // //   VectorXd dists = (h2s_matrix.colwise() - h2s_matrix.col(old_state)).cwiseAbs().colwise().sum();
// // //   VectorXd indices(dists.size());
// // //   int count = 0;
// // //   for(int i = 0; i < dists.size(); i++){
// // //     if(dists[i] < step_size & dists[i] > 0) {
// // //       indices[count] = i;
// // //       count++;
// // //     }
// // //   }
// // //   if(count == 0) {  // return all indices as candidates
// // //     for(int i = 0; i < dists.size(); i++){
// // //       indices[count] = i;
// // //       count++;
// // //     }
// // //   }
// // //   return indices.head(count);
// // // }
// // //
// // // // [[Rcpp::export()]]
// // // VectorXi sample_h2s_discrete_MH_fast_c(
// // //     Map<MatrixXd> UtEta,
// // //     Map<VectorXd> tot_Eta_prec,
// // //     Map<VectorXd> discrete_priors,
// // //     VectorXi h2_index,
// // //     Map<MatrixXd> h2s_matrix,
// // //     Map<VectorXd> s,
// // //     Map<VectorXd> r_draws,
// // //     Map<VectorXd> state_draws,
// // //     double step_size,
// // //     int grainSize
// // // ){
// // //
// // //   int p = UtEta.cols();
// // //
// // //   struct sampleColumn : public RcppParallel::Worker {
// // //     const MatrixXd UtEta;
// // //     const MatrixXd h2s_matrix;
// // //     const VectorXd s;
// // //     const VectorXd tot_Eta_prec;
// // //     const VectorXd discrete_priors;
// // //     const VectorXd r_draws;
// // //     const VectorXd state_draws;
// // //     const VectorXi h2_index;
// // //     const double step_size;
// // //     VectorXi &new_index;
// // //
// // //     sampleColumn(const MatrixXd UtEta,
// // //                  const MatrixXd h2s_matrix,
// // //                  const VectorXd s,
// // //                  const VectorXd tot_Eta_prec,
// // //                  const VectorXd discrete_priors,
// // //                  const VectorXd r_draws,
// // //                  const VectorXd state_draws,
// // //                  const VectorXi h2_index,
// // //                  const double step_size,
// // //                  VectorXi &new_index):
// // //
// // //       UtEta(UtEta), h2s_matrix(h2s_matrix),s(s),
// // //       tot_Eta_prec(tot_Eta_prec),  discrete_priors(discrete_priors),
// // //       r_draws(r_draws),state_draws(state_draws),h2_index(h2_index),step_size(step_size),new_index(new_index) {}
// // //
// // //     void operator()(std::size_t begin, std::size_t end) {
// // //       int n = UtEta.rows();
// // //       for(std::size_t j = begin; j < end; j++){
// // //         int old_state = h2_index[j] - 1;
// // //         double old_h2 = h2s_matrix.coeffRef(0,old_state-1);
// // //         double old_log_p = log_prob_h2_fast_c(UtEta.col(j),s,old_h2,n,tot_Eta_prec[j],discrete_priors[old_state]);
// // //
// // //         VectorXd candidate_new_states = find_candidate_states(h2s_matrix,step_size,old_state);
// // //         int r = state_draws[j] * (candidate_new_states.size());
// // //         int proposed_state = candidate_new_states[r];
// // //
// // //         double new_h2 = h2s_matrix.coeffRef(0,proposed_state-1);
// // //         double new_log_p = log_prob_h2_fast_c(UtEta.col(j),s,new_h2,n,tot_Eta_prec[j],discrete_priors[proposed_state]);
// // //
// // //         VectorXd candidate_states_from_new_state = find_candidate_states(h2s_matrix,step_size,proposed_state);
// // //
// // //         double forward_prob = 1.0 / candidate_new_states.size();
// // //         double back_prob = 1.0 / candidate_states_from_new_state.size();
// // //
// // //         double log_MH_ratio = new_log_p - old_log_p + log(forward_prob) - log(back_prob);
// // //
// // //         if(log(r_draws[j]) < log_MH_ratio) {
// // //           new_index[j] = proposed_state;
// // //         } else {
// // //           new_index[j] = old_state;
// // //         }
// // //       }
// // //     }
// // //   };
// // //   VectorXi new_index(p);
// // //
// // //   sampleColumn sampler(UtEta,h2s_matrix,s,tot_Eta_prec,discrete_priors,r_draws,state_draws,h2_index,step_size,new_index);
// // //   RcppParallel::parallelFor(0,p,sampler,grainSize);
// // //   return new_index;
// // // }

// // [[Rcpp::export()]]
// MatrixXd test_mult1(MatrixXd X, MatrixXd Y){
//   MatrixXd res = X * Y;
//   return res;
// }
// // [[Rcpp::export()]]
// MatrixXd test_mult1b(Map<MatrixXd> X, Map<MatrixXd> Y){
//   MatrixXd res = X * Y;
//   return res;
// }
// // [[Rcpp::export()]]
// MatrixXd test_mult2(MatrixXd X, MatrixXd Y,IntegerVector rows){
//   MatrixXd res(rows.length(),Y.cols());
//   for(int i = 0; i < rows.length(); i++){
//     res.row(i) = X.row(rows(i)-1) * Y;
//   }
//   return res;
// }
// // [[Rcpp::export()]]
// MatrixXd test_mult3(MatrixXd X, MatrixXd Y,IntegerVector rows){
//   MatrixXd res(Y.cols(),rows.length());
//   for(int i = 0; i < rows.length(); i++){
//     res.col(i) = (X.row(rows(i)-1) * Y).transpose();
//   }
//   return res.transpose();
// }
//
// // [[Rcpp::export()]]
// MatrixXd test_mult4(Map<MatrixXd> X, Map<MatrixXd> Y,IntegerVector rows){
//   Eigen::Matrix<double,Dynamic,Dynamic,RowMajor> res(rows.length(),Y.cols());
//   for(int i = 0; i < rows.length(); i++){
//     res.row(i) = X.row(rows(i)) * Y;
//   }
//   return res;
// }
// MatrixXd subset_rows(MatrixXd X, IntegerVector rows){
//   Eigen::Matrix<double,Dynamic,Dynamic,RowMajor> X_sub(rows.length(),X.cols());
//   for(int i = 0; i < rows.length(); i++){
//     X_sub.row(i) = X.row(rows(i)-1);
//   }
//   return X_sub;
// }
// // [[Rcpp::export()]]
// MatrixXd test_mult5(Map<MatrixXd> X, Map<MatrixXd> Y,IntegerVector rows){
//   MatrixXd X_sub = subset_rows(X,rows);
//   MatrixXd res = X_sub * Y;
//   return res;
// }
//
// // [[Rcpp::export()]]
// MatrixXd test_mult5b(Map<MatrixXd> X, Map<MatrixXd> Y,IntegerVector rows){
//   Eigen::Matrix<double,Dynamic,Dynamic,RowMajor> X_sub(rows.length(),X.cols());
//   for(int i = 0; i < rows.length(); i++){
//     X_sub.row(i) = X.row(rows(i)-1);
//   }
//   MatrixXd res = X_sub * Y;
//   return res;
// }
