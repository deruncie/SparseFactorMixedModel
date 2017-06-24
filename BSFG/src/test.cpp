// #include <math.h>
// #include <iostream>
// #include "BSFG_types.h"
//
// // [[Rcpp::depends(RcppEigen)]]
// using namespace Eigen;
// using namespace RcppParallel;
//
// VectorXd sample_MME_single_diagK(
//     VectorXd y,
//     SpMat W,
//     VectorXd prior_mean,
//     VectorXd prior_prec,
//     SpMat chol_R,
//     double tot_Eta_prec,
//     VectorXd randn_theta,
//     VectorXd randn_e
// ){
//
//   VectorXd theta_star = randn_theta.array() / prior_prec.cwiseSqrt().array();
//   theta_star += prior_mean;
//   VectorXd e_star = chol_R * randn_e / sqrt(tot_Eta_prec);
//   MatrixXd W_theta_star = W * theta_star;
//   VectorXd y_resid = y - W_theta_star - e_star;
//
//   MatrixXd W_mat = W;
//   MatrixXd RinvSqW = chol_R.transpose().triangularView<Lower>().solve(W_mat);
//   VectorXd WtRinvy = RinvSqW.transpose() * chol_R.transpose().triangularView<Lower>().solve(y_resid) * tot_Eta_prec;
//
//   VectorXd theta_tilda;
//
//   if(W.cols() < W.rows()) {
//     MatrixXd C = RinvSqW.transpose() * RinvSqW * tot_Eta_prec;
//     C.diagonal() = C.diagonal() + prior_prec;
//     theta_tilda = C.householderQr().solve(WtRinvy);
//   } else{
//     MatrixXd R = chol_R.transpose() * chol_R / tot_Eta_prec;
//     MatrixXd AiU = (W * prior_prec.cwiseInverse().asDiagonal()).transpose();
//     MatrixXd R_VAiU = R + W * AiU;
//     MatrixXd inner = AiU * R_VAiU.householderQr().solve(AiU.transpose());
//     theta_tilda = WtRinvy.array()/prior_prec.array();
//     theta_tilda -= inner * WtRinvy;
//   }
//
//   VectorXd theta = theta_star + theta_tilda;
//
//   return theta;
// }
//
// // [[Rcpp::export()]]
// Rcpp::List sample_MME_fixedEffects_cis_c(
//     Map<MatrixXd> Y,
//     Map<MatrixXd> W,
//     Rcpp::List cis_genotypes,
//     Rcpp::List Sigma_Choleskys,
//     Rcpp::IntegerVector h2s_index,
//     Map<VectorXd> tot_Eta_prec,
//     Map<MatrixXd> prior_mean,
//     Map<MatrixXd> prior_prec,
//     Map<MatrixXd> randn_theta,
//     Map<MatrixXd> randn_e,
//     Map<VectorXd> randn_cis,
//     Map<VectorXd> cis_effect_index,
//     int grainSize) {
//
//   int b = randn_theta.rows();
//   int p = randn_theta.cols();
//
//   std::vector<MSpMat> chol_R_list;
//   for(int i = 0; i < max(h2s_index); i++){
//     Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);
//     chol_R_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
//   }
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
//     MatrixXd Y;
//     MatrixXd W;
//     std::vector<MatrixXd> cis_X;
//     MatrixXd prior_mean, prior_prec, randn_theta, randn_e;
//     VectorXd randn_cis;
//     const std::vector<MSpMat> chol_R_list;
//     RVector<int> h2s_index;
//     VectorXd cis_effect_index;
//     VectorXd tot_Eta_prec;
//     MatrixXd &coefs;
//     VectorXd &cis_effects;
//
//     sampleColumn(MatrixXd Y,
//                  MatrixXd W,
//                  std::vector<MatrixXd> cis_X,
//                  MatrixXd prior_mean,
//                  MatrixXd prior_prec,
//                  const std::vector<MSpMat> chol_R_list,
//                  const Rcpp::IntegerVector h2s_index,
//                  VectorXd cis_effect_index,
//                  VectorXd tot_Eta_prec,
//                  MatrixXd randn_theta,
//                  MatrixXd randn_e,
//                  VectorXd randn_cis,
//                  MatrixXd &coefs,
//                  VectorXd &cis_effects):
//       Y(Y), W(W), cis_X(cis_X),prior_mean(prior_mean), prior_prec(prior_prec),
//       randn_theta(randn_theta), randn_e(randn_e), randn_cis(randn_cis),
//       chol_R_list(chol_R_list), h2s_index(h2s_index), cis_effect_index(cis_effect_index),tot_Eta_prec(tot_Eta_prec),
//       coefs(coefs), cis_effects(cis_effects) {}
//
//     void operator()(std::size_t begin, std::size_t end) {
//       int n = W.rows();
//       int b = W.cols();
//       for(std::size_t j = begin; j < end; j++){
//         int h2_index = h2s_index[j] - 1;
//         SpMat chol_R = chol_R_list[h2_index];
//
//         int b_cis = cis_X[j].cols();
//         MatrixXd W_cisj(n,b+b_cis);
//         W_cisj << W, cis_X[j];
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
//         VectorXd result = sample_MME_single_diagK(Y.col(j), W_cisj.sparseView(), prior_mean_j, prior_prec_j, chol_R, tot_Eta_prec[j], randn_theta_j,randn_e.col(j));
//         coefs.col(j) = result.head(b);
//         cis_effects.segment(cis_effect_index[j],b_cis) = result.tail(b_cis);
//       }
//     }
//   };
//
//   sampleColumn sampler(Y,W,cis_X,prior_mean,prior_prec,chol_R_list,h2s_index,cis_effect_index,tot_Eta_prec,randn_theta,randn_e, randn_cis,coefs,cis_effects);
//   RcppParallel::parallelFor(0,p,sampler,grainSize);
//
//   return(Rcpp::List::create(coefs,cis_effects));
// }
