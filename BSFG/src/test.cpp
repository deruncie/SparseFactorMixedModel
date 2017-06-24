// #include <math.h>
// #include <iostream>
// #include "BSFG_types.h"
//
// // [[Rcpp::depends(RcppEigen)]]
// using namespace Eigen;
// using namespace RcppParallel;
//
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
//     C.diagonal() = C.diagonal() + prior_prec;
//     theta_tilda = C.householderQr().solve(WtURinvy);
//   } else{
//     MatrixXd VAi = UtW * prior_prec.cwiseInverse().asDiagonal();
//     MatrixXd inner = VAi*UtW.transpose();
//     for(int i = 0; i < n; i++) {
//       inner(i,i) += (h2 * s(i) + (1.0-h2)) / tot_Eta_prec;
//     }
//     VectorXd VAiWtURinvy = VAi * WtURinvy;
//     VectorXd outerWtURinvy = VAi.transpose() * inner.householderQr().solve(VAiWtURinvy);
//     theta_tilda = WtURinvy.array() / prior_prec.array();
//     theta_tilda -= outerWtURinvy;
//   }
//
//   VectorXd coefs = theta_tilda + theta_star;
//   return coefs;
// }
//
// // [[Rcpp::export()]]
// Rcpp::List sample_cis_coefs_parallel_sparse_c_Eigen(
//     MSpMat Ut,
//     Map<MatrixXd> Eta,
//     Map<MatrixXd> W,
//     Rcpp::List cis_genotypes,
//     Map<VectorXd> h2,
//     Map<VectorXd> tot_Eta_prec,
//     Map<VectorXd> s,
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
//   MatrixXd UtEta = Ut*Eta;
//   MatrixXd UtW = Ut*W;
//
//   int p = UtEta.cols();
//   int b = UtW.cols();
//   int n = UtW.rows();
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
//     MatrixXd UtW, UtEta;
//     SpMat Ut;
//     std::vector<MatrixXd> cis_X;
//     MatrixXd prior_mean, prior_prec, randn_theta, randn_e;
//     VectorXd randn_cis;
//     VectorXd h2, tot_Eta_prec, cis_effect_index,s;
//     int b,n;
//     MatrixXd &coefs;
//     VectorXd &cis_effects;
//
//     sampleColumn(MatrixXd UtW, MatrixXd UtEta, SpMat Ut, std::vector<MatrixXd> cis_X,
//                  MatrixXd prior_mean, MatrixXd prior_prec, VectorXd h2,
//                  VectorXd tot_Eta_prec, VectorXd cis_effect_index,
//                  MatrixXd randn_theta, MatrixXd randn_e,VectorXd randn_cis,
//                  VectorXd s, int b, int n,
//                  MatrixXd &coefs, VectorXd &cis_effects) :
//       UtW(UtW), UtEta(UtEta),
//       Ut(Ut), cis_X(cis_X),
//       prior_mean(prior_mean), prior_prec(prior_prec),
//       randn_theta(randn_theta), randn_e(randn_e), randn_cis(randn_cis),
//       h2(h2), tot_Eta_prec(tot_Eta_prec), cis_effect_index(cis_effect_index),
//       s(s), b(b), n(n),
//       coefs(coefs), cis_effects(cis_effects) {}
//
//     void operator()(std::size_t begin, std::size_t end) {
//       for(std::size_t j = begin; j < end; j++){
//         int b_cis = cis_X[j].cols();
//         MatrixXd UtW_cisj(n,b+b_cis);
//         UtW_cisj << UtW, Ut*cis_X[j];
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
//         VectorXd result = sample_coefs_single(UtEta.col(j), UtW_cisj, prior_mean_j, prior_prec_j, h2(j), tot_Eta_prec(j), randn_theta_j,randn_e.col(j),s,b,n);
//         coefs.col(j) = result.head(b);
//         cis_effects.segment(cis_effect_index[j],b_cis) = result.tail(b_cis);
//       }
//     }
//   };
//
//   sampleColumn sampler(UtW, UtEta, Ut, cis_X, prior_mean, prior_prec, h2, tot_Eta_prec, cis_effect_index, randn_theta,randn_e, randn_cis, s, b, n, coefs,cis_effects);
//   RcppParallel::parallelFor(0,p,sampler,grainSize);
//
//   return(Rcpp::List::create(coefs,cis_effects));
// }




//
// // [[Rcpp::export()]]
// List LDLt(MSpMat A,double tol) {
//   Eigen::SimplicialLDLT<SpMat> chol_A;
//   chol_A.compute(A);
//   SpMat Pinv = chol_A.permutationPinv().toDenseMatrix().cast<double>().sparseView();
//   SpMat PiL = Pinv * chol_A.matrixL();
//   return(List::create(
//                       Named("PiL") =  PiL,
//                       Named("d") = chol_A.vectorD()));
// }


// //
// // // [[Rcpp::export()]]
// // MatrixXd sample_randomEffects_parallel_sparse_c_Eigen2 (
// //     Map<MatrixXd> Eta,
// //     MSpMat Z,
// //     Map<ArrayXd> tot_prec,
// //     Map<ArrayXd> h2,
// //     List invert_aZZt_Kinv,
// //     Map<ArrayXXd> randn_draws,
// //     int grainSize) {
// //   //samples genetic effects on factors (F_a) conditional on the factor scores F:
// //   // F_i = F_{a_i} + E_i, E_i~N(0,s2*(1-h2)*I) for each latent trait i
// //   // U_i = zeros(r,1) if h2_i = 0
// //   // it is assumed that s2 = 1 because this scaling factor is absorbed in
// //   // Lambda
// //   // invert_aZZt_Kinv has parameters to diagonalize a*Z*Z' + b*K for fast
// //   // inversion:
// //
// //   ArrayXd a_prec = tot_prec / h2;
// //   ArrayXd e_prec = tot_prec / (1.0-h2);
// //
// //   SpMat U = as<MSpMat>(invert_aZZt_Kinv["U"]);
// //   VectorXd s1 = as<VectorXd>(invert_aZZt_Kinv["s1"]);
// //   VectorXd s2 = as<VectorXd>(invert_aZZt_Kinv["s2"]);
// //
// //   int p = Eta.cols();
// //   int r = Z.cols();
// //   // MatrixXd b = U.transpose() * Z.transpose() * (Eta * e_prec.matrix().asDiagonal());
// //   MatrixXd b = (Z*U).transpose() * (Eta * e_prec.matrix().asDiagonal());
// //
// //   // MatrixXd z = randn(r,p);
// //
// //   MatrixXd effects(r,p);
// //   Rcout << r << " " << p << std::endl;
// //   Rcout << a_prec(p-1) << std::endl;
// //   Rcout << e_prec(p-1) << std::endl;
// //   Rcout << b.rows() << " " << b.cols() << std::endl;
// //
// //   struct sampleColumn : public Worker {
// //     VectorXd s1, s2;
// //     ArrayXd a_prec, e_prec;
// //     SpMat U;
// //     MatrixXd b;
// //     ArrayXXd randn_draws;
// //     MatrixXd &effects;
// //
// //     sampleColumn(VectorXd s1, VectorXd s2, ArrayXd a_prec, ArrayXd e_prec, SpMat U, MatrixXd b, ArrayXXd randn_draws, MatrixXd &effects)
// //       : s1(s1), s2(s2), a_prec(a_prec), e_prec(e_prec), U(U), b(b), randn_draws(randn_draws), effects(effects) {}
// //
// //     void operator()(std::size_t begin, std::size_t end) {
// //       for(std::size_t j = begin; j < end; j++){
// //         ArrayXd d = s2*a_prec(j) + s1*e_prec(j);
// //         ArrayXd mlam = b.col(j).array() / d;
// //         effects.col(j) = U * (mlam + randn_draws.col(j) / sqrt(d)).matrix();
// //       }
// //     }
// //   };
// //   for(int j = 0; j < p; j++){
// //     ArrayXd d = s2*a_prec(j) + s1*e_prec(j);
// //     Rcout << j << d.size() <<  std::endl;
// //     ArrayXd mlam = b.col(j).array() / d;
// //     Rcout << mlam + randn_draws.col(j) / sqrt(d) << std::endl;
// //     effects.col(j) = U * (mlam + randn_draws.col(j) / sqrt(d)).matrix();
// //   }
// //
// //   sampleColumn sampler(s1, s2, a_prec, e_prec, U, b, randn_draws, effects);
// //   RcppParallel::parallelFor(0,p,sampler,grainSize);
// //
// //   return(effects);
// // }

// // [[Rcpp::export()]]
// VectorXd tot_prec_scores_withX_c2 (
//     Map<MatrixXd> UtEta,
//     Map<MatrixXd> B_F,
//     Map<VectorXd> h2,
//     Map<VectorXd> s,
//     Map<MatrixXd> prec_B_F
// ) {

//   int p = UtEta.cols();

//   VectorXd scores(p);

//   for(int i = 0; i < p; i++){
//     ArrayXd Sigma_sqrt = sqrt(h2(i) * s.array() + (1.0 - h2(i)));
//     VectorXd SiUtEta_i = UtEta.col(i).array() / Sigma_sqrt;
//     VectorXd b_std = B_F.col(i).cwiseProduct(prec_B_F.cwiseSqrt());
//     scores(i) = SiUtEta_i.dot(SiUtEta_i) + b_std.dot(b_std);
//   }
//   return scores;
// }
