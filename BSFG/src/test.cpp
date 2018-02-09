// #include <math.h>
// #include <iostream>
// #include "BSFG_types.h"
//
// // [[Rcpp::depends(RcppEigen)]]
// using namespace Eigen;
// using namespace RcppParallel;
//
// // [[Rcpp::export()]]
// MatrixXd rstdnorm_mat2(int n,int p) {  // returns nxp matrix
//   VectorXd X_vec = as<VectorXd>(rnorm(n*p));
//   Map<MatrixXd> X_mat(X_vec.data(),n,p);
//   return(X_mat);
// }
//
// VectorXd sample_MME_single_hierarchical_diagK2(  // returns b x 1 vector
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
//   // Using algorithm from Bhattacharya et al 2016 Biometrika. https://academic.oup.com/biomet/article/103/4/985/2447851
//   // Phi = sqrt(tot_Eta_prec) * chol_R_invT * Z * X
//   // Phi = U * X = Vt * X, where U = Vt = sqrt(tot_Eta_prec) * chol_R_invT * Z
//   MatrixXd U = chol_R.transpose().triangularView<Lower>().solve(Z.toDense() * sqrt(tot_Eta_prec));
//   MatrixXd Phi = U * X;
//   MatrixXd V = U.transpose();
//   VectorXd alpha = chol_R.transpose().triangularView<Lower>().solve(y * sqrt(tot_Eta_prec));
//
//   VectorXd u = randn_theta.array() / prior_prec.cwiseSqrt().array();
//   u += prior_mean;
//   VectorXd v = Phi * u + randn_e;
//   VectorXd alpha_v = alpha-v;
//
// // Binomial inverse theorm: (A + UBV)^-1 = A^-1 - A^-1 * U (I + BVU)^-1*BVA^-1
//   MatrixXd B = X * prior_prec.cwiseInverse().asDiagonal() * X.transpose();
//   MatrixXd UB = U*B;
//   MatrixXd inner = B + UB.transpose()*UB;
//   VectorXd w = alpha_v - UB*inner.ldlt().solve(UB.transpose() * alpha_v);
//
//
//   // MatrixXd I = MatrixXd::Identity(B.rows(),B.cols());
//   // MatrixXd Binv = B.ldlt().solve(I);
//   // MatrixXd inner = Binv + V * U;
//   // MatrixXd inner = B + B*V*U*B;
//   // VectorXd w = alpha_v - U*B*inner.ldlt().solve(B*V * alpha_v);
//   // MatrixXd inner = B*V*U;
//   // inner.diagonal() += VectorXd::Ones(B.rows());
//   // VectorXd w = alpha_v - U*inner.ldlt().solve(B*V * alpha_v);
//
//   // MatrixXd BV = X * prior_prec.cwiseInverse().asDiagonal() * Phi.transpose();
//   // MatrixXd I = MatrixXd::Identity(BV.rows(),BV.cols());
//   // MatrixXd inner = BV*U;
//   // inner.diagonal() += VectorXd::Ones(inner.rows());
//   // VectorXd w = alpha_v - U*inner.ldlt().solve(BV * alpha_v);
//
//   // MatrixXd Sigma = Phi * prior_prec.cwiseInverse().asDiagonal() * Phi.transpose();
//   // Sigma.diagonal() += VectorXd::Ones(Sigma.rows());
//   // VectorXd w = Sigma.ldlt().solve(alpha_v);
//
//   VectorXd theta = u + prior_prec.cwiseInverse().asDiagonal() * (Phi.transpose() * w);
//
//   return(theta);
// }
//
// struct sample_MME_single_hierarchical_diagK_worker2 : public RcppParallel::Worker {
//   MatrixXd Y;
//   MSpMat Z;
//   MatrixXd X;
//   MatrixXd prior_mean, prior_prec, randn_theta, randn_e;
//   const std::vector<MSpMat> chol_R_list;
//   VectorXi h2s_index;
//   VectorXd tot_Eta_prec;
//   MatrixXd &coefs;
//
//   sample_MME_single_hierarchical_diagK_worker2(
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
//       coefs.col(j) = sample_MME_single_hierarchical_diagK2(Y.col(j), Z, X, prior_mean.col(j), prior_prec.col(j), chol_R, tot_Eta_prec[j], randn_theta.col(j),randn_e.col(j));
//     }
//   }
// };
//
// // Samples from B in model:
// // Y = ZXB + E
// // [[Rcpp::export()]]
// MatrixXd sample_MME_fixedEffects_hierarchical_c2(  // returns bxp matrix
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
//   sample_MME_single_hierarchical_diagK_worker2 sampler(Y,Z,X,prior_mean,prior_prec,chol_R_list,h2s_index,tot_Eta_prec,randn_theta,randn_e, coefs);
//   RcppParallel::parallelFor(0,p,sampler,grainSize);
//
//   return(coefs);
// }
//
