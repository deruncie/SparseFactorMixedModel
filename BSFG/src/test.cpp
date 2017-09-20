// #include <math.h>
// #include <iostream>
// #include "BSFG_types.h"
//
// // [[Rcpp::depends(RcppEigen)]]
// using namespace Eigen;
// using namespace RcppParallel;
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
//     MatrixXd RinvSqX = chol_Ri * (X * sqrt(tot_Eta_prec));
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
