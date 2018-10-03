// #include <math.h>
// #include <iostream>
// #include "BSFG_types.h"
//
// // [[Rcpp::depends(RcppEigen)]]
// using namespace Eigen;
// using namespace RcppParallel;
//
// // [[Rcpp::export()]]
// VectorXd rgamma2(int n,double shape,double scale){
//   VectorXd res(n);
//   for(int i = 0; i< n; i++){
//     res[i] = R::rgamma(shape,scale);
//   }
//   return(res);
// }
//

//
//
// // [[Rcpp::export()]]
// Rcpp::List sample_trunc_delta_omega_c_Eigen(
//     VectorXd delta,
//     VectorXd tauh,
//     double omega2,
//     double xi,
//     Map<VectorXd> scores,
//     Map<VectorXd> shapes,
//     double delta_1_rate,
//     double delta_2_rate,
//     Map<MatrixXd> randu_draws,
//     double trunc_point
// ) {
//   int times = randu_draws.rows();
//   int k = tauh.size();
//
//   double rate,delta_old, u,p;
//   for(int i = 0; i < times; i++){
//
//     rate = 1.0/xi + tauh.dot(scores);
//     u = randu_draws(i,0); // don't truncate delta(0)
//     omega2 = 1.0 / R::qgamma(u,shapes(0),1.0/rate,1,0);
//
//     rate = 1.0 + 1.0 / omega2;
//     u = randu_draws(i,1); // don't truncate delta(0)
//     xi = 1.0 / R::qgamma(u,shapes(1),1.0/rate,1,0);
//
//     VectorXd std_scores = scores / omega2;
//
//     delta_old = delta(0);
//     rate = delta_1_rate + (1/delta(0)) * tauh.dot(std_scores);
//     u = randu_draws(i,2); // don't truncate delta(0)
//     delta(0) = R::qgamma(u,shapes(2),1.0/rate,1,0);
//     // tauh = cumprod(delta);
//     tauh *= delta(0)/delta_old;   // replaces re-calculating cumprod
//
//     for(int h = 1; h < k; h++) {
//       delta_old = delta(h);
//       rate = delta_2_rate + (1/delta(h))*tauh.tail(k-h).dot(std_scores.tail(k-h));
//       p = R::pgamma(trunc_point,shapes(2+h),1.0/rate,1,0);  // left-tuncate delta(h) at trunc_point
//       if(p > 0.999) p = 0.999;  // prevent over-flow.
//       u = p + (1.0-p)*randu_draws(i,2+h);
//       delta(h) = R::qgamma(u,shapes(2+h),1.0/rate,1,0);
//       // tauh = cumprod(delta);
//       tauh.tail(k-h) *= delta(h)/delta_old; // replaces re-calculating cumprod
//       // Rcout << (tauh - cumprod(delta)).sum() << std::endl;
//     }
//   }
//   return(Rcpp::List::create(omega2,xi,delta));
// }
//
//
//
//
//
//
// // [[Rcpp::export()]]
// VectorXd sample_delta_c_Eigen_v2(
//     VectorXd delta,
//     Map<VectorXd> scores,
//     double delta_rate,
//     Map<MatrixXd> randg_draws  // all done with rate = 1;
// ) {
//   int times = randg_draws.rows();
//   int k = delta.size();
//
//   VectorXd cumprod_delta = cumprod(delta);
//   double rate,delta_old;
//   for(int i = 0; i < times; i++){
//     for(int h = 0; h < k; h++){
//       rate = delta_rate + (1/delta[h]) * cumprod_delta.tail(k-h).dot(scores.tail(k-h));
//       delta[h] = randg_draws(i,h) / rate;
//       cumprod_delta = cumprod(delta);
//     }
//   }
//   return(delta);
// }
//
//
// // [[Rcpp::export()]]
// VectorXd sample_trunc_delta_c_Eigen_v2(
//     VectorXd delta,
//     Map<VectorXd> scores,
//     Map<VectorXd> shapes,
//     double delta_rate,
//     Map<MatrixXd> randu_draws,
//     double trunc_point
// ) {
//   int times = randu_draws.rows();
//   int k = delta.size();
//   double p,u;
//   VectorXd cumprod_delta = cumprod(delta);
//
//   double rate,delta_old;
//   for(int i = 0; i < times; i++){
//     for(int h = 0; h < k; h++) {
//       rate = delta_rate + (1/delta[h]) * cumprod_delta.tail(k-h).dot(scores.tail(k-h));
//       p = R::pgamma(trunc_point,shapes(h),1.0/rate,1,0);  // left-tuncate delta(h) at trunc_point
//       if(p > 0.999) p = 0.999;  // prevent over-flow.
//       u = p + (1.0-p)*randu_draws(i,h);
//       delta(h) = R::qgamma(u,shapes(h),1.0/rate,1,0);
//       cumprod_delta = cumprod(delta);
//     }
//   }
//   return(delta);
// }


//
// // [[Rcpp::depends(RcppZiggurat)]]
//
// // [[Rcpp::export()]]
// VectorXd diag_crossprod(Map<MatrixXd> X){
//   return((X.transpose() * X).diagonal());
// }
//
// // [[Rcpp::export()]]
// VectorXd diag_crossprod2(Map<MatrixXd> X){
//   MatrixXd Y = X.transpose() * X;
//   return(Y.diagonal());
// }
//
// // [[Rcpp::export()]]
// VectorXd diag_crossprod3(Map<MatrixXd> X){
//   VectorXd res(X.cols());
//   for(int i = 0; i < X.cols(); i++){
//     res[i] = X.col(i).dot(X.col(i));
//   }
//   return(res);
// }
//
// // // [[Rcpp::export()]]
// // List LDLt(SEXP A_) {
// //   // returns P, L, d st PtLDLtP = A_
// //   if(Rf_isMatrix(A_)){
// //     Map<MatrixXd> A = as<Map<MatrixXd> >(A_);
// //     Eigen::LDLT<MatrixXd> ldlt_A;
// //     ldlt_A.compute(A);
// //     MatrixXd I = MatrixXd::Identity(ldlt_A.rows(), ldlt_A.rows());
// //     MatrixXd P = ldlt_A.transpositionsP() * I;
// //     VectorXd d = ldlt_A.vectorD();
// //     MatrixXd L = ldlt_A.matrixL();
// //     return(List::create(
// //         Named("P") = P.sparseView(),
// //         Named("L") = L.sparseView(),
// //         Named("d") = d
// //     ));
// //   } else{
// //     MSpMat A = as<MSpMat>(A_);
// //     Eigen::SimplicialLDLT<SpMat> ldlt_A;
// //     ldlt_A.compute(A);
// //     MatrixXd I = MatrixXd::Identity(ldlt_A.rows(), ldlt_A.rows());
// //     MatrixXd P = ldlt_A.permutationP() * I;
// //     return(List::create(
// //         Named("P") = P.sparseView(),
// //         Named("L") = ldlt_A.matrixL(),
// //         Named("d") = ldlt_A.vectorD()
// //     ));
// //   }
// // }
// //
// // SpMat make_chol_K_inv(Rcpp::List chol_Ki_mats, VectorXd h2s,double tol){
// //   int h = h2s.size();
// //   VectorXd sizes(h);
// //   int total_size = 0;
// //   for(int i = 0; i < h; i++){
// //     MSpMat chol_Ki = as<MSpMat>(chol_Ki_mats[i]);
// //     sizes(i) = chol_Ki.rows();
// //     total_size += sizes(i);
// //   }
// //   MatrixXd chol_K_inv_dense(total_size,total_size);
// //   chol_K_inv_dense.setZero();
// //   int curr_row = 0;
// //   int curr_col = 0;
// //   for(int i = 0; i < h; i++){
// //     if(h2s[i] == 0) {
// //       chol_K_inv_dense.block(curr_row,curr_col,sizes(i),sizes(i)).diagonal().setOnes();
// //       chol_K_inv_dense.block(curr_row,curr_col,sizes(i),sizes(i)).diagonal() /= 0;
// //     } else{
// //       MSpMat chol_Ki = as<MSpMat>(chol_Ki_mats[i]);
// //       chol_K_inv_dense.block(curr_row,curr_col,sizes(i),sizes(i)) = chol_Ki;
// //       chol_K_inv_dense.block(curr_row,curr_col,sizes(i),sizes(i)) /= sqrt(h2s[i]);
// //     }
// //     curr_row += sizes(i);
// //     curr_col += sizes(i);
// //   }
// //   SpMat chol_K_inv = chol_K_inv_dense.sparseView(0,tol);
// //   return chol_K_inv;
// // }
// //
// // struct make_chol_ZtZ_Kinv_worker : public RcppParallel::Worker {
// //   const Rcpp::List chol_Ki_mats;
// //   const Map<MatrixXd> h2s_matrix;
// //   const MSpMat ZtZ;
// //   const double tol;
// //   std::vector<SpMat>& chol_ZtZ_Kinv_list;
// //
// //   make_chol_ZtZ_Kinv_worker(
// //     const Rcpp::List chol_Ki_mats,
// //     const Map<MatrixXd> h2s_matrix,
// //     const MSpMat ZtZ,
// //     const double tol,
// //     std::vector<SpMat>& chol_ZtZ_Kinv_list
// //   ):
// //     chol_Ki_mats(chol_Ki_mats), h2s_matrix(h2s_matrix), ZtZ(ZtZ), tol(tol),
// //     chol_ZtZ_Kinv_list(chol_ZtZ_Kinv_list) {}
// //
// //   void operator()(std::size_t begin, std::size_t end) {
// //     for(std::size_t j = begin; j < end; j++){
// //       VectorXd h2s = h2s_matrix.col(j);
// //       SpMat chol_K_inv = make_chol_K_inv(chol_Ki_mats,h2s,tol);
// //       MatrixXd ZtZ_Kinv = 1.0/(1.0 - h2s.sum()) * ZtZ + chol_K_inv.transpose() * chol_K_inv;
// //       Eigen::LLT<MatrixXd> chol_ZtZ_Kinv(ZtZ_Kinv);
// //       MatrixXd chol_ZtZ_Kinv_R = chol_ZtZ_Kinv.matrixU();
// //       chol_ZtZ_Kinv_list[j] = chol_ZtZ_Kinv_R.sparseView(0,tol);
// //     }
// //   }
// // };
// //
// // // [[Rcpp::export()]]
// // Rcpp::List make_chol_ZtZ_Kinv_list(Rcpp::List chol_Ki_mats,
// //                                    Map<MatrixXd> h2s_matrix,
// //                                    MSpMat ZtZ,
// //                                    double drop0_tol,
// //                                    SEXP pb, Function setTxtProgressBar, Function getTxtProgressBar,
// //                                    int ncores) {
// //   int s = h2s_matrix.cols();
// //
// //   std::vector<SpMat> chol_ZtZ_Kinv_list;
// //   chol_ZtZ_Kinv_list.reserve(s);
// //   for(int i = 0; i < s; i++){
// //     chol_ZtZ_Kinv_list.push_back(SpMat(0,0));
// //   }
// //
// //   int n_groups = s/ncores;
// //   for(int i = 0; i <= n_groups; i++){
// //     make_chol_ZtZ_Kinv_worker make_chols(chol_Ki_mats,h2s_matrix,ZtZ,drop0_tol,chol_ZtZ_Kinv_list);
// //     int start = ncores*i;
// //     int end = std::min(static_cast<double>(s),ncores*(i+1.0));
// //     RcppParallel::parallelFor(start,end,make_chols,1);
// //     int pb_state = as<int>(getTxtProgressBar(pb));
// //     setTxtProgressBar(pb,pb_state+(end-start));
// //   }
// //
// //   Rcpp::List chol_ZtZ_Kinv_list_out(s);
// //   for(int i = 0; i < s; i++){
// //     chol_ZtZ_Kinv_list_out[i] = chol_ZtZ_Kinv_list[i];
// //   }
// //
// //   return(chol_ZtZ_Kinv_list_out);
// // }
// //
// //
// //
// // SpMat make_chol_R(Rcpp::List ZKZts, const VectorXd h2s, const double tol){  //std::vector<Map<MatrixXd> > ZKZts
// //   Map<MatrixXd> ZKZts_0 = as<Map<MatrixXd> >(ZKZts[0]);
// //   int n = ZKZts_0.rows();
// //   int h = h2s.size();
// //   MatrixXd R(n,n);
// //   R.setZero();
// //   for(int i = 0; i < h; i++){
// //     Map<MatrixXd> ZKZts_i = as<Map<MatrixXd> >(ZKZts[i]);
// //     R += h2s[i] * ZKZts_i;
// //   }
// //   R.diagonal().array() += (1.0-h2s.sum());
// //
// //   Eigen::LLT<MatrixXd> chol_R(R);
// //   MatrixXd chol_R_U = chol_R.matrixU();
// //   return chol_R_U.sparseView(0,tol);
// // }
// //
// //
// // struct make_chol_R_worker : public RcppParallel::Worker {
// //   // const std::vector<Map<MatrixXd> > ZKZts;
// //   const Rcpp::List ZKZts;
// //   const Map<MatrixXd> h2s_matrix;
// //   const double tol;
// //   std::vector<SpMat>& chol_R_list;
// //
// //   make_chol_R_worker(
// //     const Rcpp::List ZKZts, //const std::vector<Map<MatrixXd> > ZKZts,
// //     const Map<MatrixXd> h2s_matrix,
// //     const double tol,
// //     std::vector<SpMat>& chol_R_list
// //   ):
// //     ZKZts(ZKZts), h2s_matrix(h2s_matrix), tol(tol),
// //     chol_R_list(chol_R_list) {}
// //
// //   void operator()(std::size_t begin, std::size_t end) {
// //     for(std::size_t j = begin; j < end; j++){
// //       chol_R_list[j] = make_chol_R(ZKZts, h2s_matrix.col(j), tol);
// //     }
// //   }
// // };
// //
// // // [[Rcpp::export()]]
// // Rcpp::List make_chol_R_list(Rcpp::List ZKZts,
// //                             Map<MatrixXd> h2s_matrix,
// //                             double drop0_tol,
// //                             SEXP pb, Function setTxtProgressBar, Function getTxtProgressBar,
// //                             int ncores) {
// //   int s = h2s_matrix.cols();
// //
// //   std::vector<SpMat> chol_R_list;
// //   chol_R_list.reserve(s);
// //   for(int i = 0; i < s; i++){
// //     chol_R_list.push_back(SpMat(0,0));
// //   }
// //
// //   int n_groups = s/ncores;
// //   for(int i = 0; i <= n_groups; i++){
// //     make_chol_R_worker make_chols(ZKZts,h2s_matrix,drop0_tol,chol_R_list);
// //     int start = ncores*i;
// //     int end = std::min(static_cast<double>(s),ncores*(i+1.0));
// //     RcppParallel::parallelFor(start,end,make_chols,1);
// //     int pb_state = as<int>(getTxtProgressBar(pb));
// //     setTxtProgressBar(pb,pb_state+(end-start));
// //   }
// //
// //   Rcpp::List chol_R_list_out(s);
// //   for(int i = 0; i < s; i++){
// //     chol_R_list_out[i] = chol_R_list[i];
// //   }
// //
// //   return(chol_R_list_out);
// // }
//
// // // [[Rcpp::export()]]
// // MatrixXd matrix_multiply_toDense(SEXP X_, SEXP Y_){
// //   if(Rf_isNull(X_)) return(as<Map<MatrixXd> >(Y_));
// //   if(Rf_isMatrix(X_)) {
// //     Map<MatrixXd> X = as<Map<MatrixXd> >(X_);
// //     if(Rf_isMatrix(Y_)) {
// //       Map<MatrixXd> Y = as<Map<MatrixXd> >(Y_);
// //       if(X.cols() != Y.rows()) stop("Wrong dimensions of matrices");
// //       return(X*Y);
// //     } else{
// //       MSpMat Y = as<MSpMat>(Y_);
// //       if(X.cols() != Y.rows()) stop("Wrong dimensions of matrices");
// //       return(X*Y);
// //     }
// //   }
// //   else {
// //     MSpMat X = as<MSpMat>(X_);
// //     if(Rf_isMatrix(Y_)) {
// //       Map<MatrixXd> Y = as<Map<MatrixXd> >(Y_);
// //       if(X.cols() != Y.rows()) stop("Wrong dimensions of matrices");
// //       return(X*Y);
// //     } else{
// //       MSpMat Y = as<MSpMat>(Y_);
// //       if(X.cols() != Y.rows()) stop("Wrong dimensions of matrices");
// //       return(X*Y);
// //     }
// //   }
// // }
// //
// //
// // // [[Rcpp::depends(RcppZiggurat)]]
// // #include <Ziggurat.h>
// // #include <ZigguratR.h>
// // using namespace Rcpp;
// static Ziggurat::Ziggurat::Ziggurat zigg;
//
// // [[Rcpp::export()]]
// MatrixXd rstdnorm_mat2(int n,int p) {  // returns nxp matrix
//   MatrixXd X_mat = as<VectorXd>(rnorm(n*p));
//   // Map<MatrixXd> X_mat(X_vec.data(),n,p);
//   X_mat.resize(n,p);
//   return(X_mat);
// }
// //
// // // [[Rcpp::export()]]
// // MatrixXd rstdnorm_mat3(int n,int p) {  // returns nxp matrix
// //   VectorXd X_vec(n*p);
// //   for(int i = 0; i < n*p; i++){
// //     X_vec[i] = zigg.norm();
// //   }
// //   MatrixXd X_mat = Map<MatrixXd>(X_vec.data(),n,p);
// //   return(X_mat);
// // }
// //
// // // static Ziggurat::R::ZigguratR ziggr;
// // // [[Rcpp::export()]]
// // MatrixXd rstdnorm_mat4(int n,int p) {  // returns nxp matrix
// //   VectorXd X_vec(n*p);
// //   for(int i = 0; i < n*p; i++){
// //     X_vec[i] = ziggr.norm();
// //   }
// //   MatrixXd X_mat = Map<MatrixXd>(X_vec.data(),n,p);
// //   return(X_mat);
// // }
// //
// //
// //
// //
// // // [[Rcpp::export()]]
// // unsigned long int myzigg_getSeed(){
// //   return(zigg.getSeed());
// // }
// // // [[Rcpp::export()]]
// // void myzigg_setSeed(unsigned long int s){
// //   zigg.setSeed(s);
// // }
// //
// // // [[Rcpp::export()]]
// // void test_rstdnorm2(int n,int p) {
// //   MatrixXd x = rstdnorm_mat2(n,p);
// // }
// // // [[Rcpp::export()]]
// // void test_rstdnorm3(int n,int p) {
// //   MatrixXd x = rstdnorm_mat3(n,p);
// // }
// //
// // // // [[Rcpp::export]]
// // // uint32_t test(uint32_t jsr) {
// // //   uint32_t jz = jsr;
// // //   jsr^=(jsr<<13);
// // //   jsr^=(jsr>>17);
// // //   jsr^=(jsr<<5);
// // //   return(jz+jsr);
// // // }
// //
// //
// // // [[Rcpp::export()]]
// // unsigned long int ce(unsigned long int i,unsigned long int j){
// //   i ^= j;
// //   return(i);
// // }
// //
// // // -------------------------------------------- //
// // // ---------- sample_MME_fixedEffects --------- //
// // // -------------------------------------------- //
// //
// // VectorXd sample_MME_single_diagK2(  // returns b x 1 vector
// //     const Ref<const VectorXd>& y,           // nx1
// //     const MatrixXd& X,                      // nxb
// //     const Ref<const VectorXd>& prior_mean,  // bx1
// //     const Ref<const VectorXd>& prior_prec,  // bx1
// //     const SEXP& chol_R_,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
// //     double tot_Eta_prec,                    // double
// //     const Ref<const VectorXd>& randn_theta, // bx1
// //     const Ref<const VectorXd>& randn_e      // 0x1 or nx1. 0x1 if b<n
// // ){
// //   if(randn_e.size() == 0){
// //     MatrixXd RinvSqX;
// //     VectorXd XtRinvy;
// //     if(Rf_isMatrix(chol_R_)) {
// //       Map<MatrixXd> chol_R = as<Map<MatrixXd> >(chol_R_);
// //       RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
// //       XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y * sqrt(tot_Eta_prec));
// //     } else{
// //       MSpMat chol_R = as<MSpMat>(chol_R_);
// //       RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
// //       XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y * sqrt(tot_Eta_prec));
// //     }
// //
// //     VectorXd XtRinvy_std_mu = XtRinvy + prior_prec.asDiagonal()*prior_mean;
// //     MatrixXd C = RinvSqX.transpose() * RinvSqX;
// //     C.diagonal() += prior_prec;
// //     LLT<MatrixXd> C_llt;
// //     C_llt.compute(C);
// //     MatrixXd chol_C = C_llt.matrixU();
// //
// //     VectorXd b = chol_C.transpose().triangularView<Lower>().solve(XtRinvy_std_mu);
// //     b += randn_theta;
// //     b = chol_C.triangularView<Upper>().solve(b);
// //     return(b);
// //   } else {
// //     // Using algorithm from Bhattacharya et al 2016 Biometrika. https://academic.oup.com/biomet/article/103/4/985/2447851
// //     MatrixXd Phi;
// //     VectorXd alpha;
// //     if(Rf_isMatrix(chol_R_)) {
// //       Map<MatrixXd> chol_R = as<Map<MatrixXd> >(chol_R_);
// //       Phi = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
// //       alpha = chol_R.transpose().triangularView<Lower>().solve(y * sqrt(tot_Eta_prec));
// //     } else{
// //       MSpMat chol_R = as<MSpMat>(chol_R_);
// //       Phi = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
// //       alpha = chol_R.transpose().triangularView<Lower>().solve(y * sqrt(tot_Eta_prec));
// //     }
// //
// //     VectorXd u = randn_theta.array() / prior_prec.cwiseSqrt().array();
// //     u += prior_mean;
// //     VectorXd v = Phi * u + randn_e;
// //     VectorXd alpha_v = alpha-v;
// //
// //     MatrixXd D_PhiT = prior_prec.cwiseInverse().asDiagonal() * Phi.transpose();
// //     MatrixXd cov = Phi * D_PhiT;
// //     cov.diagonal().array() += 1.0;
// //
// //     VectorXd w = cov.ldlt().solve(alpha_v);
// //
// //     VectorXd theta = u + D_PhiT * w;
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
// // struct sample_MME_single_diagK_worker2 : public RcppParallel::Worker {
// //   const Map<MatrixXd>& Y;
// //   const Map<MatrixXd>& X;
// //   const Map<MatrixXd>& prior_mean, prior_prec;
// //   const MatrixXd& randn_theta, randn_e;
// //   // const std::vector<MSpMat> chol_R_list;
// //   const Rcpp::List &chol_R_list;
// //   const VectorXi h2s_index;
// //   const VectorXd tot_Eta_prec;
// //   MatrixXd &coefs;
// //
// //   sample_MME_single_diagK_worker2(
// //     const Map<MatrixXd> &Y,           // nxp
// //     const Map<MatrixXd>& X,           // nxb
// //     const Map<MatrixXd>& prior_mean,  // bxp
// //     const Map<MatrixXd>& prior_prec,  // bxp    // const std::vector<MSpMat> &chol_R_list, // each element a MSpMat nxn upper-triangular cholesky decomposition of R
// //     const Rcpp::List &chol_R_list,
// //     const VectorXi h2s_index,   // px1, 1-based index
// //     const VectorXd tot_Eta_prec,// px1
// //     const MatrixXd& randn_theta, // bxp
// //     const MatrixXd& randn_e,     // 0xp or nxp. The former if b<n, the latter if b >= n
// //     MatrixXd &coefs       // bxp
// //   ):
// //     Y(Y), X(X), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
// //     chol_R_list(chol_R_list), h2s_index(h2s_index), tot_Eta_prec(tot_Eta_prec), coefs(coefs) {}
// //
// //   void operator()(std::size_t begin, std::size_t end) {
// //     for(std::size_t j = begin; j < end; j++){
// //       int h2_index = h2s_index[j] - 1;
// //       const Rcpp::List Sigma_Choleskys_i = chol_R_list[h2_index];
// //       const SEXP chol_R = Sigma_Choleskys_i["chol_Sigma"];
// //       coefs.col(j) = sample_MME_single_diagK2(Y.col(j), X, prior_mean.col(j), prior_prec.col(j), chol_R, tot_Eta_prec[j], randn_theta.col(j),randn_e.col(j));
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
// // MatrixXd sample_MME_fixedEffects_c2(  // returns bxp matrix
// //     Map<MatrixXd> Y,              // nxp
// //     Map<MatrixXd> X,              // nxb
// //     Rcpp::List chol_R_list,   // list of Cholesky decompositions of residuals. nxn upper triangular, dgCMatrix
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
// //   // Map<MatrixXd> randn_theta(as<VectorXd>(rnorm(b*p)).data(),b,p);
// //   // // MatrixXd randn_theta = rstdnorm_mat2(b,p);
// //   // Map<MatrixXd> randn_e(as<VectorXd>(rnorm(0)).data(),0,p);
// //   // if(b < n) {
// //   //   // Rcout << "b<n" << std::endl;
// //   //   randn_e = rstdnorm_mat2(0,p);
// //   // } else{
// //   //   // Rcout << "b>=n" << std::endl;
// //   //   new (&randn_e) Map<MatrixXd>(as<VectorXd>(rnorm(n*p)).data(),n,p);
// //   // }
// //   MatrixXd randn_theta = rstdnorm_mat2(b,p);
// //   MatrixXd randn_e;
// //   if(b < n) {
// //     randn_e = rstdnorm_mat2(0,p);
// //   } else{
// //     randn_e = rstdnorm_mat2(n,p);
// //   }
// //   if(randn_theta.rows() != b || randn_theta.cols() != p) stop("wrong dimension of randn_theta");
// //   if(randn_e.cols() != p) stop("wrong dimension of randn_e");
// //
// //   // std::vector<MSpMat> chol_R_list;
// //   // for(int i = 0; i < h2s_index.maxCoeff(); i++){
// //   //   Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);
// //   //   chol_R_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
// //   // }
// //
// //   MatrixXd coefs(b,p);
// //
// //   sample_MME_single_diagK_worker2 sampler(Y,X,prior_mean,prior_prec,chol_R_list,h2s_index,tot_Eta_prec,randn_theta,randn_e, coefs);
// //   RcppParallel::parallelFor(0,p,sampler,grainSize);
// //   return(coefs);
// // }
// //
// // struct sample_MME_single_diagK_worker3 : public RcppParallel::Worker {
// //   MatrixXd Y;
// //   MatrixXd X;
// //   MatrixXd prior_mean, prior_prec;
// //   MatrixXd randn_theta, randn_e;
// //   // const std::vector<MSpMat> chol_R_list;
// //   Rcpp::List &chol_R_list;
// //   VectorXi h2s_index;
// //   VectorXd tot_Eta_prec;
// //   MatrixXd &coefs;
// //
// //   sample_MME_single_diagK_worker3(
// //     MatrixXd Y,           // nxp
// //     MatrixXd X,           // nxb
// //     MatrixXd prior_mean,  // bxp
// //     MatrixXd prior_prec,  // bxp    // const std::vector<MSpMat> &chol_R_list, // each element a MSpMat nxn upper-triangular cholesky decomposition of R
// //     Rcpp::List &chol_R_list,
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
// //       Rcpp::List Sigma_Choleskys_i = chol_R_list[h2_index];
// //       SEXP chol_R = Sigma_Choleskys_i["chol_Sigma"];
// //       coefs.col(j) = sample_MME_single_diagK2(Y.col(j), X, prior_mean.col(j), prior_prec.col(j), chol_R, tot_Eta_prec[j], randn_theta.col(j),randn_e.col(j));
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
// // MatrixXd sample_MME_fixedEffects_c3(  // returns bxp matrix
// //     Map<MatrixXd> Y,              // nxp
// //     Map<MatrixXd> X,              // nxb
// //     Rcpp::List chol_R_list,   // list of Cholesky decompositions of residuals. nxn upper triangular, dgCMatrix
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
// //   // Map<MatrixXd> randn_theta(as<VectorXd>(rnorm(b*p)).data(),b,p);
// //   // // MatrixXd randn_theta = rstdnorm_mat2(b,p);
// //   // Map<MatrixXd> randn_e(as<VectorXd>(rnorm(0)).data(),0,p);
// //   // if(b < n) {
// //   //   // Rcout << "b<n" << std::endl;
// //   //   randn_e = rstdnorm_mat2(0,p);
// //   // } else{
// //   //   // Rcout << "b>=n" << std::endl;
// //   //   new (&randn_e) Map<MatrixXd>(as<VectorXd>(rnorm(n*p)).data(),n,p);
// //   // }
// //   MatrixXd randn_theta = rstdnorm_mat4(b,p);
// //   MatrixXd randn_e;
// //   if(b < n) {
// //     randn_e = rstdnorm_mat4(0,p);
// //   } else{
// //     randn_e = rstdnorm_mat4(n,p);
// //   }
// //   if(randn_theta.rows() != b || randn_theta.cols() != p) stop("wrong dimension of randn_theta");
// //   if(randn_e.cols() != p) stop("wrong dimension of randn_e");
// //
// //   // std::vector<MSpMat> chol_R_list;
// //   // for(int i = 0; i < h2s_index.maxCoeff(); i++){
// //   //   Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);
// //   //   chol_R_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
// //   // }
// //
// //   MatrixXd coefs(b,p);
// //
// //   sample_MME_single_diagK_worker3 sampler(Y,X,prior_mean,prior_prec,chol_R_list,h2s_index,tot_Eta_prec,randn_theta,randn_e, coefs);
// //   RcppParallel::parallelFor(0,p,sampler,grainSize);
// //   return(coefs);
// // }
// //
// //
// //
// //
// // //// New strategy
// //
// // VectorXd sample_MME_single_diagK4(  // returns b x 1 vector
// //     const Ref<const VectorXd>& y,           // nx1
// //     const MatrixXd& RinvSqX,                // nxb
// //     const MatrixXd& C_,                     // bxb
// //     const Ref<const VectorXd>& prior_mean,  // bx1
// //     const Ref<const VectorXd>& prior_prec,  // bx1
// //     const SEXP& chol_R_,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
// //     double tot_Eta_prec,                    // double
// //     const Ref<const VectorXd>& randn_theta, // bx1
// //     const Ref<const VectorXd>& randn_e      // 0x1 or nx1. 0x1 if b<n
// // ){
// //   if(randn_e.size() == 0){
// //     // MatrixXd RinvSqX;
// //     VectorXd XtRinvy;
// //     // RinvSqX
// //     if(Rf_isMatrix(chol_R_)) {
// //       Map<MatrixXd> chol_R = as<Map<MatrixXd> >(chol_R_);
// //       // RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
// //       XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y * tot_Eta_prec);
// //     } else{
// //       MSpMat chol_R = as<MSpMat>(chol_R_);
// //       // RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
// //       XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y * tot_Eta_prec);
// //     }
// //
// //     VectorXd XtRinvy_std_mu = XtRinvy + prior_prec.asDiagonal()*prior_mean;
// //     // MatrixXd C = RinvSqX.transpose() * RinvSqX * tot_Eta_prec;
// //     MatrixXd C = C_ * tot_Eta_prec;
// //     C.diagonal() += prior_prec;
// //     LLT<MatrixXd> C_llt;
// //     C_llt.compute(C);
// //     MatrixXd chol_C = C_llt.matrixU();
// //
// //     VectorXd b = chol_C.transpose().triangularView<Lower>().solve(XtRinvy_std_mu);
// //     b += randn_theta;
// //     b = chol_C.triangularView<Upper>().solve(b);
// //     return(b);
// //   } else {
// //     // Using algorithm from Bhattacharya et al 2016 Biometrika. https://academic.oup.com/biomet/article/103/4/985/2447851
// //     MatrixXd Phi = RinvSqX * sqrt(tot_Eta_prec);
// //     VectorXd alpha;
// //     if(Rf_isMatrix(chol_R_)) {
// //       Map<MatrixXd> chol_R = as<Map<MatrixXd> >(chol_R_);
// //       // Phi = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
// //       alpha = chol_R.transpose().triangularView<Lower>().solve(y * sqrt(tot_Eta_prec));
// //     } else{
// //       MSpMat chol_R = as<MSpMat>(chol_R_);
// //       // Phi = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
// //       alpha = chol_R.transpose().triangularView<Lower>().solve(y * sqrt(tot_Eta_prec));
// //     }
// //
// //     VectorXd u = randn_theta.array() / prior_prec.cwiseSqrt().array();
// //     u += prior_mean;
// //     VectorXd v = Phi * u + randn_e;
// //     VectorXd alpha_v = alpha-v;
// //
// //     MatrixXd D_PhiT = prior_prec.cwiseInverse().asDiagonal() * Phi.transpose();
// //     MatrixXd cov = Phi * D_PhiT;
// //     cov.diagonal().array() += 1.0;
// //
// //     VectorXd w = cov.ldlt().solve(alpha_v);
// //
// //     VectorXd theta = u + D_PhiT * w;
// //
// //     return(theta);
// //   }
// // }
// //
// //
// // struct sample_MME_single_diagK_worker4 : public RcppParallel::Worker {
// //   MatrixXd Y;
// //   MatrixXd X;
// //   MatrixXd prior_mean, prior_prec;
// //   MatrixXd randn_theta, randn_e;
// //   // const std::vector<MSpMat> chol_R_list;
// //   Rcpp::List &chol_R_list;
// //   VectorXi h2s_index;
// //   VectorXd tot_Eta_prec;
// //   MatrixXd &coefs;
// //
// //   sample_MME_single_diagK_worker4(
// //     MatrixXd Y,           // nxp
// //     MatrixXd X,           // nxb
// //     MatrixXd prior_mean,  // bxp
// //     MatrixXd prior_prec,  // bxp    // const std::vector<MSpMat> &chol_R_list, // each element a MSpMat nxn upper-triangular cholesky decomposition of R
// //     Rcpp::List &chol_R_list,
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
// //     for(std::size_t i = begin; i < end; i++){
// //       int h2_index = i;
// //       int num_columns = 0;
// //       MatrixXd RinvSqX;
// //       MatrixXd C;
// //       SEXP chol_R_;
// //       for(int j = 0; j < h2s_index.size(); j++){
// //         if(h2s_index[j] == h2_index) {
// //           if(num_columns == 0) {
// //             num_columns = 1;
// //             Rcpp::List Sigma_Choleskys_i = chol_R_list[h2_index - 1];
// //             chol_R_ = Sigma_Choleskys_i["chol_Sigma"];
// //             if(Rf_isMatrix(chol_R_)) {
// //               Map<MatrixXd> chol_R = as<Map<MatrixXd> >(chol_R_);
// //               RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
// //             } else{
// //               MSpMat chol_R = as<MSpMat>(chol_R_);
// //               RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
// //             }
// //             C = RinvSqX.transpose() * RinvSqX;
// //           }
// //           coefs.col(j) = sample_MME_single_diagK4(Y.col(j), RinvSqX,C, prior_mean.col(j), prior_prec.col(j), chol_R_, tot_Eta_prec[j], randn_theta.col(j),randn_e.col(j));
// //         }
// //       }
// //     }
// //   }
// // };
// //
// //
// //
// // // basic version - all regression have same dimensions
// // // [[Rcpp::export()]]
// // MatrixXd sample_MME_fixedEffects_c4(  // returns bxp matrix
// //     Map<MatrixXd> Y,              // nxp
// //     Map<MatrixXd> X,              // nxb
// //     Rcpp::List chol_R_list,   // list of Cholesky decompositions of residuals. nxn upper triangular, dgCMatrix
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
// //   MatrixXd randn_theta = rstdnorm_mat4(b,p);
// //   MatrixXd randn_e;
// //   if(b < n) {
// //     randn_e = rstdnorm_mat4(0,p);
// //   } else{
// //     randn_e = rstdnorm_mat4(n,p);
// //   }
// //   if(randn_theta.rows() != b || randn_theta.cols() != p) stop("wrong dimension of randn_theta");
// //   if(randn_e.cols() != p) stop("wrong dimension of randn_e");
// //
// //   MatrixXd coefs(b,p);
// //
// //   sample_MME_single_diagK_worker4 sampler(Y,X,prior_mean,prior_prec,chol_R_list,h2s_index,tot_Eta_prec,randn_theta,randn_e, coefs);
// //   RcppParallel::parallelFor(h2s_index.minCoeff(),h2s_index.maxCoeff()+1,sampler,grainSize);
// //   return(coefs);
// // }
// //
// // // new2
// //
// // struct sample_MME_single_diagK_worker5 : public RcppParallel::Worker {
// //   const MatrixXd& Y;
// //   const MatrixXd& X;
// //   const MatrixXd& RinvSqX;
// //   const MatrixXd& C;
// //   SEXP chol_R_;
// //   const IntegerVector& cols;
// //   const MatrixXd& prior_mean, prior_prec;
// //   const MatrixXd& randn_theta, randn_e;
// //   // const std::vector<MSpMat> chol_R_list;
// //   // Rcpp::List &chol_R_list;
// //   // VectorXi h2s_index;
// //   const VectorXd& tot_Eta_prec;
// //   MatrixXd &coefs;
// //
// //   sample_MME_single_diagK_worker5(
// //     const MatrixXd& Y,           // nxp
// //     const MatrixXd& X,           // nxb
// //     const MatrixXd& RinvSqX,
// //     const MatrixXd& C,
// //     SEXP chol_R_,
// //     const IntegerVector& cols,
// //     const MatrixXd& prior_mean,  // bxp
// //     const MatrixXd& prior_prec,  // bxp    // const std::vector<MSpMat> &chol_R_list, // each element a MSpMat nxn upper-triangular cholesky decomposition of R
// //     const VectorXd& tot_Eta_prec,// px1
// //     const MatrixXd& randn_theta, // bxp
// //     const MatrixXd& randn_e,     // 0xp or nxp. The former if b<n, the latter if b >= n
// //     MatrixXd &coefs       // bxp
// //   ):
// //     Y(Y), X(X), RinvSqX(RinvSqX), C(C), chol_R_(chol_R_), cols(cols),prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
// //     tot_Eta_prec(tot_Eta_prec), coefs(coefs) {}
// //
// //   void operator()(std::size_t begin, std::size_t end) {
// //     for(std::size_t i = begin; i < end; i++){
// //       int j = cols[i];
// //       coefs.col(j) = sample_MME_single_diagK4(Y.col(j), RinvSqX,C, prior_mean.col(j), prior_prec.col(j), chol_R_, tot_Eta_prec[j], randn_theta.col(j),randn_e.col(j));
// //     }
// //   }
// // };
// //
// //
// // // [[Rcpp::export]]
// // Rcpp::IntegerVector which2(Rcpp::LogicalVector x) {
// //   Rcpp::IntegerVector v = Rcpp::seq(0, x.size()-1);
// //   return v[x];
// // }
// //
// // // basic version - all regression have same dimensions
// // // [[Rcpp::export()]]
// // MatrixXd sample_MME_fixedEffects_c5(  // returns bxp matrix
// //     Map<MatrixXd> Y,              // nxp
// //     Map<MatrixXd> X,              // nxb
// //     Rcpp::List chol_R_list,   // list of Cholesky decompositions of residuals. nxn upper triangular, dgCMatrix
// //     IntegerVector h2s_index,           // px1 index of Cholesky matrix for each column
// //     Map<VectorXd> tot_Eta_prec,   // px1
// //     Map<MatrixXd> prior_mean,     // bxp
// //     Map<MatrixXd> prior_prec,     // bxp
// //     int grainSize) {
// //
// //   int b = X.cols();
// //   int p = Y.cols();
// //   int n = Y.rows();
// //
// //   MatrixXd randn_theta = rstdnorm_mat4(b,p);
// //   MatrixXd randn_e;
// //   if(b < n) {
// //     randn_e = rstdnorm_mat4(0,p);
// //   } else{
// //     randn_e = rstdnorm_mat4(n,p);
// //   }
// //   if(randn_theta.rows() != b || randn_theta.cols() != p) stop("wrong dimension of randn_theta");
// //   if(randn_e.cols() != p) stop("wrong dimension of randn_e");
// //
// //   MatrixXd coefs(b,p);
// //
// //   for(int i = min(h2s_index); i <= max(h2s_index); i++) {
// //     int h2_index = i;
// //     IntegerVector cols = which2(h2s_index == h2_index);
// //     if(cols.size() > 0){
// //       MatrixXd RinvSqX;
// //       MatrixXd C;
// //       Rcpp::List Sigma_Choleskys_i = chol_R_list[h2_index - 1];
// //       SEXP chol_R_ = Sigma_Choleskys_i["chol_Sigma"];
// //       if(Rf_isMatrix(chol_R_)) {
// //         Map<MatrixXd> chol_R = as<Map<MatrixXd> >(chol_R_);
// //         RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
// //       } else{
// //         MSpMat chol_R = as<MSpMat>(chol_R_);
// //         RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
// //       }
// //       C = RinvSqX.transpose() * RinvSqX;
// //       sample_MME_single_diagK_worker5 sampler(Y,X,RinvSqX,C,chol_R_,cols,prior_mean,prior_prec,tot_Eta_prec,randn_theta,randn_e, coefs);
// //       RcppParallel::parallelFor(0,cols.size(),sampler,grainSize);
// //     }
// //   }
// //   return(coefs);
// // }
// //
// //
// // // [[Rcpp::export()]]
// // VectorXd sample_MME_block(  // returns b x 1 vector
// //     VectorXd y,           // nx1
// //     MatrixXd X,           // nxb
// //     VectorXd prior_mean,  // bx1
// //     VectorXd prior_prec,  // bx1
// //     MSpMat chol_R,        // nxn upper triangular Cholesky decomposition of R. Format: dgCMatrix
// //     VectorXd randn_theta, // bx1
// //     VectorXd randn_e,      // 0x1 or nx1. 0x1 if b<n
// //     double rgamma_1
// // ){
// //   if(randn_e.size() == 0){
// //     MatrixXd RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X);
// //     VectorXd RinvSqy = chol_R.transpose().triangularView<Lower>().solve(y);
// //     VectorXd XtRinvy = RinvSqX.transpose() * RinvSqy;
// //
// //     MatrixXd C = RinvSqX.transpose() * RinvSqX;
// //     C.diagonal() += prior_prec;
// //     LLT<MatrixXd> C_llt;
// //     C_llt.compute(C);
// //     MatrixXd chol_C = C_llt.matrixU();
// //
// //     VectorXd prod1 = chol_C.transpose().triangularView<Lower>().solve(XtRinvy);
// //     double score = (RinvSqy.dot(RinvSqy) - prod1.dot(prod1))/2;
// //     double tot_Eta_prec = rgamma_1/score;
// //
// //     VectorXd XtRinvy_std_mu = XtRinvy*tot_Eta_prec + prior_prec.asDiagonal()*prior_mean;
// //     VectorXd b = chol_C.transpose().triangularView<Lower>().solve(XtRinvy_std_mu) / sqrt(tot_Eta_prec);
// //     b += randn_theta;
// //     b = chol_C.triangularView<Upper>().solve(b) / sqrt(tot_Eta_prec);
// //
// //     VectorXd result(b.size() + 1);
// //     result << tot_Eta_prec,b;
// //     return(result);
// //   } else {
// //     // Using algorithm from Bhattacharya et al 2016 Biometrika. https://academic.oup.com/biomet/article/103/4/985/2447851
// //     MatrixXd Phi = chol_R.transpose().triangularView<Lower>().solve(X);
// //     VectorXd alpha = chol_R.transpose().triangularView<Lower>().solve(y);
// //
// //     MatrixXd D_PhiT = prior_prec.cwiseInverse().asDiagonal() * Phi.transpose();
// //     MatrixXd cov = Phi * D_PhiT;
// //     cov.diagonal().array() += 1.0;
// //
// //     LDLT<MatrixXd> cov_ldlt;
// //     cov_ldlt.compute(cov);
// //
// //     VectorXd prod1 = cov_ldlt.solve(alpha);
// //     double score = alpha.dot(prod1)/2;
// //     double tot_Eta_prec = rgamma_1/score;
// //
// //     Phi *= sqrt(tot_Eta_prec);
// //     alpha *= sqrt(tot_Eta_prec);
// //     D_PhiT /= sqrt(tot_Eta_prec);
// //
// //     VectorXd u = randn_theta.array() / (prior_prec * tot_Eta_prec).cwiseSqrt().array();
// //     u += prior_mean;
// //     VectorXd v = Phi * u + randn_e;
// //     VectorXd alpha_v = alpha-v;
// //
// //     VectorXd w = cov_ldlt.solve(alpha_v);
// //
// //     VectorXd theta = u + D_PhiT * w;
// //
// //     VectorXd result(theta.size() + 1);
// //     result << tot_Eta_prec,theta;
// //     return(result);
// //   }
// // }
// //
// //
// //
// // MatrixXd get_RinvSqX(SEXP chol_R_, MatrixXd X){
// //   MatrixXd RinvSqX;
// //   if(Rf_isMatrix(chol_R_)) {
// //     Map<MatrixXd> chol_R = as<Map<MatrixXd> >(chol_R_);
// //     RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
// //   } else{
// //     MSpMat chol_R = as<MSpMat>(chol_R_);
// //     RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
// //   }
// //   return(RinvSqX);
// // }
// //
// //
// // VectorXd regression_sampler_v1(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b < n
// //     const Ref<const VectorXd>& y,           // nx1
// //     const MatrixXd& W,           // nxa
// //     const MatrixXd& RinvSqX,                // nxb
// //     const MatrixXd& C,                     // bxb
// //     const VectorXd& prior_prec_alpha, // ax 1
// //     const Ref<const VectorXd>& prior_mean_beta,  // bx1
// //     const Ref<const VectorXd>& prior_prec_beta,  // bx1
// //     const SEXP chol_R_,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
// //     double Y_prec,                    // double
// //     const VectorXd& randn_alpha,
// //     const Ref<const VectorXd>& randn_beta,
// //     const double rgamma_1,
// //     const double Y_prec_b0
// // ){
// //   int n = y.size();
// //   int a = W.cols();
// //   int b = RinvSqX.cols();
// //
// //   // Check inputs
// //   if(W.rows() != n) stop("Wrong dimension of W");
// //   if(RinvSqX.rows() != n) stop("Wrong dimension of X");
// //   if(C.rows() != b || C.cols() != b) stop("Wrong dimension of C");
// //   if(prior_prec_alpha.size() != a) stop("Wrong length of prior_prec_alpha");
// //   if(prior_mean_beta.size() != b) stop("Wrong length of prior_mean_beta");
// //   if(prior_prec_beta.size() != b) stop("Wrong length of prior_prec_beta");
// //   if(randn_alpha.size() != a) stop("Wrong length of randn_alpha");
// //   if(randn_beta.size() != b) stop("Wrong length of randn_beta");
// //
// //   // Calculate cholesky of A_beta
// //   // C = Xt(RtR)^-1X
// //   MatrixXd C_beta = C;
// //   C_beta.diagonal() += prior_prec_beta;
// //   LLT<MatrixXd> A_beta_llt;
// //   A_beta_llt.compute(C_beta);
// //   MatrixXd chol_A_beta = A_beta_llt.matrixU();
// //   // Y_prec * chol_A_beta^\T * chol_A_beta = A_beta
// //
// //   // Step 1
// //   VectorXd alpha(a);
// //   VectorXd y_tilde = y;
// //   if(a > 0) {
// //     // Sample alpha
// //     // Calculate A_alpha = Y_prec*W^T*Sigma_beta^{-1}*W + D_alpha^{-1}
// //     // We don't need to actually calculate Sigma_beta^{-1} directly.
// //     MatrixXd RinvSqW = get_RinvSqX(chol_R_,W);  // n*n*a -> n x a
// //     MatrixXd WtRinvX = RinvSqW.transpose() * RinvSqX; // a*n*b -> a*b
// //     MatrixXd invSqAbXtRinvW = chol_A_beta.transpose().triangularView<Lower>().solve(WtRinvX.transpose()); // b*b*a -> b x a
// //
// //     MatrixXd A_alpha = Y_prec * (RinvSqW.transpose() * RinvSqW - invSqAbXtRinvW.transpose() * invSqAbXtRinvW);
// //     A_alpha.diagonal() += prior_prec_alpha;
// //
// //     LLT<MatrixXd> A_alpha_llt;
// //     A_alpha_llt.compute(A_alpha);
// //     MatrixXd chol_A_alpha = A_alpha_llt.matrixU();
// //
// //     VectorXd Rinvsqy = get_RinvSqX(chol_R_,y); // n*n*q -> n x 1;
// //     VectorXd XtRinvy = RinvSqX.transpose() * Rinvsqy; // b*n*1 >- b x 1
// //     VectorXd invSqAbXtRinvy = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy); // b*b*1 -> b*1
// //
// //     VectorXd WtSbinvy = RinvSqW.transpose() * Rinvsqy - invSqAbXtRinvW.transpose() * invSqAbXtRinvy;
// //
// //     alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(WtSbinvy) * Y_prec + randn_alpha;
// //     alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
// //     y_tilde = y - W * alpha;
// //   }
// //
// //   // Step 2 - sample Y_prec
// //   // We don't need to actually calculate Sigma_beta^{-1} directly.
// //   VectorXd RinvSqy = get_RinvSqX(chol_R_,y_tilde);
// //   VectorXd XtRinvy = RinvSqX.transpose() * RinvSqy;
// //   VectorXd prod1 = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy);
// //   double score = Y_prec_b0 + (RinvSqy.dot(RinvSqy) - prod1.dot(prod1))/2;
// //   Y_prec = rgamma_1/score;
// //
// //   // Step 3 - sample beta
// //   VectorXd XtRinvy_std_mu = XtRinvy*Y_prec + prior_prec_beta.asDiagonal()*prior_mean_beta;
// //   VectorXd beta = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy_std_mu) / sqrt(Y_prec) + randn_beta;
// //   beta = chol_A_beta.triangularView<Upper>().solve(beta) / sqrt(Y_prec);
// //
// //   VectorXd result(1+a+b);
// //   result << Y_prec,alpha,beta;
// //
// //   return(result);
// // }
// //
// //
// // VectorXd regression_sampler_v2(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b > n
// //     const Ref<const VectorXd>& y,           // nx1
// //     const MatrixXd& W,           // nxa
// //     const MatrixXd& X,           // nxm or nxb
// //     const VectorXd& prior_prec_alpha, // ax 1
// //     const Ref<const VectorXd>& prior_mean_beta,  // bx1
// //     const Ref<const VectorXd>& prior_prec_beta,  // bx1
// //     const SEXP chol_R_,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
// //     const MatrixXd& R,
// //     double Y_prec,                    // double
// //     const VectorXd& randn_alpha,
// //     const Ref<const VectorXd>& randn_beta,
// //     const Ref<const VectorXd>& randn_e,
// //     const double rgamma_1,
// //     const double Y_prec_b0
// // ){
// //   int n = y.size();
// //   int a = W.cols();
// //   int b = X.cols();
// //
// //   // Check inputs
// //   if(W.rows() != n) stop("Wrong dimension of W");
// //   if(X.rows() != n) stop("Wrong dimension of X");
// //   if(prior_prec_alpha.size() != a) stop("Wrong length of prior_prec_alpha");
// //   if(prior_mean_beta.size() != b) stop("Wrong length of prior_mean_beta");
// //   if(prior_prec_beta.size() != b) stop("Wrong length of prior_prec_beta");
// //   if(randn_alpha.size() != a) stop("Wrong length of randn_alpha");
// //   if(randn_beta.size() != b) stop("Wrong length of randn_beta");
// //
// //   // Calculate inverse of Sigma_beta
// //   MatrixXd DXt = prior_prec_beta.cwiseInverse().asDiagonal() * X.transpose();
// //   MatrixXd Sigma_beta = X * DXt + R;
// //   LDLT<MatrixXd> Sigma_beta_ldlt;
// //   Sigma_beta_ldlt.compute(Sigma_beta);
// //
// //   // Step 1
// //   VectorXd alpha(a);
// //   VectorXd y_tilde = y;
// //   if(a > 0) {
// //     // Sample alpha
// //     MatrixXd SbinvW = Sigma_beta_ldlt.solve(W);
// //     MatrixXd A_alpha = Y_prec * SbinvW.transpose() * W;
// //     A_alpha.diagonal() += prior_prec_alpha;
// //
// //     LLT<MatrixXd> A_alpha_llt;
// //     A_alpha_llt.compute(A_alpha);
// //     MatrixXd chol_A_alpha = A_alpha_llt.matrixU();
// //
// //     VectorXd WtSbinvy = SbinvW.transpose() * y;
// //     alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(WtSbinvy) * Y_prec + randn_alpha;
// //     alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
// //     y_tilde = y - W * alpha;
// //   }
// //
// //   // Step 2 - sample Y_prec
// //   VectorXd e2 = y_tilde.transpose() * Sigma_beta_ldlt.solve(y_tilde);
// //   double score = Y_prec_b0 + e2[0]/2;
// //   Y_prec = rgamma_1/score;
// //
// //   // Step 3 - sample beta
// //   // what about prior mean?
// //   VectorXd u = randn_beta.array() / (prior_prec_beta * Y_prec).cwiseSqrt().array() + prior_mean_beta.array();
// //   VectorXd v = sqrt(Y_prec) * X * u;
// //   if(Rf_isMatrix(chol_R_)) {
// //     Map<MatrixXd> chol_R = as<Map<MatrixXd> >(chol_R_);
// //     v += chol_R.transpose().triangularView<Lower>() * randn_e;
// //   } else{
// //     MSpMat chol_R = as<MSpMat>(chol_R_);
// //     v += chol_R.transpose().triangularView<Lower>() * randn_e;
// //   }
// //   VectorXd w = Sigma_beta_ldlt.solve(y_tilde * sqrt(Y_prec) - v);
// //   VectorXd beta = u + DXt * w / sqrt(Y_prec);
// //
// //   VectorXd result(1+a+b);
// //   result << Y_prec,alpha,beta;
// //
// //   return(result);
// // }
// //
// //
// // VectorXd regression_sampler_v3(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b > n > m
// //     const Ref<const VectorXd>& y,           // nx1
// //     const MatrixXd& W,           // nxa
// //     const MatrixXd& U,           // nxm or nxb
// //     const MatrixXd& V,           // mxb
// //     const VectorXd& prior_prec_alpha, // ax 1
// //     const Ref<const VectorXd>& prior_mean_beta,  // bx1
// //     const Ref<const VectorXd>& prior_prec_beta,  // bx1
// //     const SEXP chol_R_,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
// //     const MatrixXd& Rinv,
// //     const MatrixXd& RinvU,
// //     const MatrixXd& UtRinvU,
// //     double Y_prec,                    // double
// //     const VectorXd& randn_alpha,
// //     const Ref<const VectorXd>& randn_beta,
// //     const Ref<const VectorXd>& randn_e,
// //     const double rgamma_1,
// //     const double Y_prec_b0
// // ){
// //   int n = y.size();
// //   int a = W.cols();
// //   if(V.rows() != U.cols()) stop("Wrong dimensions of V");
// //   MatrixXd X = U*V;
// //   int b = X.cols();
// //
// //   // Check inputs
// //   if(W.rows() != n) stop("Wrong dimension of W");
// //   if(X.rows() != n) stop("Wrong dimension of X");
// //   if(prior_prec_alpha.size() != a) stop("Wrong length of prior_prec_alpha");
// //   if(prior_mean_beta.size() != b) stop("Wrong length of prior_mean_beta");
// //   if(prior_prec_beta.size() != b) stop("Wrong length of prior_prec_beta");
// //   if(randn_alpha.size() != a) stop("Wrong length of randn_alpha");
// //   if(randn_beta.size() != b) stop("Wrong length of randn_beta");
// //
// //   // Calculate inverse of Sigma_beta
// //   // MatrixXd Sigma_beta_inv;
// //   MatrixXd DVt = prior_prec_beta.cwiseInverse().asDiagonal() * V.transpose();
// //   if(RinvU.rows() != n) stop("Wrong dimensions of RinvU");
// //   MatrixXd inner = (V * DVt).inverse() + UtRinvU;
// //   LDLT<MatrixXd> inner_ldlt;
// //   inner_ldlt.compute(inner);
// //   // Sigma_beta_inv = Rinv - RinvU * inner.ldlt().solve(RinvU.transpose());  // Don't actually calculate this. Stay in mxm space
// //
// //   // Step 1
// //   VectorXd alpha(a);
// //   VectorXd y_tilde = y;
// //   if(a > 0) {
// //     // Sample alpha
// //     // MatrixXd SbinvW = Sigma_beta_inv * W;
// //     // MatrixXd A_alpha = Y_prec * SbinvW.transpose() * W;
// //     MatrixXd RinvW = Rinv * W;
// //     MatrixXd UtRinvW = U.transpose() * RinvW;
// //     MatrixXd A_alpha = Y_prec * (W.transpose() * RinvW - UtRinvW.transpose() * inner_ldlt.solve(UtRinvW));
// //     A_alpha.diagonal() += prior_prec_alpha;
// //
// //     LLT<MatrixXd> A_alpha_llt;
// //     A_alpha_llt.compute(A_alpha);
// //     MatrixXd chol_A_alpha = A_alpha_llt.matrixU();
// //
// //     // VectorXd WtSbinvy = SbinvW.transpose() * y;
// //     VectorXd UtRinvy = RinvU.transpose() * y;
// //     VectorXd WtSbinvy = RinvW.transpose() * y - UtRinvW.transpose() * inner_ldlt.solve(UtRinvy);
// //
// //     alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(WtSbinvy) * Y_prec + randn_alpha;
// //     alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
// //     y_tilde = y - W * alpha;
// //   }
// //
// //   // Step 2 - sample Y_prec
// //   // VectorXd e2 = y_tilde.transpose() * Sigma_beta_inv * y_tilde;
// //   VectorXd Rinv_y = Rinv * y_tilde;
// //   VectorXd UtRinvy = RinvU.transpose() * y_tilde;
// //   VectorXd e2 = y_tilde.transpose() * Rinv_y - UtRinvy.transpose() * inner_ldlt.solve(UtRinvy);
// //
// //   double score = Y_prec_b0 + e2[0]/2;
// //   Y_prec = rgamma_1/score;
// //
// //   // Step 3 - sample beta
// //   // what about prior mean?
// //   VectorXd u = randn_beta.array() / (prior_prec_beta * Y_prec).cwiseSqrt().array() + prior_mean_beta.array();
// //   VectorXd v = sqrt(Y_prec) * X * u;
// //   if(Rf_isMatrix(chol_R_)) {
// //     Map<MatrixXd> chol_R = as<Map<MatrixXd> >(chol_R_);
// //     v += chol_R.transpose().triangularView<Lower>() * randn_e;
// //   } else{
// //     MSpMat chol_R = as<MSpMat>(chol_R_);
// //     v += chol_R.transpose().triangularView<Lower>() * randn_e;
// //   }
// //   // VectorXd w = Sigma_beta_inv * (y_tilde * sqrt(Y_prec) - v);
// //   VectorXd e = y_tilde * sqrt(Y_prec) - v;
// //   VectorXd UtRinve = RinvU.transpose() * e;
// //   VectorXd w = Rinv * e - RinvU * inner_ldlt.solve(UtRinve);
// //
// //   VectorXd beta = u + DVt * (U.transpose() * w) / sqrt(Y_prec); //b*b*1 + b*n*1
// //
// //   VectorXd result(1+a+b);
// //   result << Y_prec,alpha,beta;
// //
// //   return(result);
// // }
// //
// //
// //
// // struct regression_sampler_worker : public RcppParallel::Worker {
// //   const int sampler;  // which sampler to use?
// //   const IntegerVector& trait_set;
// //   const Map<MatrixXd> Y;           // nx1
// //   const Map<MatrixXd> W_base;           // nxa
// //   const Rcpp::List W_list;           // nxa
// //   const Map<MatrixXd> X_U;   // could be X or U.
// //   const Map<MatrixXd> V;
// //   const MatrixXd RinvSqX;                // nxb
// //   const MatrixXd& C;                     // bxb
// //   SEXP chol_R_;                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
// //   const MatrixXd& R;
// //   const MatrixXd& Rinv;
// //   const MatrixXd& RinvU;
// //   const MatrixXd& UtRinvU;
// //   const Map<MatrixXd> prior_prec_alpha1; // a1 x p matrix of prior precisions for alpha1
// //   const VectorXd& prior_prec_alpha2;     // p-vector of precision of alpha2s for each trait
// //   const Map<MatrixXd> prior_mean_beta; // b x p matrix of prior means of beta
// //   const Map<MatrixXd> prior_prec_beta; // b x p matrix of prior precisions of beta
// //   double Y_prec_b0;
// //
// //   const MatrixXd& randn_alpha1;
// //   const std::vector<VectorXd>& randn_alpha2;
// //   const MatrixXd& randn_beta;
// //   const MatrixXd& randn_e;
// //   const VectorXd& rgamma_1;
// //
// //   MatrixXd& alpha1;
// //   Rcpp::List alpha2;
// //   MatrixXd& beta;
// //   VectorXd& Y_prec;
// //
// //   regression_sampler_worker(
// //     const int sampler,
// //     const Rcpp::IntegerVector& trait_set,
// //     const Map<MatrixXd> Y,           // nx1
// //     const Map<MatrixXd> W_base,           // nxa
// //     const Rcpp::List W_list,           // nxa
// //     const Map<MatrixXd> X_U,
// //     const Map<MatrixXd> V,
// //     const MatrixXd RinvSqX,                // nxb
// //     const MatrixXd& C,                     // bxb
// //     SEXP chol_R_,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
// //     const MatrixXd& R,
// //     const MatrixXd& Rinv,
// //     const MatrixXd& RinvU,
// //     const MatrixXd& UtRinvU,
// //     const Map<MatrixXd> prior_prec_alpha1, // a1 x p matrix of prior precisions for alpha1
// //     const VectorXd& prior_prec_alpha2,     // p-vector of precision of alpha2s for each trait
// //     const Map<MatrixXd> prior_mean_beta, // b x p matrix of prior means of beta
// //     const Map<MatrixXd> prior_prec_beta, // b x p matrix of prior precisions of beta
// //     double Y_prec_b0,
// //     const MatrixXd& randn_alpha1,
// //     const std::vector<VectorXd>& randn_alpha2,
// //     const MatrixXd& randn_beta,
// //     const MatrixXd& randn_e,
// //     const VectorXd& rgamma_1,
// //     MatrixXd& alpha1,
// //     Rcpp::List alpha2,
// //     MatrixXd& beta,
// //     VectorXd& Y_prec
// //   ):
// //     sampler(sampler), trait_set(trait_set),
// //     Y(Y), W_base(W_base), W_list(W_list), X_U(X_U), V(V), RinvSqX(RinvSqX),
// //     C(C), chol_R_(chol_R_), R(R), Rinv(Rinv), RinvU(RinvU), UtRinvU(UtRinvU),
// //     prior_prec_alpha1(prior_prec_alpha1), prior_prec_alpha2(prior_prec_alpha2), prior_mean_beta(prior_mean_beta), prior_prec_beta(prior_prec_beta),
// //     Y_prec_b0(Y_prec_b0),
// //     randn_alpha1(randn_alpha1), randn_alpha2(randn_alpha2), randn_beta(randn_beta), randn_e(randn_e), rgamma_1(rgamma_1),
// //     alpha1(alpha1), alpha2(alpha2), beta(beta), Y_prec(Y_prec)
// //     {}
// //
// //   void operator()(std::size_t begin, std::size_t end) {
// //     int n = Y.rows();
// //     int a1 = W_base.cols();
// //
// //     for(std::size_t i = begin; i < end; i++){
// //       int j = trait_set[i];
// //       MatrixXd W;
// //       int a;
// //       int a2 = 0;
// //       int b;
// //       VectorXd prior_prec_alpha;
// //       VectorXd randn_alpha;
// //       if(W_list.length() == 0) {
// //         W = W_base;
// //         a = a1;
// //         prior_prec_alpha = prior_prec_alpha1.col(j);
// //         randn_alpha = randn_alpha1.col(j);
// //       } else{
// //         Map<MatrixXd> W2 = as<Map<MatrixXd> >(W_list[j]);
// //         a2 = W2.cols();
// //         a = a1+a2;
// //         W = MatrixXd(n,a);
// //         W << W_base,W2;
// //         prior_prec_alpha = VectorXd(a);
// //         prior_prec_alpha.head(a1) = prior_prec_alpha1.col(j);
// //         prior_prec_alpha.tail(a2).array() = prior_prec_alpha2[j];
// //         randn_alpha = VectorXd(a);
// //         randn_alpha.head(a1) = randn_alpha1.col(j);
// //         randn_alpha.tail(a2) = randn_alpha2[j];
// //       }
// //
// //       VectorXd samples;
// //       if(sampler == 1) {
// //         b = RinvSqX.cols();
// //         samples = regression_sampler_v1(Y.col(j), W, RinvSqX, C, prior_prec_alpha, prior_mean_beta.col(j),
// //                                                prior_prec_beta.col(j), chol_R_, Y_prec[j], randn_alpha,
// //                                                randn_beta.col(j), rgamma_1[j],Y_prec_b0);
// //       } else if(sampler == 2) {
// //         b = X_U.cols();
// //         samples = regression_sampler_v2(Y.col(j), W, X_U, prior_prec_alpha, prior_mean_beta.col(j),
// //                                         prior_prec_beta.col(j), chol_R_, R, Y_prec[j], randn_alpha,
// //                                         randn_beta.col(j), randn_e.col(j),rgamma_1[j],Y_prec_b0);
// //       } else if(sampler == 3) {
// //         b = V.cols();
// //         samples = regression_sampler_v3(Y.col(j), W, X_U, V, prior_prec_alpha, prior_mean_beta.col(j),
// //                                         prior_prec_beta.col(j), chol_R_, Rinv, RinvU, UtRinvU, Y_prec[j], randn_alpha,
// //                                         randn_beta.col(j), randn_e.col(j),rgamma_1[j],Y_prec_b0);
// //
// //       } else {
// //         stop("sampler not implemented");
// //       }
// //
// //       // extract samples
// //       Y_prec[j] = samples[0];
// //       if(a1 > 0) alpha1.col(j) = samples.segment(1,a1);
// //       if(a2 > 0) alpha2[j] = samples.segment(1+a1,a2);
// //       if(b > 0) beta.col(j) = samples.tail(b);
// //     }
// //   }
// // };
// //
// //
// //
// // // [[Rcpp::export]]
// // Rcpp::IntegerVector which(Rcpp::LogicalVector x) {
// //   Rcpp::IntegerVector v = Rcpp::seq(0, x.size()-1);
// //   return v[x];
// // }
// //
// // // [[Rcpp::export]]
// // Rcpp::List regression_sampler_parallel(
// //   Map<MatrixXd> Y,               // n x p matrix of observations
// //   Map<MatrixXd> W_base,          // n x a1 matrix of W covariates common to all p. Can be NULL
// //   Rcpp::List W_list,             // p-list of n x a2 matrices of W covariates unique to each p. Can be NULL
// //   Map<MatrixXd> X,               // either X, a n x b matrix, or U, a n x m matrix. If U, then V must be non-NULL
// //   SEXP V_,                       // m x b matrix if X is U
// //   Rcpp::IntegerVector h2s_index, // p-vector of indices for appropriate V of each trait
// //   Rcpp::List chol_V_list,        // list of cholesky decompositions of V: RtR (each nxn). Can be either dense or sparse
// //   VectorXd Y_prec,               // p-vector of Y current precisions
// //   double Y_prec_a0,
// //   double Y_prec_b0,
// //   Map<MatrixXd> prior_prec_alpha1, // a1 x p matrix of prior precisions for alpha1
// //   VectorXd prior_prec_alpha2,     // p-vector of precision of alpha2s for each trait
// //   Map<MatrixXd> prior_mean_beta, // b x p matrix of prior means of beta
// //   Map<MatrixXd> prior_prec_beta, // b x p matrix of prior precisions of beta
// //   int grainSize) {
// //
// //   int n = Y.rows();
// //   int p = Y.cols();
// //
// //   // W_base
// //   if(W_base.rows() != n) stop("Wrong dimension of W_base");
// //   int a1 = W_base.cols();
// //
// //   // W_list
// //   if(W_list.size() > 0) {
// //     if(W_list.size() != p) stop("Wrong length of W_list");
// //   }
// //
// //   // X or U and V
// //   Map<MatrixXd> U = X;
// //   MatrixXd z = MatrixXd::Zero(0,0);
// //   Map<MatrixXd> V(z.data(),0,0);
// //   int b = X.cols();
// //   if(X.rows() != n) stop("Wrong dimension of X");
// //   if(Rf_isMatrix(V_)) {
// //     // Map<MatrixXd> V__ = as<Map<MatrixXd> >(V_);
// //     // new (&v) Map<MatrixXd> (V__,V__.rows(),V__.cols());
// //     new (&V) Map<MatrixXd> (as<Map<MatrixXd> >(V_));
// //     if(U.cols() != V.rows()) stop("X and V_ have incompatible dimensions");
// //     b = V.cols();
// //   }
// //
// //   // chol_V_list
// //   if(max(h2s_index) > chol_V_list.size()) {
// //     stop("max(h2s_index) > length(chol_V_list)");
// //   }
// //
// //   // priors
// //   if(Y_prec.size() != p) {
// //     stop("Wrong length of Y_prec");
// //   }
// //   if(prior_prec_alpha1.rows() != a1 || prior_prec_alpha1.cols() != p) stop("Wrong dimensions of prior_prec_alpha1");
// //   if(prior_prec_alpha2.size() != p) {
// //     stop("Wrong length of prior_prec_alpha2");
// //   }
// //   if(prior_mean_beta.rows() != b || prior_mean_beta.cols() != p) stop("Wrong dimensions of prior_mean_beta");
// //   if(prior_prec_beta.rows() != b || prior_prec_beta.cols() != p) stop("Wrong dimensions of prior_prec_beta");
// //
// //   // generate random numbers
// //   MatrixXd randn_alpha1 = rstdnorm_mat2(a1,p);
// //   std::vector<VectorXd> randn_alpha2;
// //   if(W_list.size() > 0){
// //     for(int i = 0; i < p; i++){
// //       Map<MatrixXd> W2 = as<Map<MatrixXd> >(W_list[i]);
// //       randn_alpha2.push_back(rstdnorm_mat2(W2.cols(),1));
// //     }
// //   }
// //   MatrixXd randn_beta = rstdnorm_mat2(b,p);
// //   MatrixXd randn_e;
// //   if(b > n) {
// //     randn_e = rstdnorm_mat2(n,p);
// //   }
// //   VectorXd rgamma_1 = as<VectorXd>(rgamma(p,Y_prec_a0 + n/2.0,1.0));
// //
// //   // Results structures
// //   MatrixXd alpha1(a1,p);
// //   Rcpp::List alpha2(W_list.size());
// //   MatrixXd beta(b,p);
// //
// //   // go through h2s indices and sample columns with same index as a set
// //   for(int i = min(h2s_index); i <= max(h2s_index); i++) {
// //     int h2_index = i;
// //     IntegerVector trait_set = which(h2s_index == h2_index);  // list of traits with same h2_index
// //
// //     if(trait_set.size() > 0){
// //       // prepare matrices for sampler
// //       MatrixXd RinvSqX, C, R, Rinv, RinvU, UtRinvU;
// //       SEXP chol_R_ = chol_V_list[h2_index - 1];
// //       int which_sampler;
// //       // Decide which sampler to use
// //       if(b <= n) {
// //         // use regression_sampler_v1
// //         which_sampler = 1;
// //         if(Rf_isMatrix(chol_R_)) {
// //           Map<MatrixXd> chol_R = as<Map<MatrixXd> >(chol_R_);
// //           RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X);
// //         } else{
// //           MSpMat chol_R = as<MSpMat>(chol_R_);
// //           RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X);
// //         }
// //         C = RinvSqX.transpose() * RinvSqX;
// //       }
// //       else if(V.cols() == 0) {
// //         // use regression_sampler_v2
// //         which_sampler = 2;
// //         if(Rf_isMatrix(chol_R_)) {
// //           Map<MatrixXd> chol_R = as<Map<MatrixXd> >(chol_R_);
// //           R = chol_R.transpose().triangularView<Lower>() * chol_R;
// //         } else{
// //           MSpMat chol_R = as<MSpMat>(chol_R_);
// //           R = chol_R.transpose().triangularView<Lower>() * chol_R;
// //         }
// //       } else {
// //         // use regression_sampler_v3
// //         which_sampler = 3;
// //         if(Rf_isMatrix(chol_R_)) {
// //           Map<MatrixXd> chol_R = as<Map<MatrixXd> >(chol_R_);
// //           Rinv = chol_R.triangularView<Upper>().solve(chol_R.transpose().triangularView<Lower>().solve(MatrixXd::Identity(n,n)));
// //           RinvU = chol_R.triangularView<Upper>().solve(chol_R.transpose().triangularView<Lower>().solve(U));
// //         } else{
// //           MSpMat chol_R = as<MSpMat>(chol_R_);
// //           Rinv = chol_R.triangularView<Upper>().solve(chol_R.transpose().triangularView<Lower>().solve(MatrixXd::Identity(n,n)));
// //           RinvU = chol_R.triangularView<Upper>().solve(chol_R.transpose().triangularView<Lower>().solve(U));
// //         }
// //         UtRinvU = U.transpose() * RinvU;
// //       }
// //       regression_sampler_worker sampler(which_sampler, trait_set,Y,W_base,W_list,X,V,RinvSqX,C,
// //                                         chol_R_,R,Rinv,RinvU,UtRinvU,
// //                                         prior_prec_alpha1,prior_prec_alpha2,prior_mean_beta,prior_prec_beta,Y_prec_b0,
// //                                         randn_alpha1,randn_alpha2,randn_beta,randn_e, rgamma_1,
// //                                         alpha1,alpha2,beta,Y_prec);
// //       RcppParallel::parallelFor(0,trait_set.size(),sampler,grainSize);
// //     }
// //   }
// //   return(Rcpp::List::create(
// //       Named("alpha1") = alpha1,
// //       Named("alpha2") = alpha2,
// //       Named("beta") = beta,
// //       Named("Y_prec") = Y_prec
// //            ));
// // }
// //
// //
// //
// // // // [[Rcpp::export()]]
// // // VectorXd regression_sampler_v2a(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b > n
// // //     VectorXd y,           // nx1
// // //     MatrixXd W,           // nxa
// // //     MatrixXd RinvSqX,                // nxb
// // //     VectorXd prior_prec_alpha, // ax 1
// // //     VectorXd prior_mean_beta,  // bx1
// // //     VectorXd prior_prec_beta,  // bx1
// // //     SEXP chol_R_,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
// // //     double Y_prec,                    // double
// // //     VectorXd randn_alpha,
// // //     VectorXd randn_beta,
// // //     VectorXd randn_e,
// // //     double rgamma_1,
// // //     double Y_prec_b0
// // // ){
// // //   int n = y.size();
// // //   int a = W.cols();
// // //   int b = RinvSqX.cols();
// // //
// // //   // Check inputs
// // //   if(W.rows() != n) stop("Wrong dimension of W");
// // //   if(RinvSqX.rows() != n) stop("Wrong dimension of X");
// // //   if(prior_prec_alpha.size() != a) stop("Wrong length of prior_prec_alpha");
// // //   if(prior_mean_beta.size() != b) stop("Wrong length of prior_mean_beta");
// // //   if(prior_prec_beta.size() != b) stop("Wrong length of prior_prec_beta");
// // //   if(randn_alpha.size() != a) stop("Wrong length of randn_alpha");
// // //   if(randn_beta.size() != b) stop("Wrong length of randn_beta");
// // //
// // //   // Calculate cholesky of Sigma_beta
// // //
// // //   // form LinvXDbXtLit + I -> LDLt -> Sbinv
// // //   MatrixXd RinvSqX_D = RinvSqX * prior_prec_beta.cwiseInverse().asDiagonal();
// // //   MatrixXd A_w = RinvSqX_D * RinvSqX.transpose();
// // //   A_w.diagonal().array() += 1.0;
// // //   LDLT<MatrixXd> A_w_ldlt;
// // //   A_w_ldlt.compute(A_w);
// // //   // Sigma_beta_inv = chol_R * A_w_ldlt.solve(chol_R.transpose())
// // //
// // //   // Step 1
// // //   VectorXd alpha(a);
// // //   VectorXd y_tilde = y;
// // //   if(a > 0) {
// // //     // Sample alpha
// // //     // Calculate A_alpha = Y_prec*W^T*Sigma_beta^{-1}*W + D_alpha^{-1}
// // //     // We don't need to actually calculate Sigma_beta^{-1} directly.
// // //     MatrixXd RinvsqW = get_RinvSqX(chol_R_,W); // n*n*a -> n x a;
// // //     VectorXd Rinvsqy = get_RinvSqX(chol_R_,y); // n*n*1 -> n x 1;
// // //     MatrixXd A_alpha = Y_prec * RinvsqW.transpose() * A_w_ldlt.solve(RinvsqW);
// // //     A_alpha.diagonal() += prior_prec_alpha;
// // //
// // //     LLT<MatrixXd> A_alpha_llt;
// // //     A_alpha_llt.compute(A_alpha);
// // //     MatrixXd chol_A_alpha = A_alpha_llt.matrixU();
// // //
// // //     VectorXd WtSbinvy = RinvsqW.transpose() * A_w_ldlt.solve(Rinvsqy);
// // //
// // //     alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(WtSbinvy) * Y_prec + randn_alpha;
// // //     alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
// // //     y_tilde = y - W * alpha;
// // //   }
// // //
// // //   // Step 2 - sample Y_prec
// // //   // We don't need to actually calculate Sigma_beta^{-1} directly.
// // //   VectorXd Rinvsqy = get_RinvSqX(chol_R_,y_tilde); // n*n*1 -> n x 1;
// // //   VectorXd e2 = Rinvsqy.transpose() * A_w_ldlt.solve(Rinvsqy);
// // //   double score = Y_prec_b0 + e2[0]/2;
// // //   Y_prec = rgamma_1/score;
// // //
// // //   // Step 3 - sample beta
// // //   VectorXd u = randn_beta.array() / (prior_prec_beta * Y_prec).cwiseSqrt().array();
// // //   u += prior_mean_beta;
// // //   VectorXd v = sqrt(Y_prec)*RinvSqX * u + randn_e;
// // //   VectorXd w = A_w_ldlt.solve(Rinvsqy * sqrt(Y_prec) - v);
// // //   VectorXd beta = u + RinvSqX_D.transpose() * w / sqrt(Y_prec);
// // //
// // //   // // actually, probably don't need to form Sigma_beta explicitly
// // //   // MatrixXd Sigma_beta;
// // //   // if(Rf_isMatrix(chol_R_)) {
// // //   //   Map<MatrixXd> chol_R = as<Map<MatrixXd> >(chol_R_);
// // //   //   Sigma_beta = chol_R.triangularView<Upper>() * A_w_ldlt.solve(chol_R.transpose().triangularView<Lower>());
// // //   // } else{
// // //   //   MSpMat chol_R = as<MSpMat>(chol_R_);
// // //   //   Sigma_beta = chol_R.triangularView<Upper>() * A_w_ldlt.solve(chol_R.transpose().triangularView<Lower>());
// // //   // }
// // //
// // //
// // //   VectorXd result(1+a+b);
// // //   result << Y_prec,alpha,beta;
// // //
// // //   return(result);
// // // }
// //
// //
// //
// // void mod_list(Rcpp::List l,int i){
// //   l[i] = i+1;
// // }
// //
// // // [[Rcpp::export()]]
// // Rcpp::List makeList(int n){
// //   Rcpp::List l(n);
// //   for(int i = 0; i < n; i++){
// //     mod_list(l,i);
// //   }
// //   return(l);
// // }
// //
// // // }
// //
// //

// //
// // // [[Rcpp::export()]]
// // NumericVector my_gamma(int n, NumericVector shape, NumericVector scale) {
// //   NumericVector result(n);
// //   if(shape.size() == 1) shape = NumericVector::create(n,shape[0]);
// //   if(shape.size() != n) stop("Wrong length of shape");
// //   for(int i = 0; i < n; i++){
// //     // int shape_i;
// //     // int scale_i;
// //     // if(shape.size() == 1) {
// //     //   shape_i = shape[0];
// //     // } else {
// //     //   shape_i = shape[i];
// //     // }
// //     // if(scale.size() == 1) {
// //     //   scale_i = scale[0];
// //     // } else {
// //     //   scale_i = scale[i];
// //     // }
// //     result[i] = Rcpp::rgamma(1,shape[i],1.0)[0];
// //   }
// //   if(scale.size() == 1) return(result * scale[0]);
// //   if(scale.size() == n) return(result * scale);
// //   stop("Wrong length of scale");
// // }
