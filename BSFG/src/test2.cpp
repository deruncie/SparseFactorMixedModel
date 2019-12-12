// #include <math.h>
// #include <iostream>
// #include <RcppEigen.h>
// // [[Rcpp::plugins(openmp)]]
// #include <omp.h>
// #include "BSFG_types.h"
//
// using namespace Rcpp;
// using namespace Eigen;
//
//
// VectorXd cumprod2(const VectorXd& x) {
//   int n = x.size();
//   VectorXd res(n);
//   res[0] = x[0];
//   if(n > 1) {
//     for(int i = 1; i < n; i++){
//       res[i] = res[i-1]*x[i];
//     }
//   }
//   return(res);
// }
//
// // [[Rcpp::export()]]
// Rcpp::List sample_tau2_delta_c_Eigen_v3(
//     double tau2,
//     double xi,
//     VectorXd delta,
//     Map<VectorXd> scores,
//     double tau_0,
//     double delta_shape,
//     double delta_scale,  // shape and scale for inverse gamma distribution (shape and rate for Gamma)
//     int p,
//     int times
// ) {
//
//   int K = scores.size();
//   if(delta.size() != K) stop("Wrong size of delta");
//   double shape;
//   double scale;
//   VectorXd cumprod_delta = cumprod2(delta);
//   for(int i = 0; i < times; i++){
//     // sample tau2
//     shape = (p*K + 1)/2.0;
//     scale = 1.0/xi + cumprod_delta.cwiseInverse().dot(scores);
//     tau2 = 1.0/R::rgamma(shape,1.0/scale);
//
//     // sample xi
//     shape = 1.0;
//     scale = 1.0/(tau_0*tau_0) + 1.0/tau2;
//     xi = 1.0/R::rgamma(shape,1.0/scale);
//
//     for(int h = 1; h < K; h++) {
//       // delta_h
//       shape = delta_shape + p*(K-h)/2.0;
//       scale = delta_scale + delta(h) * cumprod_delta.tail(K-h).cwiseInverse().dot(scores.tail(K-h)) / tau2;
//       delta[h] = 1.0/R::rgamma(shape,1.0/scale);
//       cumprod_delta = cumprod2(delta);
//     }
//   }
//
//   return(Rcpp::List::create(Named("tau2") = tau2,
//                             Named("xi") = xi,
//                             Named("delta") = delta));
// }
//
// struct R_matrix2 {
//   Map<MatrixXd> dense;
//   MSpMat sparse;
//   bool isDense;
//   R_matrix2() {
//     SpMat null_s = null_d.sparseView();
//     MSpMat M_null_s(0,0,0,null_s.outerIndexPtr(),null_s.innerIndexPtr(),null_s.valuePtr());
//   }
//   R_matrix2(Map<MatrixXd> dense_, MSpMat sparse_,bool isDense_) : dense(dense_), sparse(sparse_), isDense(isDense_) {}
// };
//
// void test(SEXP Z_) {
//   R_matrix2 Z;
//   MatrixXd null_d = MatrixXd::Zero(0,0);
//   if(Rf_isMatrix(Z_)){
//     Map<MatrixXd> Zm = as<Map<MatrixXd> >(Z_);
//     SpMat null_s = null_d.sparseView();
//     MSpMat M_null_s(0,0,0,null_s.outerIndexPtr(),null_s.innerIndexPtr(),null_s.valuePtr());
//     Z = new R_matrix2(Zm,M_null_s,true);
//     // r = Zm.cols();
//   } else{
//     MSpMat Zm = as<MSpMat>(Z_);
//     Map<MatrixXd> M_null_d(null_d.data(),0,0);
//     Z = new R_matrix2(M_null_d,Zm,false);
//     // r = Zm.cols();
//   }
//   Rcout << Z.isDense << std::endl;
// }

// Rcpp::IntegerVector which2(Rcpp::LogicalVector x) {
//   Rcpp::IntegerVector v = Rcpp::seq(0, x.size()-1);
//   return v[x];
// }
// // [[Rcpp::export()]]
// MatrixXd rstdnorm_mat2(int n,int p) {  // returns nxp matrix
//   VectorXd X_vec(n*p);
//   for(int i = 0; i < n*p; i++){
//     X_vec[i] = ziggr.norm();
//   }
//   MatrixXd X_mat = Map<MatrixXd>(X_vec.data(),n,p);
//   return(X_mat);
// }
//
// void load_R_matrices_list2(const Rcpp::List X_list, std::vector<R_matrix>& X_vector){
//   // null_matrices
//   MatrixXd null_d = MatrixXd::Zero(0,0);
//   Map<MatrixXd> M_null_d(null_d.data(),0,0);
//   SpMat null_s = null_d.sparseView();
//   MSpMat M_null_s(0,0,0,null_s.outerIndexPtr(),null_s.innerIndexPtr(),null_s.valuePtr());
//
//   int p = X_list.size();
//   X_vector.reserve(p);
//   for(int i = 0; i < p; i++){
//     SEXP Xi_ = X_list[i];
//     if(Rf_isMatrix(Xi_)){
//       Map<MatrixXd> Xi = as<Map<MatrixXd> >(Xi_);
//       R_matrix Xim(Xi,M_null_s,true);
//       X_vector.push_back(Xim);
//     } else{
//       MSpMat Xi = as<MSpMat>(Xi_);
//       R_matrix Xim(M_null_d,Xi,false);
//       X_vector.push_back(Xim);
//     }
//   }
// }
//
//
//
//
// VectorXd sample_MME_single_diagR2(
//     const Ref<const VectorXd>& y,           // nx1
//     MSpMat Z,             // nxr dgCMatrix
//     MSpMat chol_ZtZ_Kinv,       // rxr CsparseMatrix upper triangular: chol(ZtRinvZ + diag(Kinv))
//     double tot_Eta_prec,   // double
//     double pe,            // double
//     VectorXd randn_theta  // rx1
// ){
//   VectorXd b = Z.transpose() * y * pe;
//   b = chol_ZtZ_Kinv.transpose().triangularView<Lower>().solve(b / sqrt(tot_Eta_prec));
//   b += randn_theta;
//   b = chol_ZtZ_Kinv.triangularView<Upper>().solve(b / sqrt(tot_Eta_prec));
//   return(b);
// }
//
//
// // samples random effects from model:
// // Y = ZU + E
// // U[,j] ~ N(0,1/tot_Eta_prec[j] * h2[j] * K)
// // E[,j] ~ N(0,1/tot_Eta_prec[j] * (1-h2[j]) * I_n)
// // For complete data, ie no missing obs.
// // [[Rcpp::export()]]
// MatrixXd sample_MME_ZKZts_c2(
//     Map<MatrixXd> Y,                    // nxp
//     MSpMat Z,                           // nxr
//     Map<VectorXd> tot_Eta_prec,         // px1
//     Rcpp::List chol_ZtZ_Kinv_list_,      // List or R st RtR = ZtZ_Kinv
//     Map<MatrixXd> h2s,                  // n_RE x p
//     VectorXi h2s_index                 // px1
// ) {
//
//   int p = Y.cols();
//   int r = Z.cols();
//
//   MatrixXd randn_theta = rstdnorm_mat2(r,p);
//
//   std::vector<R_matrix> chol_ZtZ_Kinv_list;
//   load_R_matrices_list2(chol_ZtZ_Kinv_list_, chol_ZtZ_Kinv_list);
//
//   MatrixXd U(r,p);
//   ArrayXd h2_e = 1.0 - h2s.colwise().sum().array();
//   ArrayXd pes = tot_Eta_prec.array() / h2_e.array();
//
//   // #pragma omp parallel for
//   // for(int i = min(h2s_index); i <= max(h2s_index); i++) {
//   //   int h2_index = i;
//   //   VectorXi trait_set = as<VectorXi>(which(h2s_index == h2_index));  // list of traits with same h2_index
//   // }
//
//
// #pragma omp parallel for
//   for(std::size_t j = 0; j < p; j++){
//     int h2_index = h2s_index[j] - 1;
//     // ZtZ_Kinv needs to be scaled by tot_Eta_prec[j].
//     U.col(j) = sample_MME_single_diagR2(Y.col(j), Z, chol_ZtZ_Kinv_list[h2_index].sparse, tot_Eta_prec[j], pes[j],randn_theta.col(j));
//   }
//
//   return(U);
// }
//
//
//
// VectorXd sample_MME_single_diagR3(
//     VectorXd y,           // nx1
//     MSpMat Z,             // nxr dgCMatrix
//     MSpMat chol_ZtZ_Kinv,       // rxr CsparseMatrix upper triangular: chol(ZtRinvZ + diag(Kinv))
//     double tot_Eta_prec,   // double
//     double pe,            // double
//     VectorXd randn_theta  // rx1
// ){
//   VectorXd b = Z.transpose() * y * pe;
//   b = chol_ZtZ_Kinv.transpose().triangularView<Lower>().solve(b / sqrt(tot_Eta_prec));
//   b += randn_theta;
//   b = chol_ZtZ_Kinv.triangularView<Upper>().solve(b / sqrt(tot_Eta_prec));
//   return(b);
// }
//
//
// // samples random effects from model:
// // Y = ZU + E
// // U[,j] ~ N(0,1/tot_Eta_prec[j] * h2[j] * K)
// // E[,j] ~ N(0,1/tot_Eta_prec[j] * (1-h2[j]) * I_n)
// // For complete data, ie no missing obs.
// // [[Rcpp::export()]]
// MatrixXd sample_MME_ZKZts_c3(
//     Map<MatrixXd> Y,                    // nxp
//     MSpMat Z,                           // nxr
//     Map<VectorXd> tot_Eta_prec,         // px1
//     MSpMat chol_ZtZ_Kinv,      // List or R st RtR = ZtZ_Kinv
//     Map<MatrixXd> h2s,                  // n_RE x p
//     VectorXi h2s_index                 // px1
// ) {
//
//   int p = Y.cols();
//   int r = Z.cols();
//
//   MatrixXd randn_theta = rstdnorm_mat2(r,p);
//
//   // std::vector<R_matrix> chol_ZtZ_Kinv_list;
//   // load_R_matrices_list2(chol_ZtZ_Kinv_list_, chol_ZtZ_Kinv_list);
//
//   // MatrixXd U(r,p);
//   ArrayXd h2_e = 1.0 - h2s.colwise().sum().array();
//   ArrayXd pes = tot_Eta_prec.array() / h2_e.array();
//
//   // #pragma omp parallel for
//   // for(int i = min(h2s_index); i <= max(h2s_index); i++) {
//   //   int h2_index = i;
//   //   VectorXi trait_set = as<VectorXi>(which(h2s_index == h2_index));  // list of traits with same h2_index
//   // }
//
//   MatrixXd b = Z.transpose() * Y * pes.matrix().asDiagonal();
//   b = chol_ZtZ_Kinv.transpose().triangularView<Lower>().solve(b * tot_Eta_prec.cwiseInverse().cwiseSqrt().asDiagonal());
//   b += randn_theta;
//   b = chol_ZtZ_Kinv.triangularView<Upper>().solve(b * tot_Eta_prec.cwiseInverse().cwiseSqrt().asDiagonal());
//   return(b);
//
//
// // #pragma omp parallel for
// //   for(std::size_t j = 0; j < p; j++){
// //     int h2_index = h2s_index[j] - 1;
// //     // ZtZ_Kinv needs to be scaled by tot_Eta_prec[j].
// //     U.col(j) = sample_MME_single_diagR3(Y.col(j), Z, chol_ZtZ_Kinv_list[h2_index].sparse, tot_Eta_prec[j], pes[j],randn_theta.col(j));
// //   }
//
//   // return(U);
// }


//
// // [[Rcpp::export()]]
// MatrixXd my_mult(Map<MatrixXd> X,Map<MatrixXd> Y,int n_threads) {
//   Eigen::setNbThreads(n_threads);
//   Rcout << Eigen::nbThreads( ) << std::endl;
//   return(X * Y);
// }


//
// // Attempt to use MatrixXf instead of MatriXd. Still slower because of data copy.
// #include <math.h>
// #include <iostream>
// #include <RcppEigen.h>
// #include <RcppParallel.h>
// // [[Rcpp::depends("RcppParallel")]]
// #include <ZigguratR.h>
// // [[Rcpp::depends("RcppZiggurat")]]
//
// using Eigen::Map;               	      // 'Eigen::Maps' rather than copies
// using Eigen::MatrixXd;                  // variable size matrix, double precision
// using Eigen::VectorXd;                  // variable size vector, double precision
// using Eigen::VectorXi;                  // variable size vector, int
// using Eigen::ArrayXXd;                  // variable size array, double precision
// using Eigen::ArrayXd;                  // variable size array, double precision
// using Eigen::Upper;
// using Eigen::Lower;
// typedef Eigen::SparseMatrix<double> SpMat;
// typedef Eigen::Map<SpMat> MSpMat;
//
// using namespace Rcpp;
// using namespace RcppParallel;
// using namespace Eigen;
//
// static Ziggurat::R::ZigguratR ziggr;
//
// VectorXd find_candidate_states(MatrixXd, double, int);
// // MatrixXd uncorrelated_prec_mat(VectorXd,VectorXd,VectorXd);
// MatrixXd rstdnorm_mat(int n,int p);
// struct R_matrix {
//   Map<MatrixXd> dense;
//   MSpMat sparse;
//   bool isDense;
//   R_matrix(Map<MatrixXd> dense_, MSpMat sparse_,bool isDense_) : dense(dense_), sparse(sparse_), isDense(isDense_) {}
// };
// void load_R_matrices_list(const Rcpp::List X_list, std::vector<R_matrix>& X_vector);
//
//
//
// // #include "BSFG_types.h"
// // // #include <Rinternals.h>
// // // #include <R.h>
// //
// // [[Rcpp::depends(RcppEigen)]]
// using namespace Eigen;
// using namespace RcppParallel;
//
// // [[Rcpp::plugins(openmp)]]
// #include <omp.h>
//
//
// // [[Rcpp::export]]
// void set_nthreads(int threads) {
//   if ( threads > 0 )
//     omp_set_num_threads( threads );
//   REprintf("Number of threads=%i\\n", omp_get_max_threads());
// }
//
// // #define FLOAT(x) ((float*) INTEGER(x))
// //
// // // [[Rcpp::export()]]
// // SEXP R_add1(SEXP x_) {
// //   SEXP ret;
// //   PROTECT(ret = Rf_allocVector(INTSXP , 1));
// //
// //   float *x = FLOAT(x_);
// //   FLOAT(ret)[0] = x[0] + 1.0f;
// //
// //   UNPROTECT(1);
// //   return ret;
// // }
// //
// // // [[Rcpp::export()]]
// // MatrixXf R_add1b(SEXP x_) {
// //   MatrixXf Xf = as<MatrixXf>(FLOAT(x_));
// //   Xf(0,0) += 1.0f;
// //
// //   return Xf;
// // }
// // [[Rcpp::export()]]
// MatrixXd rstdnorm_mata(int n,int p) {  // returns nxp matrix
//   VectorXd X_vec(n*p);
//   for(int i = 0; i < n*p; i++){
//     X_vec[i] = ziggr.norm();
//   }
//   MatrixXd X_mat = Map<MatrixXd>(X_vec.data(),n,p);
//   return(X_mat);
// }
//
//
// void load_R_matrices_lista(const Rcpp::List X_list, std::vector<R_matrix>& X_vector){
//   // null_matrices
//   MatrixXd null_d = MatrixXd::Zero(0,0);
//   Map<MatrixXd> M_null_d(null_d.data(),0,0);
//   SpMat null_s = null_d.sparseView();
//   MSpMat M_null_s(0,0,0,null_s.outerIndexPtr(),null_s.innerIndexPtr(),null_s.valuePtr());
//
//   int p = X_list.size();
//   X_vector.reserve(p);
//   for(int i = 0; i < p; i++){
//     SEXP Xi_ = X_list[i];
//     if(Rf_isMatrix(Xi_)){
//       Map<MatrixXd> Xi = as<Map<MatrixXd> >(Xi_);
//       R_matrix Xim(Xi,M_null_s,true);
//       X_vector.push_back(Xim);
//     } else{
//       MSpMat Xi = as<MSpMat>(Xi_);
//       R_matrix Xim(M_null_d,Xi,false);
//       X_vector.push_back(Xim);
//     }
//   }
// }
//
// MatrixXd get_RinvSqXa(const R_matrix& chol_R, MatrixXd X){
//   MatrixXd RinvSqX;
//   if(chol_R.isDense) {
//     RinvSqX = chol_R.dense.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
//   } else{
//     RinvSqX = chol_R.sparse.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
//   }
//   return(RinvSqX);
// }
// Rcpp::IntegerVector whicha(Rcpp::LogicalVector x) {
//   Rcpp::IntegerVector v = Rcpp::seq(0, x.size()-1);
//   return v[x];
// }
//
//
// VectorXd regression_sampler_v1a(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b < n
//     const Ref<const VectorXd>& y,           // nx1
//     const MatrixXd& W,           // nxa
//     const MatrixXd& RinvSqX,                // nxb
//     const MatrixXd& C,                     // bxb
//     const VectorXd& prior_prec_alpha, // ax 1
//     const Ref<const VectorXd>& prior_mean_beta,  // bx1
//     const Ref<const VectorXd>& prior_prec_beta,  // bx1
//     const R_matrix& chol_R,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
//     double Y_prec,                    // double
//     const VectorXd& randn_alpha,
//     const Ref<const VectorXd>& randn_beta,
//     const double rgamma_1,
//     const double Y_prec_b0
// ){
//   int n = y.size();
//   int a = W.cols();
//   int b = RinvSqX.cols();
//
//   // Check inputs
//   if(W.rows() != n) stop("Wrong dimension of W");
//   if(RinvSqX.rows() != n) stop("Wrong dimension of X");
//   if(C.rows() != b || C.cols() != b) stop("Wrong dimension of C");
//   if(prior_prec_alpha.size() != a) stop("Wrong length of prior_prec_alpha");
//   if(prior_mean_beta.size() != b) stop("Wrong length of prior_mean_beta");
//   if(prior_prec_beta.size() != b) stop("Wrong length of prior_prec_beta");
//   if(randn_alpha.size() != a) stop("Wrong length of randn_alpha");
//   if(randn_beta.size() != b) stop("Wrong length of randn_beta");
//
//   // Calculate cholesky of A_beta
//   // C = Xt(RtR)^-1X
//   MatrixXd C_beta = C;
//   C_beta.diagonal() += prior_prec_beta;
//   LLT<MatrixXd> A_beta_llt;
//   A_beta_llt.compute(C_beta);
//   MatrixXd chol_A_beta = A_beta_llt.matrixU();
//   // Y_prec * chol_A_beta^\T * chol_A_beta = A_beta
//
//   // Step 1
//   VectorXd alpha(a);
//   VectorXd y_tilde = y;
//   if(a > 0) {
//     // Sample alpha
//     // Calculate A_alpha = Y_prec*W^T*Sigma_beta^{-1}*W + D_alpha^{-1}
//     // We don't need to actually calculate Sigma_beta^{-1} directly.
//     MatrixXd RinvSqW = get_RinvSqXa(chol_R,W);  // n*n*a -> n x a
//     MatrixXd WtRinvX = RinvSqW.transpose() * RinvSqX; // a*n*b -> a*b
//     MatrixXd invSqAbXtRinvW = chol_A_beta.transpose().triangularView<Lower>().solve(WtRinvX.transpose()); // b*b*a -> b x a
//
//     MatrixXd A_alpha = Y_prec * (RinvSqW.transpose() * RinvSqW - invSqAbXtRinvW.transpose() * invSqAbXtRinvW);
//     A_alpha.diagonal() += prior_prec_alpha;
//
//     VectorXd Rinvsqy = get_RinvSqXa(chol_R,y); // n*n*q -> n x 1;
//     VectorXd XtRinvy = RinvSqX.transpose() * Rinvsqy; // b*n*1 >- b x 1
//     VectorXd invSqAbXtRinvy = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy); // b*b*1 -> b*1
//
//     VectorXd WtSbinvy = RinvSqW.transpose() * Rinvsqy - invSqAbXtRinvW.transpose() * invSqAbXtRinvy;
//
//     LLT<MatrixXd> A_alpha_llt;
//     A_alpha_llt.compute(A_alpha);
//     MatrixXd chol_A_alpha = A_alpha_llt.matrixU();
//
//     alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(WtSbinvy) * Y_prec + randn_alpha;
//     alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
//
//     y_tilde = y - W * alpha;
//   }
//
//   // Step 2 - sample Y_prec
//   // We don't need to actually calculate Sigma_beta^{-1} directly.
//   VectorXd RinvSqy = get_RinvSqXa(chol_R,y_tilde);
//   VectorXd XtRinvy = RinvSqX.transpose() * RinvSqy;
//   VectorXd prod1 = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy);
//   double score = Y_prec_b0 + (RinvSqy.dot(RinvSqy) - prod1.dot(prod1))/2;
//   Y_prec = rgamma_1/score;
//
//   // Step 3 - sample beta
//   VectorXd XtRinvy_std_mu = XtRinvy*Y_prec + prior_prec_beta.asDiagonal()*prior_mean_beta;
//   VectorXd beta = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy_std_mu) / sqrt(Y_prec) + randn_beta;
//   beta = chol_A_beta.triangularView<Upper>().solve(beta) / sqrt(Y_prec);
//
//   VectorXd result(1+a+b);
//   result << Y_prec,alpha,beta;
//
//   return(result);
// }
//
//
// VectorXd regression_sampler_v2a(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b > n
//     const Ref<const VectorXd>& y,           // nx1
//     const MatrixXd& W,           // nxa
//     const MatrixXd& X,           // nxm or nxb
//     const VectorXd& prior_prec_alpha, // ax 1
//     const Ref<const VectorXd>& prior_mean_beta,  // bx1
//     const Ref<const VectorXd>& prior_prec_beta,  // bx1
//     const R_matrix& chol_R,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
//     const MatrixXd& R,
//     double Y_prec,                    // double
//     const VectorXd& randn_alpha,
//     const Ref<const VectorXd>& randn_beta,
//     const Ref<const VectorXd>& randn_e,
//     const double rgamma_1,
//     const double Y_prec_b0
// ){
//   int n = y.size();
//   int a = W.cols();
//   int b = X.cols();
//
//   // Check inputs
//   if(W.rows() != n) stop("Wrong dimension of W");
//   if(X.rows() != n) stop("Wrong dimension of X");
//   if(prior_prec_alpha.size() != a) stop("Wrong length of prior_prec_alpha");
//   if(prior_mean_beta.size() != b) stop("Wrong length of prior_mean_beta");
//   if(prior_prec_beta.size() != b) stop("Wrong length of prior_prec_beta");
//   if(randn_alpha.size() != a) stop("Wrong length of randn_alpha");
//   if(randn_beta.size() != b) stop("Wrong length of randn_beta");
//
//   // Calculate inverse of Sigma_beta
//   MatrixXd DXt = prior_prec_beta.cwiseInverse().asDiagonal() * X.transpose();
//   MatrixXd Sigma_beta = X * DXt + R;
//   LDLT<MatrixXd> Sigma_beta_ldlt;
//   Sigma_beta_ldlt.compute(Sigma_beta);
//
//   // Step 1
//   VectorXd alpha(a);
//   VectorXd y_tilde = y;
//   if(a > 0) {
//     // Sample alpha
//     MatrixXd SbinvW = Sigma_beta_ldlt.solve(W);
//     MatrixXd A_alpha = Y_prec * SbinvW.transpose() * W;
//     A_alpha.diagonal() += prior_prec_alpha;
//
//     LLT<MatrixXd> A_alpha_llt;
//     A_alpha_llt.compute(A_alpha);
//     MatrixXd chol_A_alpha = A_alpha_llt.matrixU();
//
//     VectorXd WtSbinvy = SbinvW.transpose() * y;
//     alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(WtSbinvy) * Y_prec + randn_alpha;
//     alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
//     y_tilde = y - W * alpha;
//   }
//
//   // Step 2 - sample Y_prec
//   VectorXd e2 = y_tilde.transpose() * Sigma_beta_ldlt.solve(y_tilde);
//   double score = Y_prec_b0 + e2[0]/2;
//   Y_prec = rgamma_1/score;
//
//   // Step 3 - sample beta
//   // what about prior mean?
//   VectorXd u = randn_beta.array() / (prior_prec_beta * Y_prec).cwiseSqrt().array() + prior_mean_beta.array();
//   VectorXd v = sqrt(Y_prec) * X * u;
//   if(chol_R.isDense) {
//     v += chol_R.dense.transpose().triangularView<Lower>() * randn_e;
//   } else{
//     v += chol_R.sparse.transpose().triangularView<Lower>() * randn_e;
//   }
//   VectorXd w = Sigma_beta_ldlt.solve(y_tilde * sqrt(Y_prec) - v);
//   VectorXd beta = u + DXt * w / sqrt(Y_prec);
//
//   VectorXd result(1+a+b);
//   result << Y_prec,alpha,beta;
//
//   return(result);
// }
//
//
// VectorXd regression_sampler_v3a(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b > n > m
//     const Ref<const VectorXd>& y,           // nx1
//     const MatrixXd& W,           // nxa
//     const MatrixXd& U,           // nxm or nxb
//     const MatrixXd& V,           // mxb
//     const VectorXd& prior_prec_alpha, // ax 1
//     const Ref<const VectorXd>& prior_mean_beta,  // bx1
//     const Ref<const VectorXd>& prior_prec_beta,  // bx1
//     const R_matrix& chol_R,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
//     const MatrixXd& Rinv,
//     const MatrixXd& RinvU,
//     const MatrixXd& UtRinvU,
//     double Y_prec,                    // double
//     const VectorXd& randn_alpha,
//     const Ref<const VectorXd>& randn_beta,
//     const Ref<const VectorXd>& randn_e,
//     const double rgamma_1,
//     const double Y_prec_b0
// ){
//   int n = y.size();
//   int a = W.cols();
//   if(V.rows() != U.cols()) stop("Wrong dimensions of V");
//   MatrixXd X = U*V;
//   int b = X.cols();
//
//   // Check inputs
//   if(W.rows() != n) stop("Wrong dimension of W");
//   if(X.rows() != n) stop("Wrong dimension of X");
//   if(prior_prec_alpha.size() != a) stop("Wrong length of prior_prec_alpha");
//   if(prior_mean_beta.size() != b) stop("Wrong length of prior_mean_beta");
//   if(prior_prec_beta.size() != b) stop("Wrong length of prior_prec_beta");
//   if(randn_alpha.size() != a) stop("Wrong length of randn_alpha");
//   if(randn_beta.size() != b) stop("Wrong length of randn_beta");
//   if(randn_e.size() != n) stop("Wrong length of randn_e");
//
//   // Calculate inverse of Sigma_beta
//   // MatrixXd Sigma_beta_inv;
//   // Using Ainv - Ainv * U * (I + BVAinvU)inv * BVAinv in case B = VDVt is singular
//   MatrixXd DVt = prior_prec_beta.cwiseInverse().asDiagonal() * V.transpose();
//   MatrixXd VDVt = V * DVt;
//   if(RinvU.rows() != n) stop("Wrong dimensions of RinvU");
//   MatrixXd inner = VDVt * UtRinvU;
//   inner.diagonal().array() += 1.0;
//   LDLT<MatrixXd> inner_ldlt;
//   inner_ldlt.compute(inner);
//   // Sigma_beta_inv = Rinv - RinvU * inner.ldlt().solve(VDVt * RinvU.transpose());  // Don't actually calculate this. Stay in mxm space
//
//   // Step 1
//   VectorXd alpha(a);
//   VectorXd y_tilde = y;
//   if(a > 0) {
//     // Sample alpha
//     // MatrixXd SbinvW = Sigma_beta_inv * W;
//     // MatrixXd A_alpha = Y_prec * SbinvW.transpose() * W;
//     MatrixXd RinvW = Rinv * W;
//     MatrixXd UtRinvW = U.transpose() * RinvW;
//     MatrixXd A_alpha = Y_prec * (W.transpose() * RinvW - UtRinvW.transpose() * inner_ldlt.solve(VDVt * UtRinvW));
//     A_alpha.diagonal() += prior_prec_alpha;
//
//     // VectorXd WtSbinvy = SbinvW.transpose() * y;
//     VectorXd UtRinvy = RinvU.transpose() * y;
//     VectorXd WtSbinvy = RinvW.transpose() * y - UtRinvW.transpose() * inner_ldlt.solve(VDVt * UtRinvy);
//
//     LLT<MatrixXd> A_alpha_llt;
//     A_alpha_llt.compute(A_alpha);
//     MatrixXd chol_A_alpha = A_alpha_llt.matrixU();
//
//     alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(WtSbinvy) * Y_prec + randn_alpha;
//     alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
//
//     y_tilde = y - W * alpha;
//   }
//
//   // Step 2 - sample Y_prec
//   // VectorXd e2 = y_tilde.transpose() * Sigma_beta_inv * y_tilde;
//   VectorXd Rinv_y = Rinv * y_tilde;
//   VectorXd UtRinvy = RinvU.transpose() * y_tilde;
//   VectorXd e2 = y_tilde.transpose() * Rinv_y - UtRinvy.transpose() * inner_ldlt.solve(VDVt * UtRinvy);
//
//   double score = Y_prec_b0 + e2[0]/2;
//   Y_prec = rgamma_1/score;
//
//   // Step 3 - sample beta
//   // what about prior mean?
//   VectorXd u = randn_beta.array() / (prior_prec_beta * Y_prec).cwiseSqrt().array() + prior_mean_beta.array();
//   VectorXd v = sqrt(Y_prec) * X * u;
//   if(chol_R.isDense) {
//     v += chol_R.dense.transpose().triangularView<Lower>() * randn_e;
//   } else{
//     v += chol_R.sparse.transpose().triangularView<Lower>() * randn_e;
//   }
//   // VectorXd w = Sigma_beta_inv * (y_tilde * sqrt(Y_prec) - v);
//   VectorXd e = y_tilde * sqrt(Y_prec) - v;
//   VectorXd UtRinve = RinvU.transpose() * e;
//   VectorXd w = Rinv * e - RinvU * inner_ldlt.solve(VDVt * UtRinve);
//
//   VectorXd beta = u + DVt * (U.transpose() * w) / sqrt(Y_prec); //b*b*1 + b*n*1
//
//   VectorXd result(1+a+b);
//   result << Y_prec,alpha,beta;
//
//   return(result);
// }
//
//
//
// struct regression_sampler_workera : public RcppParallel::Worker {
//   const int sampler;  // whicha sampler to use?
//   const VectorXi& trait_set;
//   const Map<MatrixXd> Y;           // nx1
//   const Map<MatrixXd> W_base;           // nxa
//   const std::vector<R_matrix>& W_list;           // nxa
//   const Map<MatrixXd> X_U;   // could be X or U.
//   const Map<MatrixXd> V;
//   const MatrixXd RinvSqX;                // nxb
//   const MatrixXd& C;                     // bxb
//   const R_matrix& chol_R;                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
//   const MatrixXd& R;
//   const MatrixXd& Rinv;
//   const MatrixXd& RinvU;
//   const MatrixXd& UtRinvU;
//   const Map<MatrixXd> prior_prec_alpha1; // a1 x p matrix of prior precisions for alpha1
//   const VectorXd& prior_prec_alpha2;     // p-vector of precision of alpha2s for each trait
//   const Map<MatrixXd> prior_mean_beta; // b x p matrix of prior means of beta
//   const Map<MatrixXd> prior_prec_beta; // b x p matrix of prior precisions of beta
//   double Y_prec_b0;
//
//   const MatrixXd& randn_alpha1;
//   const std::vector<VectorXd>& randn_alpha2;
//   const MatrixXd& randn_beta;
//   const MatrixXd& randn_e;
//   const VectorXd& rgamma_1;
//
//   MatrixXd& alpha1;
//   std::vector<VectorXd>& alpha2;
//   MatrixXd& beta;
//   VectorXd& Y_prec;
//
//   regression_sampler_workera(
//     const int sampler,
//     const VectorXi& trait_set,
//     const Map<MatrixXd> Y,           // nx1
//     const Map<MatrixXd> W_base,           // nxa
//     const std::vector<R_matrix>& W_list,           // nxa
//     const Map<MatrixXd> X_U,
//     const Map<MatrixXd> V,
//     const MatrixXd RinvSqX,                // nxb
//     const MatrixXd& C,                     // bxb
//     const R_matrix& chol_R,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
//     const MatrixXd& R,
//     const MatrixXd& Rinv,
//     const MatrixXd& RinvU,
//     const MatrixXd& UtRinvU,
//     const Map<MatrixXd> prior_prec_alpha1, // a1 x p matrix of prior precisions for alpha1
//     const VectorXd& prior_prec_alpha2,     // p-vector of precision of alpha2s for each trait
//     const Map<MatrixXd> prior_mean_beta, // b x p matrix of prior means of beta
//     const Map<MatrixXd> prior_prec_beta, // b x p matrix of prior precisions of beta
//     double Y_prec_b0,
//     const MatrixXd& randn_alpha1,
//     const std::vector<VectorXd>& randn_alpha2,
//     const MatrixXd& randn_beta,
//     const MatrixXd& randn_e,
//     const VectorXd& rgamma_1,
//     MatrixXd& alpha1,
//     std::vector<VectorXd>& alpha2,
//     MatrixXd& beta,
//     VectorXd& Y_prec
//   ):
//     sampler(sampler), trait_set(trait_set),
//     Y(Y), W_base(W_base), W_list(W_list), X_U(X_U), V(V), RinvSqX(RinvSqX),
//     C(C), chol_R(chol_R), R(R), Rinv(Rinv), RinvU(RinvU), UtRinvU(UtRinvU),
//     prior_prec_alpha1(prior_prec_alpha1), prior_prec_alpha2(prior_prec_alpha2), prior_mean_beta(prior_mean_beta), prior_prec_beta(prior_prec_beta),
//     Y_prec_b0(Y_prec_b0),
//     randn_alpha1(randn_alpha1), randn_alpha2(randn_alpha2), randn_beta(randn_beta), randn_e(randn_e), rgamma_1(rgamma_1),
//     alpha1(alpha1), alpha2(alpha2), beta(beta), Y_prec(Y_prec)
//   {}
//
//   void operator()(std::size_t begin, std::size_t end) {
//     int n = Y.rows();
//     int a1 = W_base.cols();
//
//     for(std::size_t i = begin; i < end; i++){
//       int j = trait_set[i];
//       MatrixXd W;
//       int a;
//       int a2 = 0;
//       int b;
//       VectorXd prior_prec_alpha;
//       VectorXd randn_alpha;
//       if(W_list.size() == 0) {
//         W = W_base;
//         a = a1;
//         prior_prec_alpha = prior_prec_alpha1.col(j);
//         randn_alpha = randn_alpha1.col(j);
//       } else{
//         Map<MatrixXd> W2 = W_list[j].dense;
//         a2 = W2.cols();
//         a = a1+a2;
//         W = MatrixXd(n,a);
//         W << W_base,W2;
//         prior_prec_alpha = VectorXd(a);
//         prior_prec_alpha.head(a1) = prior_prec_alpha1.col(j);
//         prior_prec_alpha.tail(a2).array() = prior_prec_alpha2[j];
//         randn_alpha = VectorXd(a);
//         randn_alpha.head(a1) = randn_alpha1.col(j);
//         randn_alpha.tail(a2) = randn_alpha2[j];
//       }
//
//       VectorXd samples;
//       if(sampler == 1) {
//         b = RinvSqX.cols();
//         samples = regression_sampler_v1a(Y.col(j), W, RinvSqX, C, prior_prec_alpha, prior_mean_beta.col(j),
//                                         prior_prec_beta.col(j), chol_R, Y_prec[j], randn_alpha,
//                                         randn_beta.col(j), rgamma_1[j],Y_prec_b0);
//       } else if(sampler == 2) {
//         b = X_U.cols();
//         samples = regression_sampler_v2a(Y.col(j), W, X_U, prior_prec_alpha, prior_mean_beta.col(j),
//                                         prior_prec_beta.col(j), chol_R, R, Y_prec[j], randn_alpha,
//                                         randn_beta.col(j), randn_e.col(j),rgamma_1[j],Y_prec_b0);
//       } else if(sampler == 3) {
//         b = V.cols();
//         samples = regression_sampler_v3a(Y.col(j), W, X_U, V, prior_prec_alpha, prior_mean_beta.col(j),
//                                         prior_prec_beta.col(j), chol_R, Rinv, RinvU, UtRinvU, Y_prec[j], randn_alpha,
//                                         randn_beta.col(j), randn_e.col(j),rgamma_1[j],Y_prec_b0);
//
//       } else {
//         stop("sampler not implemented");
//       }
//
//       // extract samples
//       Y_prec[j] = samples[0];
//       if(a1 > 0) alpha1.col(j) = samples.segment(1,a1);
//       if(a2 > 0) alpha2[j] = samples.segment(1+a1,a2);
//       if(b > 0) beta.col(j) = samples.tail(b);
//     }
//   }
// };
//
//
// // [[Rcpp::export]]
// Rcpp::List regression_sampler_parallela(
//     Map<MatrixXd> Y,               // n x p matrix of observations
//     Map<MatrixXd> W_base,          // n x a1 matrix of W covariates common to all p. Can be NULL
//     Rcpp::List W_list_,             // p-list of n x a2 matrices of W covariates unique to each p. Can be NULL
//     Map<MatrixXd> X,               // either X, a n x b matrix, or U, a n x m matrix. If U, then V must be non-NULL
//     SEXP V_,                       // m x b matrix if X is U
//     Rcpp::IntegerVector h2s_index, // p-vector of indices for appropriate V of each trait
//     Rcpp::List chol_V_list_,        // list of cholesky decompositions of V: RtR (each nxn). Can be either dense or sparse
//     VectorXd Y_prec,               // p-vector of Y current precisions
//     double Y_prec_a0,
//     double Y_prec_b0,
//     Map<MatrixXd> prior_prec_alpha1, // a1 x p matrix of prior precisions for alpha1
//     VectorXd prior_prec_alpha2,     // p-vector of precision of alpha2s for each trait
//     Map<MatrixXd> prior_mean_beta, // b x p matrix of prior means of beta
//     Map<MatrixXd> prior_prec_beta, // b x p matrix of prior precisions of beta
//     int grainSize) {
//
//   int n = Y.rows();
//   int p = Y.cols();
//
//   // W_base
//   if(W_base.rows() != n) stop("Wrong dimension of W_base");
//   int a1 = W_base.cols();
//
//   // W_list
//   std::vector<R_matrix> W_list;
//   load_R_matrices_lista(W_list_, W_list);
//   if(W_list.size() > 0) {
//     if(W_list.size() != p) stop("Wrong length of W_list");
//   }
//
//   // X or U and V
//   Map<MatrixXd> U = X;
//   MatrixXd z = MatrixXd::Zero(0,0);
//   Map<MatrixXd> V(z.data(),0,0);
//   int b = X.cols();
//   if(X.rows() != n) stop("Wrong dimension of X");
//   if(Rf_isMatrix(V_)) {
//     // Map<MatrixXd> V__ = as<Map<MatrixXd> >(V_);
//     // new (&v) Map<MatrixXd> (V__,V__.rows(),V__.cols());
//     new (&V) Map<MatrixXd> (as<Map<MatrixXd> >(V_));
//     if(U.cols() != V.rows()) stop("X and V_ have incompatible dimensions");
//     b = V.cols();
//   }
//
//   // chol_V_list
//   std::vector<R_matrix> chol_V_list;
//   load_R_matrices_lista(chol_V_list_, chol_V_list);
//   if(max(h2s_index) > chol_V_list.size()) {
//     stop("max(h2s_index) > length(chol_V_list)");
//   }
//
//   // priors
//   if(Y_prec.size() != p) {
//     stop("Wrong length of Y_prec");
//   }
//   if(prior_prec_alpha1.rows() != a1 || prior_prec_alpha1.cols() != p) stop("Wrong dimensions of prior_prec_alpha1");
//   if(W_list.size() > 0 && prior_prec_alpha2.size() != p) {
//     stop("Wrong length of prior_prec_alpha2");
//   }
//   if(prior_mean_beta.rows() != b || prior_mean_beta.cols() != p) stop("Wrong dimensions of prior_mean_beta");
//   if(prior_prec_beta.rows() != b || prior_prec_beta.cols() != p) stop("Wrong dimensions of prior_prec_beta");
//
//   // generate random numbers
//   MatrixXd randn_alpha1 = rstdnorm_mata(a1,p);
//   std::vector<VectorXd> randn_alpha2;
//   if(W_list.size() > 0){
//     for(int i = 0; i < p; i++){
//       randn_alpha2.push_back(rstdnorm_mata(W_list[i].dense.cols(),1));
//     }
//   }
//   MatrixXd randn_beta = rstdnorm_mata(b,p);
//   MatrixXd randn_e;
//   if(b > n) {
//     randn_e = rstdnorm_mata(n,p);
//   }
//   VectorXd rgamma_1 = as<VectorXd>(rgamma(p,Y_prec_a0 + n/2.0,1.0));
//
//   // Results structures
//   MatrixXd alpha1(a1,p);
//   std::vector<VectorXd> alpha2;
//   alpha2.reserve(W_list.size());
//   int alpha2_size = 0;
//   if(W_list.size() > 0){
//     for(int i = 0; i < W_list.size(); i++){
//       int a2 = W_list[i].dense.cols();
//       alpha2.push_back(VectorXd::Zero(a2));
//       alpha2_size += a2;
//     }
//   }
//   MatrixXd beta(b,p);
//
//   // go through h2s indices and sample columns with same index as a set
//   for(int i = min(h2s_index); i <= max(h2s_index); i++) {
//     int h2_index = i;
//     VectorXi trait_set = as<VectorXi>(whicha(h2s_index == h2_index));  // list of traits with same h2_index
//
//     if(trait_set.size() > 0){
//       // prepare matrices for sampler
//       MatrixXd RinvSqX, C, R, Rinv, RinvU, UtRinvU;
//       R_matrix chol_R = chol_V_list[h2_index - 1];
//       int which_sampler;
//       // Decide whicha sampler to use
//       if(b <= n) {
//         // use regression_sampler_v1a
//         which_sampler = 1;
//         if(chol_R.isDense) {
//           RinvSqX = chol_R.dense.transpose().triangularView<Lower>().solve(X);
//         } else{
//           RinvSqX = chol_R.sparse.transpose().triangularView<Lower>().solve(X);
//         }
//         C = RinvSqX.transpose() * RinvSqX;
//       }
//       else if(V.cols() == 0) {
//         // use regression_sampler_v2a
//         which_sampler = 2;
//         if(chol_R.isDense) {
//           R = chol_R.dense.transpose().triangularView<Lower>() * chol_R.dense;
//         } else{
//           R = chol_R.sparse.transpose().triangularView<Lower>() * chol_R.sparse;
//         }
//       } else {
//         // use regression_sampler_v3a
//         which_sampler = 3;
//         if(chol_R.isDense) {
//           Rinv = chol_R.dense.triangularView<Upper>().solve(chol_R.dense.transpose().triangularView<Lower>().solve(MatrixXd::Identity(n,n)));
//           RinvU = chol_R.dense.triangularView<Upper>().solve(chol_R.dense.transpose().triangularView<Lower>().solve(U));
//         } else{
//           Rinv = chol_R.sparse.triangularView<Upper>().solve(chol_R.sparse.transpose().triangularView<Lower>().solve(MatrixXd::Identity(n,n)));
//           RinvU = chol_R.sparse.triangularView<Upper>().solve(chol_R.sparse.transpose().triangularView<Lower>().solve(U));
//           // RinvU = Rinv * U;
//           // Rcout << i << std::endl;
//           // Rcout << Rinv.diagonal().transpose() << std::endl;
//         }
//         UtRinvU = U.transpose() * RinvU;
//       }
//
//       // Rcout << omp_get_max_threads() << std::endl;
//       // omp_set_num_threads(4);
//       regression_sampler_workera sampler(which_sampler, trait_set,Y,W_base,W_list,X,V,RinvSqX,C,
//                                         chol_R,R,Rinv,RinvU,UtRinvU,
//                                         prior_prec_alpha1,prior_prec_alpha2,prior_mean_beta,prior_prec_beta,Y_prec_b0,
//                                         randn_alpha1,randn_alpha2,randn_beta,randn_e, rgamma_1,
//                                         alpha1,alpha2,beta,Y_prec);
//       RcppParallel::parallelFor(0,trait_set.size(),sampler,grainSize);
//     }
//   }
//
//   // collect alpha2 into a vector
//   VectorXd alpha2_vec(alpha2_size);
//   if(W_list.size() > 0){
//     int index = 0;
//     for(int i = 0; i < W_list.size(); i++){
//       alpha2_vec.segment(index,alpha2[i].size()) = alpha2[i];
//       index += alpha2[i].size();
//     }
//   }
//
//   return(Rcpp::List::create(
//       Named("alpha1") = alpha1,
//       Named("alpha2") = alpha2_vec,
//       Named("beta") = beta,
//       Named("Y_prec") = Y_prec
//   ));
// }
//
//
// // [[Rcpp::export]]
// Rcpp::List regression_sampler_parallelb(
//     Map<MatrixXd> Y,               // n x p matrix of observations
//     Map<MatrixXd> W_base,          // n x a1 matrix of W covariates common to all p. Can be NULL
//     Rcpp::List W_list_,             // p-list of n x a2 matrices of W covariates unique to each p. Can be NULL
//     Map<MatrixXd> X,               // either X, a n x b matrix, or U, a n x m matrix. If U, then V must be non-NULL
//     SEXP V_,                       // m x b matrix if X is U
//     Rcpp::IntegerVector h2s_index, // p-vector of indices for appropriate V of each trait
//     Rcpp::List chol_V_list_,        // list of cholesky decompositions of V: RtR (each nxn). Can be either dense or sparse
//     VectorXd Y_prec,               // p-vector of Y current precisions
//     double Y_prec_a0,
//     double Y_prec_b0,
//     Map<MatrixXd> prior_prec_alpha1, // a1 x p matrix of prior precisions for alpha1
//     VectorXd prior_prec_alpha2,     // p-vector of precision of alpha2s for each trait
//     Map<MatrixXd> prior_mean_beta, // b x p matrix of prior means of beta
//     Map<MatrixXd> prior_prec_beta, // b x p matrix of prior precisions of beta
//     int grainSize) {
//
//   // omp_set_num_threads(4);
//   int n = Y.rows();
//   int p = Y.cols();
//
//   // W_base
//   if(W_base.rows() != n) stop("Wrong dimension of W_base");
//   int a1 = W_base.cols();
//
//   // W_list
//   std::vector<R_matrix> W_list;
//   load_R_matrices_lista(W_list_, W_list);
//   if(W_list.size() > 0) {
//     if(W_list.size() != p) stop("Wrong length of W_list");
//   }
//
//   // X or U and V
//   Map<MatrixXd> U = X;
//   MatrixXd z = MatrixXd::Zero(0,0);
//   Map<MatrixXd> V(z.data(),0,0);
//   int b = X.cols();
//   if(X.rows() != n) stop("Wrong dimension of X");
//   if(Rf_isMatrix(V_)) {
//     // Map<MatrixXd> V__ = as<Map<MatrixXd> >(V_);
//     // new (&v) Map<MatrixXd> (V__,V__.rows(),V__.cols());
//     new (&V) Map<MatrixXd> (as<Map<MatrixXd> >(V_));
//     if(U.cols() != V.rows()) stop("X and V_ have incompatible dimensions");
//     b = V.cols();
//   }
//
//   // chol_V_list
//   std::vector<R_matrix> chol_V_list;
//   load_R_matrices_lista(chol_V_list_, chol_V_list);
//   if(max(h2s_index) > chol_V_list.size()) {
//     stop("max(h2s_index) > length(chol_V_list)");
//   }
//
//   // priors
//   if(Y_prec.size() != p) {
//     stop("Wrong length of Y_prec");
//   }
//   if(prior_prec_alpha1.rows() != a1 || prior_prec_alpha1.cols() != p) stop("Wrong dimensions of prior_prec_alpha1");
//   if(W_list.size() > 0 && prior_prec_alpha2.size() != p) {
//     stop("Wrong length of prior_prec_alpha2");
//   }
//   if(prior_mean_beta.rows() != b || prior_mean_beta.cols() != p) stop("Wrong dimensions of prior_mean_beta");
//   if(prior_prec_beta.rows() != b || prior_prec_beta.cols() != p) stop("Wrong dimensions of prior_prec_beta");
//
//   // generate random numbers
//   MatrixXd randn_alpha1 = rstdnorm_mata(a1,p);
//   std::vector<VectorXd> randn_alpha2;
//   if(W_list.size() > 0){
//     for(int i = 0; i < p; i++){
//       randn_alpha2.push_back(rstdnorm_mata(W_list[i].dense.cols(),1));
//     }
//   }
//   MatrixXd randn_beta = rstdnorm_mata(b,p);
//   MatrixXd randn_e;
//   if(b > n) {
//     randn_e = rstdnorm_mata(n,p);
//   }
//   VectorXd rgamma_1 = as<VectorXd>(rgamma(p,Y_prec_a0 + n/2.0,1.0));
//
//   // Results structures
//   MatrixXd alpha1(a1,p);
//   std::vector<VectorXd> alpha2;
//   alpha2.reserve(W_list.size());
//   int alpha2_size = 0;
//   if(W_list.size() > 0){
//     for(int i = 0; i < W_list.size(); i++){
//       int a2 = W_list[i].dense.cols();
//       alpha2.push_back(VectorXd::Zero(a2));
//       alpha2_size += a2;
//     }
//   }
//   MatrixXd beta(b,p);
//
//   // go through h2s indices and sample columns with same index as a set
//   for(int i = min(h2s_index); i <= max(h2s_index); i++) {
//     int h2_index = i;
//     VectorXi trait_set = as<VectorXi>(whicha(h2s_index == h2_index));  // list of traits with same h2_index
//
//     if(trait_set.size() > 0){
//       // prepare matrices for sampler
//       MatrixXd RinvSqX, C, R, Rinv, RinvU, UtRinvU;
//       R_matrix chol_R = chol_V_list[h2_index - 1];
//       int which_sampler;
//       // Decide whicha sampler to use
//       if(b <= n) {
//         // use regression_sampler_v1a
//         which_sampler = 1;
//         if(chol_R.isDense) {
//           RinvSqX = chol_R.dense.transpose().triangularView<Lower>().solve(X);
//         } else{
//           RinvSqX = chol_R.sparse.transpose().triangularView<Lower>().solve(X);
//         }
//         C = RinvSqX.transpose() * RinvSqX;
//       }
//       else if(V.cols() == 0) {
//         // use regression_sampler_v2a
//         which_sampler = 2;
//         if(chol_R.isDense) {
//           R = chol_R.dense.transpose().triangularView<Lower>() * chol_R.dense;
//         } else{
//           R = chol_R.sparse.transpose().triangularView<Lower>() * chol_R.sparse;
//         }
//       } else {
//         // use regression_sampler_v3a
//         which_sampler = 3;
//         if(chol_R.isDense) {
//           Rinv = chol_R.dense.triangularView<Upper>().solve(chol_R.dense.transpose().triangularView<Lower>().solve(MatrixXd::Identity(n,n)));
//           RinvU = chol_R.dense.triangularView<Upper>().solve(chol_R.dense.transpose().triangularView<Lower>().solve(U));
//         } else{
//           Rinv = chol_R.sparse.triangularView<Upper>().solve(chol_R.sparse.transpose().triangularView<Lower>().solve(MatrixXd::Identity(n,n)));
//           RinvU = chol_R.sparse.triangularView<Upper>().solve(chol_R.sparse.transpose().triangularView<Lower>().solve(U));
//           // RinvU = Rinv * U;
//           // Rcout << i << std::endl;
//           // Rcout << Rinv.diagonal().transpose() << std::endl;
//         }
//         UtRinvU = U.transpose() * RinvU;
//       }
//
//
//       int sampler = which_sampler;
//
//       // Rcout << omp_get_max_threads() << std::endl;
//       #pragma omp parallel for
//       for(int i = 0; i < trait_set.size(); i++){
//         int j = trait_set[i];
//         MatrixXd W;
//         int a;
//         int a2 = 0;
//         int b;
//         VectorXd prior_prec_alpha;
//         VectorXd randn_alpha;
//         if(W_list.size() == 0) {
//           W = W_base;
//           a = a1;
//           prior_prec_alpha = prior_prec_alpha1.col(j);
//           randn_alpha = randn_alpha1.col(j);
//         } else{
//           Map<MatrixXd> W2 = W_list[j].dense;
//           a2 = W2.cols();
//           a = a1+a2;
//           W = MatrixXd(n,a);
//           W << W_base,W2;
//           prior_prec_alpha = VectorXd(a);
//           prior_prec_alpha.head(a1) = prior_prec_alpha1.col(j);
//           prior_prec_alpha.tail(a2).array() = prior_prec_alpha2[j];
//           randn_alpha = VectorXd(a);
//           randn_alpha.head(a1) = randn_alpha1.col(j);
//           randn_alpha.tail(a2) = randn_alpha2[j];
//         }
//
//         VectorXd samples;
//         if(sampler == 1) {
//           b = RinvSqX.cols();
//           samples = regression_sampler_v1a(Y.col(j), W, RinvSqX, C, prior_prec_alpha, prior_mean_beta.col(j),
//                                            prior_prec_beta.col(j), chol_R, Y_prec[j], randn_alpha,
//                                            randn_beta.col(j), rgamma_1[j],Y_prec_b0);
//         } else if(sampler == 2) {
//           b = X.cols();
//           samples = regression_sampler_v2a(Y.col(j), W, X, prior_prec_alpha, prior_mean_beta.col(j),
//                                            prior_prec_beta.col(j), chol_R, R, Y_prec[j], randn_alpha,
//                                            randn_beta.col(j), randn_e.col(j),rgamma_1[j],Y_prec_b0);
//         } else if(sampler == 3) {
//           b = V.cols();
//           samples = regression_sampler_v3a(Y.col(j), W, X, V, prior_prec_alpha, prior_mean_beta.col(j),
//                                            prior_prec_beta.col(j), chol_R, Rinv, RinvU, UtRinvU, Y_prec[j], randn_alpha,
//                                            randn_beta.col(j), randn_e.col(j),rgamma_1[j],Y_prec_b0);
//
//         } else {
//           stop("sampler not implemented");
//         }
//
//         // extract samples
//         Y_prec[j] = samples[0];
//         if(a1 > 0) alpha1.col(j) = samples.segment(1,a1);
//         if(a2 > 0) alpha2[j] = samples.segment(1+a1,a2);
//         if(b > 0) beta.col(j) = samples.tail(b);
//       }
//
//       // regression_sampler_workera sampler(which_sampler, trait_set,Y,W_base,W_list,X,V,RinvSqX,C,
//       //                                    chol_R,R,Rinv,RinvU,UtRinvU,
//       //                                    prior_prec_alpha1,prior_prec_alpha2,prior_mean_beta,prior_prec_beta,Y_prec_b0,
//       //                                    randn_alpha1,randn_alpha2,randn_beta,randn_e, rgamma_1,
//       //                                    alpha1,alpha2,beta,Y_prec);
//       // RcppParallel::parallelFor(0,trait_set.size(),sampler,grainSize);
//     }
//   }
//
//   // collect alpha2 into a vector
//   VectorXd alpha2_vec(alpha2_size);
//   if(W_list.size() > 0){
//     int index = 0;
//     for(int i = 0; i < W_list.size(); i++){
//       alpha2_vec.segment(index,alpha2[i].size()) = alpha2[i];
//       index += alpha2[i].size();
//     }
//   }
//
//   return(Rcpp::List::create(
//       Named("alpha1") = alpha1,
//       Named("alpha2") = alpha2_vec,
//       Named("beta") = beta,
//       Named("Y_prec") = Y_prec
//   ));
// }
//
// //
// //
// //
// // typedef Eigen::SparseMatrix<float> SpMatf;
// //
// //
// // // [[Rcpp::export()]]
// // MatrixXf rstdnorm_matb(int n,int p) {  // returns nxp matrix
// //   VectorXf X_vec(n*p);
// //   for(int i = 0; i < n*p; i++){
// //     X_vec[i] = ziggr.norm();
// //   }
// //   Map<MatrixXf> X_mat(X_vec.data(),n,p);
// //   return(X_mat);
// // }
// //
// // struct R_matrixb {
// //   MatrixXf dense;
// //   SpMatf sparse;
// //   bool isDense;
// //   R_matrixb(MatrixXf dense_, SpMatf sparse_,bool isDense_) : dense(dense_), sparse(sparse_), isDense(isDense_) {}
// // };
// // // void load_R_matrices_list(const Rcpp::List X_list, std::vector<R_matrixb>& X_vector);
// //
// //
// // void load_R_matrices_listb(const Rcpp::List X_list, std::vector<R_matrixb>& X_vector){
// //   // null_matrices
// //   MatrixXf M_null_d = MatrixXf::Zero(0,0);
// //   // MatrixXf M_null_d(null_d.data(),0,0);
// //   SpMatf M_null_s = M_null_d.sparseView();
// //   // SpMatf M_null_s(0,0,0,null_s.outerIndexPtr(),null_s.innerIndexPtr(),null_s.valuePtr());
// //
// //   int p = X_list.size();
// //   X_vector.reserve(p);
// //   for(int i = 0; i < p; i++){
// //     SEXP Xi_ = X_list[i];
// //     if(Rf_isMatrix(Xi_)){
// //       MatrixXf Xi = as<MatrixXf >(Xi_);
// //       R_matrixb Xim(Xi,M_null_s,true);
// //       X_vector.push_back(Xim);
// //     } else{
// //       SpMatf Xi = as<SpMatf>(Xi_);
// //       R_matrixb Xim(M_null_d,Xi,false);
// //       X_vector.push_back(Xim);
// //     }
// //   }
// // }
// //
// // MatrixXf get_RinvSqXb(const R_matrixb& chol_R, MatrixXf X){
// //   MatrixXf RinvSqX;
// //   if(chol_R.isDense) {
// //     RinvSqX = chol_R.dense.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
// //   } else{
// //     RinvSqX = chol_R.sparse.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
// //   }
// //   return(RinvSqX);
// // }
// //
// //
// //
// // VectorXf regression_sampler_v1b(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b < n
// //     const Ref<const VectorXf>& y,           // nx1
// //     const MatrixXf& W,           // nxa
// //     const MatrixXf& RinvSqX,                // nxb
// //     const MatrixXf& C,                     // bxb
// //     const VectorXf& prior_prec_alpha, // ax 1
// //     const Ref<const VectorXf>& prior_mean_beta,  // bx1
// //     const Ref<const VectorXf>& prior_prec_beta,  // bx1
// //     const R_matrixb& chol_R,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
// //     double Y_prec,                    // double
// //     const VectorXf& randn_alpha,
// //     const Ref<const VectorXf>& randn_beta,
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
// //   MatrixXf C_beta = C;
// //   C_beta.diagonal() += prior_prec_beta;
// //   LLT<MatrixXf> A_beta_llt;
// //   A_beta_llt.compute(C_beta);
// //   MatrixXf chol_A_beta = A_beta_llt.matrixU();
// //   // Y_prec * chol_A_beta^\T * chol_A_beta = A_beta
// //
// //   // Step 1
// //   VectorXf alpha(a);
// //   VectorXf y_tilde = y;
// //   if(a > 0) {
// //     // Sample alpha
// //     // Calculate A_alpha = Y_prec*W^T*Sigma_beta^{-1}*W + D_alpha^{-1}
// //     // We don't need to actually calculate Sigma_beta^{-1} directly.
// //     MatrixXf RinvSqW = get_RinvSqXb(chol_R,W);  // n*n*a -> n x a
// //     MatrixXf WtRinvX = RinvSqW.transpose() * RinvSqX; // a*n*b -> a*b
// //     MatrixXf invSqAbXtRinvW = chol_A_beta.transpose().triangularView<Lower>().solve(WtRinvX.transpose()); // b*b*a -> b x a
// //
// //     MatrixXf A_alpha = Y_prec * (RinvSqW.transpose() * RinvSqW - invSqAbXtRinvW.transpose() * invSqAbXtRinvW);
// //     A_alpha.diagonal() += prior_prec_alpha;
// //
// //     VectorXf Rinvsqy = get_RinvSqXb(chol_R,y); // n*n*q -> n x 1;
// //     VectorXf XtRinvy = RinvSqX.transpose() * Rinvsqy; // b*n*1 >- b x 1
// //     VectorXf invSqAbXtRinvy = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy); // b*b*1 -> b*1
// //
// //     VectorXf WtSbinvy = RinvSqW.transpose() * Rinvsqy - invSqAbXtRinvW.transpose() * invSqAbXtRinvy;
// //
// //     LLT<MatrixXf> A_alpha_llt;
// //     A_alpha_llt.compute(A_alpha);
// //     MatrixXf chol_A_alpha = A_alpha_llt.matrixU();
// //
// //     alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(WtSbinvy) * Y_prec + randn_alpha;
// //     alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
// //
// //     y_tilde = y - W * alpha;
// //   }
// //
// //   // Step 2 - sample Y_prec
// //   // We don't need to actually calculate Sigma_beta^{-1} directly.
// //   VectorXf RinvSqy = get_RinvSqXb(chol_R,y_tilde);
// //   VectorXf XtRinvy = RinvSqX.transpose() * RinvSqy;
// //   VectorXf prod1 = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy);
// //   double score = Y_prec_b0 + (RinvSqy.dot(RinvSqy) - prod1.dot(prod1))/2;
// //   Y_prec = rgamma_1/score;
// //
// //   // Step 3 - sample beta
// //   VectorXf XtRinvy_std_mu = XtRinvy*Y_prec + prior_prec_beta.asDiagonal()*prior_mean_beta;
// //   VectorXf beta = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy_std_mu) / sqrt(Y_prec) + randn_beta;
// //   beta = chol_A_beta.triangularView<Upper>().solve(beta) / sqrt(Y_prec);
// //
// //   VectorXf result(1+a+b);
// //   result << Y_prec,alpha,beta;
// //
// //   return(result);
// // }
// //
// //
// // VectorXf regression_sampler_v2b(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b > n
// //     const Ref<const VectorXf>& y,           // nx1
// //     const MatrixXf& W,           // nxa
// //     const MatrixXf& X,           // nxm or nxb
// //     const VectorXf& prior_prec_alpha, // ax 1
// //     const Ref<const VectorXf>& prior_mean_beta,  // bx1
// //     const Ref<const VectorXf>& prior_prec_beta,  // bx1
// //     const R_matrixb& chol_R,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
// //     const MatrixXf& R,
// //     double Y_prec,                    // double
// //     const VectorXf& randn_alpha,
// //     const Ref<const VectorXf>& randn_beta,
// //     const Ref<const VectorXf>& randn_e,
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
// //   MatrixXf DXt = prior_prec_beta.cwiseInverse().asDiagonal() * X.transpose();
// //   MatrixXf Sigma_beta = X * DXt + R;
// //   LDLT<MatrixXf> Sigma_beta_ldlt;
// //   Sigma_beta_ldlt.compute(Sigma_beta);
// //
// //   // Step 1
// //   VectorXf alpha(a);
// //   VectorXf y_tilde = y;
// //   if(a > 0) {
// //     // Sample alpha
// //     MatrixXf SbinvW = Sigma_beta_ldlt.solve(W);
// //     MatrixXf A_alpha = Y_prec * SbinvW.transpose() * W;
// //     A_alpha.diagonal() += prior_prec_alpha;
// //
// //     LLT<MatrixXf> A_alpha_llt;
// //     A_alpha_llt.compute(A_alpha);
// //     MatrixXf chol_A_alpha = A_alpha_llt.matrixU();
// //
// //     VectorXf WtSbinvy = SbinvW.transpose() * y;
// //     alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(WtSbinvy) * Y_prec + randn_alpha;
// //     alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
// //     y_tilde = y - W * alpha;
// //   }
// //
// //   // Step 2 - sample Y_prec
// //   VectorXf e2 = y_tilde.transpose() * Sigma_beta_ldlt.solve(y_tilde);
// //   double score = Y_prec_b0 + e2[0]/2;
// //   Y_prec = rgamma_1/score;
// //
// //   // Step 3 - sample beta
// //   // what about prior mean?
// //   VectorXf u = randn_beta.array() / (prior_prec_beta * Y_prec).cwiseSqrt().array() + prior_mean_beta.array();
// //   VectorXf v = sqrt(Y_prec) * X * u;
// //   if(chol_R.isDense) {
// //     v += chol_R.dense.transpose().triangularView<Lower>() * randn_e;
// //   } else{
// //     v += chol_R.sparse.transpose().triangularView<Lower>() * randn_e;
// //   }
// //   VectorXf w = Sigma_beta_ldlt.solve(y_tilde * sqrt(Y_prec) - v);
// //   VectorXf beta = u + DXt * w / sqrt(Y_prec);
// //
// //   VectorXf result(1+a+b);
// //   result << Y_prec,alpha,beta;
// //
// //   return(result);
// // }
// //
// //
// // VectorXf regression_sampler_v3b(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b > n > m
// //     const Ref<const VectorXf>& y,           // nx1
// //     const MatrixXf& W,           // nxa
// //     const MatrixXf& U,           // nxm or nxb
// //     const MatrixXf& V,           // mxb
// //     const VectorXf& prior_prec_alpha, // ax 1
// //     const Ref<const VectorXf>& prior_mean_beta,  // bx1
// //     const Ref<const VectorXf>& prior_prec_beta,  // bx1
// //     const R_matrixb& chol_R,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
// //     const MatrixXf& Rinv,
// //     const MatrixXf& RinvU,
// //     const MatrixXf& UtRinvU,
// //     double Y_prec,                    // double
// //     const VectorXf& randn_alpha,
// //     const Ref<const VectorXf>& randn_beta,
// //     const Ref<const VectorXf>& randn_e,
// //     const double rgamma_1,
// //     const double Y_prec_b0
// // ){
// //   int n = y.size();
// //   int a = W.cols();
// //   if(V.rows() != U.cols()) stop("Wrong dimensions of V");
// //   MatrixXf X = U*V;
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
// //   if(randn_e.size() != n) stop("Wrong length of randn_e");
// //
// //   // Calculate inverse of Sigma_beta
// //   // MatrixXf Sigma_beta_inv;
// //   // Using Ainv - Ainv * U * (I + BVAinvU)inv * BVAinv in case B = VDVt is singular
// //   MatrixXf DVt = prior_prec_beta.cwiseInverse().asDiagonal() * V.transpose();
// //   MatrixXf VDVt = V * DVt;
// //   if(RinvU.rows() != n) stop("Wrong dimensions of RinvU");
// //   MatrixXf inner = VDVt * UtRinvU;
// //   inner.diagonal().array() += 1.0;
// //   LDLT<MatrixXf> inner_ldlt;
// //   inner_ldlt.compute(inner);
// //   // Sigma_beta_inv = Rinv - RinvU * inner.ldlt().solve(VDVt * RinvU.transpose());  // Don't actually calculate this. Stay in mxm space
// //
// //   // Step 1
// //   VectorXf alpha(a);
// //   VectorXf y_tilde = y;
// //   if(a > 0) {
// //     // Sample alpha
// //     // MatrixXf SbinvW = Sigma_beta_inv * W;
// //     // MatrixXf A_alpha = Y_prec * SbinvW.transpose() * W;
// //     MatrixXf RinvW = Rinv * W;
// //     MatrixXf UtRinvW = U.transpose() * RinvW;
// //     MatrixXf A_alpha = Y_prec * (W.transpose() * RinvW - UtRinvW.transpose() * inner_ldlt.solve(VDVt * UtRinvW));
// //     A_alpha.diagonal() += prior_prec_alpha;
// //
// //     // VectorXf WtSbinvy = SbinvW.transpose() * y;
// //     VectorXf UtRinvy = RinvU.transpose() * y;
// //     VectorXf WtSbinvy = RinvW.transpose() * y - UtRinvW.transpose() * inner_ldlt.solve(VDVt * UtRinvy);
// //
// //     LLT<MatrixXf> A_alpha_llt;
// //     A_alpha_llt.compute(A_alpha);
// //     MatrixXf chol_A_alpha = A_alpha_llt.matrixU();
// //
// //     alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(WtSbinvy) * Y_prec + randn_alpha;
// //     alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
// //
// //     y_tilde = y - W * alpha;
// //   }
// //
// //   // Step 2 - sample Y_prec
// //   // VectorXf e2 = y_tilde.transpose() * Sigma_beta_inv * y_tilde;
// //   VectorXf Rinv_y = Rinv * y_tilde;
// //   VectorXf UtRinvy = RinvU.transpose() * y_tilde;
// //   VectorXf e2 = y_tilde.transpose() * Rinv_y - UtRinvy.transpose() * inner_ldlt.solve(VDVt * UtRinvy);
// //
// //   double score = Y_prec_b0 + e2[0]/2;
// //   Y_prec = rgamma_1/score;
// //
// //   // Step 3 - sample beta
// //   // what about prior mean?
// //   VectorXf u = randn_beta.array() / (prior_prec_beta * Y_prec).cwiseSqrt().array() + prior_mean_beta.array();
// //   VectorXf v = sqrt(Y_prec) * X * u;
// //   if(chol_R.isDense) {
// //     v += chol_R.dense.transpose().triangularView<Lower>() * randn_e;
// //   } else{
// //     v += chol_R.sparse.transpose().triangularView<Lower>() * randn_e;
// //   }
// //   // VectorXf w = Sigma_beta_inv * (y_tilde * sqrt(Y_prec) - v);
// //   VectorXf e = y_tilde * sqrt(Y_prec) - v;
// //   VectorXf UtRinve = RinvU.transpose() * e;
// //   VectorXf w = Rinv * e - RinvU * inner_ldlt.solve(VDVt * UtRinve);
// //
// //   VectorXf beta = u + DVt * (U.transpose() * w) / sqrt(Y_prec); //b*b*1 + b*n*1
// //
// //   VectorXf result(1+a+b);
// //   result << Y_prec,alpha,beta;
// //
// //   return(result);
// // }
// //
// //
// //
// // struct regression_sampler_workerb : public RcppParallel::Worker {
// //   const int sampler;  // whicha sampler to use?
// //   const VectorXi& trait_set;
// //   const MatrixXf Y;           // nx1
// //   const MatrixXf W_base;           // nxa
// //   const std::vector<R_matrixb>& W_list;           // nxa
// //   const MatrixXf X_U;   // could be X or U.
// //   const MatrixXf V;
// //   const MatrixXf RinvSqX;                // nxb
// //   const MatrixXf& C;                     // bxb
// //   const R_matrixb& chol_R;                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
// //   const MatrixXf& R;
// //   const MatrixXf& Rinv;
// //   const MatrixXf& RinvU;
// //   const MatrixXf& UtRinvU;
// //   const MatrixXf prior_prec_alpha1; // a1 x p matrix of prior precisions for alpha1
// //   const VectorXf& prior_prec_alpha2;     // p-vector of precision of alpha2s for each trait
// //   const MatrixXf prior_mean_beta; // b x p matrix of prior means of beta
// //   const MatrixXf prior_prec_beta; // b x p matrix of prior precisions of beta
// //   double Y_prec_b0;
// //
// //   const MatrixXf& randn_alpha1;
// //   const std::vector<VectorXf>& randn_alpha2;
// //   const MatrixXf& randn_beta;
// //   const MatrixXf& randn_e;
// //   const VectorXf& rgamma_1;
// //
// //   MatrixXf& alpha1;
// //   std::vector<VectorXf>& alpha2;
// //   MatrixXf& beta;
// //   VectorXf& Y_prec;
// //
// //   regression_sampler_workerb(
// //     const int sampler,
// //     const VectorXi& trait_set,
// //     const MatrixXf Y,           // nx1
// //     const MatrixXf W_base,           // nxa
// //     const std::vector<R_matrixb>& W_list,           // nxa
// //     const MatrixXf X_U,
// //     const MatrixXf V,
// //     const MatrixXf RinvSqX,                // nxb
// //     const MatrixXf& C,                     // bxb
// //     const R_matrixb& chol_R,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
// //     const MatrixXf& R,
// //     const MatrixXf& Rinv,
// //     const MatrixXf& RinvU,
// //     const MatrixXf& UtRinvU,
// //     const MatrixXf prior_prec_alpha1, // a1 x p matrix of prior precisions for alpha1
// //     const VectorXf& prior_prec_alpha2,     // p-vector of precision of alpha2s for each trait
// //     const MatrixXf prior_mean_beta, // b x p matrix of prior means of beta
// //     const MatrixXf prior_prec_beta, // b x p matrix of prior precisions of beta
// //     double Y_prec_b0,
// //     const MatrixXf& randn_alpha1,
// //     const std::vector<VectorXf>& randn_alpha2,
// //     const MatrixXf& randn_beta,
// //     const MatrixXf& randn_e,
// //     const VectorXf& rgamma_1,
// //     MatrixXf& alpha1,
// //     std::vector<VectorXf>& alpha2,
// //     MatrixXf& beta,
// //     VectorXf& Y_prec
// //   ):
// //     sampler(sampler), trait_set(trait_set),
// //     Y(Y), W_base(W_base), W_list(W_list), X_U(X_U), V(V), RinvSqX(RinvSqX),
// //     C(C), chol_R(chol_R), R(R), Rinv(Rinv), RinvU(RinvU), UtRinvU(UtRinvU),
// //     prior_prec_alpha1(prior_prec_alpha1), prior_prec_alpha2(prior_prec_alpha2), prior_mean_beta(prior_mean_beta), prior_prec_beta(prior_prec_beta),
// //     Y_prec_b0(Y_prec_b0),
// //     randn_alpha1(randn_alpha1), randn_alpha2(randn_alpha2), randn_beta(randn_beta), randn_e(randn_e), rgamma_1(rgamma_1),
// //     alpha1(alpha1), alpha2(alpha2), beta(beta), Y_prec(Y_prec)
// //   {}
// //
// //   void operator()(std::size_t begin, std::size_t end) {
// //     int n = Y.rows();
// //     int a1 = W_base.cols();
// //
// //     for(std::size_t i = begin; i < end; i++){
// //       int j = trait_set[i];
// //       MatrixXf W;
// //       int a;
// //       int a2 = 0;
// //       int b;
// //       VectorXf prior_prec_alpha;
// //       VectorXf randn_alpha;
// //       if(W_list.size() == 0) {
// //         W = W_base;
// //         a = a1;
// //         prior_prec_alpha = prior_prec_alpha1.col(j);
// //         randn_alpha = randn_alpha1.col(j);
// //       } else{
// //         MatrixXf W2 = W_list[j].dense;
// //         a2 = W2.cols();
// //         a = a1+a2;
// //         W = MatrixXf(n,a);
// //         W << W_base,W2;
// //         prior_prec_alpha = VectorXf(a);
// //         prior_prec_alpha.head(a1) = prior_prec_alpha1.col(j);
// //         prior_prec_alpha.tail(a2).array() = prior_prec_alpha2[j];
// //         randn_alpha = VectorXf(a);
// //         randn_alpha.head(a1) = randn_alpha1.col(j);
// //         randn_alpha.tail(a2) = randn_alpha2[j];
// //       }
// //
// //       VectorXf samples;
// //       if(sampler == 1) {
// //         b = RinvSqX.cols();
// //         samples = regression_sampler_v1b(Y.col(j), W, RinvSqX, C, prior_prec_alpha, prior_mean_beta.col(j),
// //                                          prior_prec_beta.col(j), chol_R, Y_prec[j], randn_alpha,
// //                                          randn_beta.col(j), rgamma_1[j],Y_prec_b0);
// //       } else if(sampler == 2) {
// //         b = X_U.cols();
// //         samples = regression_sampler_v2b(Y.col(j), W, X_U, prior_prec_alpha, prior_mean_beta.col(j),
// //                                          prior_prec_beta.col(j), chol_R, R, Y_prec[j], randn_alpha,
// //                                          randn_beta.col(j), randn_e.col(j),rgamma_1[j],Y_prec_b0);
// //       } else if(sampler == 3) {
// //         b = V.cols();
// //         samples = regression_sampler_v3b(Y.col(j), W, X_U, V, prior_prec_alpha, prior_mean_beta.col(j),
// //                                          prior_prec_beta.col(j), chol_R, Rinv, RinvU, UtRinvU, Y_prec[j], randn_alpha,
// //                                          randn_beta.col(j), randn_e.col(j),rgamma_1[j],Y_prec_b0);
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
// // // [[Rcpp::export]]
// // Rcpp::List regression_sampler_parallelb(
// //     MatrixXf Y,               // n x p matrix of observations
// //     MatrixXf W_base,          // n x a1 matrix of W covariates common to all p. Can be NULL
// //     Rcpp::List W_list_,             // p-list of n x a2 matrices of W covariates unique to each p. Can be NULL
// //     MatrixXf X,               // either X, a n x b matrix, or U, a n x m matrix. If U, then V must be non-NULL
// //     SEXP V_,                       // m x b matrix if X is U
// //     Rcpp::IntegerVector h2s_index, // p-vector of indices for appropriate V of each trait
// //     Rcpp::List chol_V_list_,        // list of cholesky decompositions of V: RtR (each nxn). Can be either dense or sparse
// //     VectorXf Y_prec,               // p-vector of Y current precisions
// //     double Y_prec_a0,
// //     double Y_prec_b0,
// //     MatrixXf prior_prec_alpha1, // a1 x p matrix of prior precisions for alpha1
// //     VectorXf prior_prec_alpha2,     // p-vector of precision of alpha2s for each trait
// //     MatrixXf prior_mean_beta, // b x p matrix of prior means of beta
// //     MatrixXf prior_prec_beta, // b x p matrix of prior precisions of beta
// //     int grainSize) {
// //
// //   int n = Y.rows();
// //   int p = Y.cols();
// //
// //   // W_base
// //   if(W_base.rows() != n) stop("Wrong dimension of W_base");
// //   int a1 = W_base.cols();
// //
// //   // W_list
// //   std::vector<R_matrixb> W_list;
// //   load_R_matrices_listb(W_list_, W_list);
// //   if(W_list.size() > 0) {
// //     if(W_list.size() != p) stop("Wrong length of W_list");
// //   }
// //
// //   // X or U and V
// //   MatrixXf U = X;
// //   // MatrixXf z = MatrixXf::Zero(0,0);
// //   MatrixXf V = MatrixXf::Zero(0,0);
// //   // MatrixXf V(z.data(),0,0);
// //   int b = X.cols();
// //   if(X.rows() != n) stop("Wrong dimension of X");
// //   if(Rf_isMatrix(V_)) {
// //     // MatrixXf V__ = as<MatrixXf >(V_);
// //     // new (&v) MatrixXf (V__,V__.rows(),V__.cols());
// //     new (&V) MatrixXf (as<MatrixXf >(V_));
// //     if(U.cols() != V.rows()) stop("X and V_ have incompatible dimensions");
// //     b = V.cols();
// //   }
// //
// //   // chol_V_list
// //   std::vector<R_matrixb> chol_V_list;
// //   load_R_matrices_listb(chol_V_list_, chol_V_list);
// //   if(max(h2s_index) > chol_V_list.size()) {
// //     stop("max(h2s_index) > length(chol_V_list)");
// //   }
// //
// //   // priors
// //   if(Y_prec.size() != p) {
// //     stop("Wrong length of Y_prec");
// //   }
// //   if(prior_prec_alpha1.rows() != a1 || prior_prec_alpha1.cols() != p) stop("Wrong dimensions of prior_prec_alpha1");
// //   if(W_list.size() > 0 && prior_prec_alpha2.size() != p) {
// //     stop("Wrong length of prior_prec_alpha2");
// //   }
// //   if(prior_mean_beta.rows() != b || prior_mean_beta.cols() != p) stop("Wrong dimensions of prior_mean_beta");
// //   if(prior_prec_beta.rows() != b || prior_prec_beta.cols() != p) stop("Wrong dimensions of prior_prec_beta");
// //
// //   // generate random numbers
// //   MatrixXf randn_alpha1 = rstdnorm_matb(a1,p);
// //   std::vector<VectorXf> randn_alpha2;
// //   if(W_list.size() > 0){
// //     for(int i = 0; i < p; i++){
// //       randn_alpha2.push_back(rstdnorm_matb(W_list[i].dense.cols(),1));
// //     }
// //   }
// //   MatrixXf randn_beta = rstdnorm_matb(b,p);
// //   MatrixXf randn_e;
// //   if(b > n) {
// //     randn_e = rstdnorm_matb(n,p);
// //   }
// //   VectorXf rgamma_1 = as<VectorXf>(rgamma(p,Y_prec_a0 + n/2.0,1.0));
// //
// //   // Results structures
// //   MatrixXf alpha1(a1,p);
// //   std::vector<VectorXf> alpha2;
// //   alpha2.reserve(W_list.size());
// //   int alpha2_size = 0;
// //   if(W_list.size() > 0){
// //     for(int i = 0; i < W_list.size(); i++){
// //       int a2 = W_list[i].dense.cols();
// //       alpha2.push_back(VectorXf::Zero(a2));
// //       alpha2_size += a2;
// //     }
// //   }
// //   MatrixXf beta(b,p);
// //
// //   // go through h2s indices and sample columns with same index as a set
// //   for(int i = min(h2s_index); i <= max(h2s_index); i++) {
// //     int h2_index = i;
// //     VectorXi trait_set = as<VectorXi>(whicha(h2s_index == h2_index));  // list of traits with same h2_index
// //
// //     if(trait_set.size() > 0){
// //       // prepare matrices for sampler
// //       MatrixXf RinvSqX, C, R, Rinv, RinvU, UtRinvU;
// //       R_matrixb chol_R = chol_V_list[h2_index - 1];
// //       int which_sampler;
// //       // Decide whicha sampler to use
// //       if(b <= n) {
// //         // use regression_sampler_v1b
// //         which_sampler = 1;
// //         if(chol_R.isDense) {
// //           RinvSqX = chol_R.dense.transpose().triangularView<Lower>().solve(X);
// //         } else{
// //           RinvSqX = chol_R.sparse.transpose().triangularView<Lower>().solve(X);
// //         }
// //         C = RinvSqX.transpose() * RinvSqX;
// //       }
// //       else if(V.cols() == 0) {
// //         // use regression_sampler_v2b
// //         which_sampler = 2;
// //         if(chol_R.isDense) {
// //           R = chol_R.dense.transpose().triangularView<Lower>() * chol_R.dense;
// //         } else{
// //           R = chol_R.sparse.transpose().triangularView<Lower>() * chol_R.sparse;
// //         }
// //       } else {
// //         // use regression_sampler_v3b
// //         which_sampler = 3;
// //         if(chol_R.isDense) {
// //           Rinv = chol_R.dense.triangularView<Upper>().solve(chol_R.dense.transpose().triangularView<Lower>().solve(MatrixXf::Identity(n,n)));
// //           RinvU = chol_R.dense.triangularView<Upper>().solve(chol_R.dense.transpose().triangularView<Lower>().solve(U));
// //         } else{
// //           Rinv = chol_R.sparse.triangularView<Upper>().solve(chol_R.sparse.transpose().triangularView<Lower>().solve(MatrixXf::Identity(n,n)));
// //           RinvU = chol_R.sparse.triangularView<Upper>().solve(chol_R.sparse.transpose().triangularView<Lower>().solve(U));
// //           // RinvU = Rinv * U;
// //           // Rcout << i << std::endl;
// //           // Rcout << Rinv.diagonal().transpose() << std::endl;
// //         }
// //         UtRinvU = U.transpose() * RinvU;
// //       }
// //       regression_sampler_workerb sampler(which_sampler, trait_set,Y,W_base,W_list,X,V,RinvSqX,C,
// //                                          chol_R,R,Rinv,RinvU,UtRinvU,
// //                                          prior_prec_alpha1,prior_prec_alpha2,prior_mean_beta,prior_prec_beta,Y_prec_b0,
// //                                          randn_alpha1,randn_alpha2,randn_beta,randn_e, rgamma_1,
// //                                          alpha1,alpha2,beta,Y_prec);
// //       RcppParallel::parallelFor(0,trait_set.size(),sampler,grainSize);
// //     }
// //   }
// //
// //   // collect alpha2 into a vector
// //   VectorXf alpha2_vec(alpha2_size);
// //   if(W_list.size() > 0){
// //     int index = 0;
// //     for(int i = 0; i < W_list.size(); i++){
// //       alpha2_vec.segment(index,alpha2[i].size()) = alpha2[i];
// //       index += alpha2[i].size();
// //     }
// //   }
// //
// //   return(Rcpp::List::create(
// //       Named("alpha1") = alpha1,
// //       Named("alpha2") = alpha2_vec,
// //       Named("beta") = beta,
// //       Named("Y_prec") = Y_prec
// //   ));
// // }
// //
// //
