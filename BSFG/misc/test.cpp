// #include <RcppArmadillo.h>
// // [[Rcpp::depends("RcppArmadillo")]]
//
// #ifdef _OPENMP
// #include <omp.h>
// #endif
// // [[Rcpp::plugins(openmp)]]
// // [[Rcpp::depends(RcppProgress)]]
// #include <progress.hpp>
// #include <progress_bar.hpp>
//
// // [[Rcpp::export]]
// double long_computation_omp2(int nb, int threads=1) {
// #ifdef _OPENMP
//   if ( threads > 0 )
//     omp_set_num_threads( threads );
//   REprintf("Number of threads=%i\\n", omp_get_max_threads());
// #endif
//   Progress p(nb, true);
//   double sum = 0;
// #pragma omp parallel for schedule(dynamic)
//   for (int i = 0; i < nb; ++i) {
//     double thread_sum = 0;
//     if ( ! Progress::check_abort() ) {
//       p.increment(); // update progress
//       for (int j = 0; j < nb; ++j) {
//         thread_sum += R::dlnorm(i+j, 0.0, 1.0, 0);
//       }
//     }
//     sum += thread_sum;
//   }
//
//   return sum + nb;
// }

//
//
// // [[Rcpp::plugins(openmp)]]
// #include <omp.h>
#include <math.h>
#include <iostream>

// [[Rcpp::depends(RcppEigen,RcppParallel)]]
#include <RcppEigen.h>
#include <RcppParallel.h>


using Eigen::Map;               	      // 'Eigen::Maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::VectorXi;                  // variable size vector, int
using Eigen::ArrayXXd;                  // variable size array, double precision
using Eigen::ArrayXd;                  // variable size array, double precision
using Eigen::Upper;
using Eigen::Lower;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Map<SpMat> MSpMat;




// [[Rcpp::export()]]
MatrixXd rstdnorm_mat(int n,int p) {  // returns nxp matrix
  VectorXd X_vec = Rcpp::as<VectorXd>(Rcpp::rnorm(n*p));
  Map<MatrixXd> X_mat(X_vec.data(),n,p);
  return(X_mat);
}

VectorXd sample_coefs_uncorrelated(  // returns bx1
    VectorXd y,             // nx1
    MatrixXd X,             // nxb
    VectorXd prior_mean,    // bx1
    VectorXd prior_prec,    // bx1
    ArrayXd  resid_prec,    // nx1
    VectorXd randn_theta,   // bx1
    VectorXd randn_e,       // nx1
    int b,
    int n
) {
  VectorXd R_sq_diag = resid_prec.sqrt().inverse();
  VectorXd theta_star = randn_theta.array()/prior_prec.array().sqrt();
  theta_star += prior_mean;
  VectorXd e_star = randn_e.array() * R_sq_diag.array();
  MatrixXd QtW_theta_star = X * theta_star;
  VectorXd eta_resid = y - QtW_theta_star - e_star;
  MatrixXd RinvSqQtW = R_sq_diag.cwiseInverse().asDiagonal() * X;
  VectorXd eta_std = eta_resid.array()/R_sq_diag.array();
  VectorXd WtURinvy = RinvSqQtW.transpose() * eta_std;

  VectorXd theta_tilda;
  if(b < n) {
    MatrixXd C = RinvSqQtW.transpose() * RinvSqQtW;
    C.diagonal() += prior_prec;
    theta_tilda = C.llt().solve(WtURinvy);
  } else{
    MatrixXd VAi = X * prior_prec.cwiseInverse().asDiagonal();
    MatrixXd inner = VAi*X.transpose();
    inner.diagonal() += R_sq_diag.cwiseProduct(R_sq_diag);
    VectorXd VAiWtURinvy = VAi * WtURinvy;
    VectorXd outerWtURinvy = VAi.transpose() * inner.ldlt().solve(VAiWtURinvy);
    theta_tilda = WtURinvy.array() / prior_prec.array();
    theta_tilda -= outerWtURinvy;
  }

  VectorXd coefs = theta_tilda + theta_star;
  return coefs;
}

// RcppParallel struct for sampling coefficients from a set of parallel regression problems:
// Y = X %*% B + E, where
// b[i,j] ~ N(prior_mean[i,j],1/prior_prec[i,j])
// E[i,j] ~ N(0,1/resid_prec[i,j])
// Where Y is complete, so everythign has the same dimensions
struct sample_coefs_uncorrelated_worker : public RcppParallel::Worker {
  MatrixXd Y, X;
  MatrixXd resid_prec, prior_mean, prior_prec, randn_theta, randn_e;
  int b,n;
  MatrixXd &coefs;

  sample_coefs_uncorrelated_worker( MatrixXd Y,           // nxp
                                    MatrixXd X,           // nxb
                                    MatrixXd resid_prec,  // nxp
                                    MatrixXd prior_mean,  // bxp
                                    MatrixXd prior_prec,  // bxp
                                    MatrixXd randn_theta, // bxp
                                    MatrixXd randn_e,     // nxp
                                    int b, int n,
                                    MatrixXd &coefs       // bxp
  ) :
    Y(Y), X(X), resid_prec(resid_prec), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
    b(b), n(n),
    coefs(coefs) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++){
      coefs.col(j) = sample_coefs_uncorrelated(Y.col(j), X, prior_mean.col(j), prior_prec.col(j), resid_prec.col(j), randn_theta.col(j),randn_e.col(j),b,n);
    }
  }
};


// -------------------------------------------------- //
// -- Versions of the independent regression model --- //
// -------------------------------------------------- //

// basic version - all regression have same dimensions
// [[Rcpp::export()]]
MatrixXd sample_coefs_parallel_sparse_c_Eigen(  // returns bxp matrix
    Map<MatrixXd> Y,            // nxp
    Map<MatrixXd> X,            // nxb
    Map<MatrixXd> resid_prec,   // nxp
    Map<MatrixXd> prior_mean,   // bxp
    Map<MatrixXd> prior_prec,   // bxp
    int grainSize){

  // Sample regression coefficients
  // columns of matrices are independent
  // each column conditional posterior is a MVN due to conjugacy

  int p = Y.cols();
  int b = X.cols();
  int n = X.rows();

  MatrixXd randn_theta = rstdnorm_mat(b,p);
  MatrixXd randn_e = rstdnorm_mat(n,p);

  MatrixXd coefs(b,p);

  sample_coefs_uncorrelated_worker sampler(Y, X, resid_prec, prior_mean, prior_prec, randn_theta,randn_e,b, n, coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);

  return(coefs);
}

// // [[Rcpp::export()]]
// MatrixXd sample_coefs_parallel_sparse_c_Eigen_openmp(  // returns bxp matrix
//     Map<MatrixXd> Y,            // nxp
//     Map<MatrixXd> X,            // nxb
//     Map<MatrixXd> resid_prec,   // nxp
//     Map<MatrixXd> prior_mean,   // bxp
//     Map<MatrixXd> prior_prec,   // bxp
//     int grainSize){
//
//   // Sample regression coefficients
//   // columns of matrices are independent
//   // each column conditional posterior is a MVN due to conjugacy
//   int p = Y.cols();
//   int b = X.cols();
//   int n = X.rows();
//
//   MatrixXd randn_theta = rstdnorm_mat(b,p);
//   MatrixXd randn_e = rstdnorm_mat(n,p);
//
//   MatrixXd coefs(b,p);
//
//
// #pragma omp parallel for num_threads(grainSize)
//   for(int j = 0; j < p; j++){
//     VectorXd a = sample_coefs_uncorrelated(Y.col(j), X, prior_mean.col(j), prior_prec.col(j), resid_prec.col(j), randn_theta.col(j),randn_e.col(j),b,n);
//   }
//
//   return(coefs);
// }
// [[Rcpp::export()]]
MatrixXd sample_coefs_sparse_c_Eigen(  // returns bxp matrix
    Map<MatrixXd> Y,            // nxp
    Map<MatrixXd> X,            // nxb
    Map<MatrixXd> resid_prec,   // nxp
    Map<MatrixXd> prior_mean,   // bxp
    Map<MatrixXd> prior_prec,   // bxp
    int grainSize){

  // Sample regression coefficients
  // columns of matrices are independent
  // each column conditional posterior is a MVN due to conjugacy
  // omp_set_num_threads(8);
  int p = Y.cols();
  int b = X.cols();
  int n = X.rows();

  MatrixXd randn_theta = rstdnorm_mat(b,p);
  MatrixXd randn_e = rstdnorm_mat(n,p);

  MatrixXd coefs(b,p);

  for(int j = 0; j < p; j++){
    coefs.col(j) = sample_coefs_uncorrelated(Y.col(j), X, prior_mean.col(j), prior_prec.col(j), resid_prec.col(j), randn_theta.col(j),randn_e.col(j),b,n);
  }

  return(coefs);
}







// // // // // // // // // // #include<Eigen/SparseCholesky>
//
// // [[Rcpp::export()]]
// VectorXd res(VectorXd x){
//   int n = x.size();
//
//   #pragma omp parallel for shared(x, n)
//   for (size_t i=0; i<n; i++) {
//     for(int j = 0; j < 1e7; j++) {
//       x[i] += 1;
//     }
//   }
//   return(x);
// }
//
// // [[Rcpp::export()]]
// VectorXd res2(VectorXd x){
//   int n = x.size();
//
//   for (size_t i=0; i<n; i++) {
//     for(int j = 0; j < 1e7; j++) {
//       x[i] += 1;
//     }
//   }
//   return(x);
// }
//
//
// // // [[Rcpp::export()]]
// // void v1(MSpMat chol_R,MatrixXd W,VectorXd y,double tot_Eta_prec){
// //   SpMat chol_R_t = chol_R.transpose();
// //   MatrixXd RinvSqW = chol_R_t.triangularView<Lower>().solve(W);
// //   VectorXd WtRinvy = RinvSqW.transpose() * chol_R_t.triangularView<Lower>().solve(y) * tot_Eta_prec;
// // }
// // // [[Rcpp::export()]]
// // void v2(ArrayXd  resid_prec,MatrixXd X,VectorXd randn_e){
// //   VectorXd R_sq_diag = resid_prec.sqrt().inverse();
// //   MatrixXd RinvSqQtW = R_sq_diag.cwiseInverse().asDiagonal() * X;
// //   VectorXd eta_std = randn_e.array()/R_sq_diag.array();
// //   VectorXd WtURinvy = RinvSqQtW.transpose() * eta_std;
// // }
// //
// // // [[Rcpp::export()]]
// // VectorXd sample_MME_single_diagK2(
// //     VectorXd y,
// //     MatrixXd W,
// //     VectorXd prior_mean,
// //     VectorXd prior_prec,
// //     MSpMat chol_R,
// //     double tot_Eta_prec,
// //     VectorXd randn_theta,
// //     VectorXd randn_e,
// //     int b,
// //     int n
// // ){
// //
// //   chol_R *= 1/tot_Eta_prec;
// //   VectorXd theta_star = randn_theta.array() / prior_prec.cwiseSqrt().array();
// //   theta_star += prior_mean;
// //   VectorXd e_star = chol_R * randn_e;
// //   MatrixXd W_theta_star = W * theta_star;
// //   VectorXd y_resid = y - W_theta_star - e_star;
// //
// //   MatrixXd RinvSqW = chol_R.transpose().triangularView<Lower>().solve(W);
// //   VectorXd WtRinvy = RinvSqW.transpose() * chol_R.transpose().triangularView<Lower>().solve(y_resid);
// //
// //   VectorXd theta_tilda;
// //
// //   if(b < n) {
// //     MatrixXd C = RinvSqW.transpose() * RinvSqW;
// //     C.diagonal() += prior_prec;
// //     theta_tilda = C.llt().solve(WtRinvy);
// //   } else{
// //     MatrixXd VAi = W * prior_prec.cwiseInverse().asDiagonal();
// //     MatrixXd inner = VAi*W.transpose() + chol_R.transpose() * chol_R;
// //     VectorXd VAiWtURinvy = VAi * WtRinvy;
// //     VectorXd outerWtURinvy = VAi.transpose() * inner.ldlt().solve(VAiWtURinvy);
// //     theta_tilda = WtRinvy.array() / prior_prec.array();
// //     theta_tilda -= outerWtURinvy;
// //   }
// //
// //   VectorXd theta = theta_star + theta_tilda;
// //
// //   return theta;
// // }
// //
// // // [[Rcpp::export()]]
// // VectorXd sample_coefs_uncorrelated2(  // returns bx1
// //     VectorXd y,             // nx1
// //     MatrixXd X,             // nxb
// //     VectorXd prior_mean,    // bx1
// //     VectorXd prior_prec,    // bx1
// //     ArrayXd  resid_prec,    // nx1
// //     VectorXd randn_theta,   // bx1
// //     VectorXd randn_e,       // nx1
// //     int b,
// //     int n
// // ) {
// //   VectorXd R_sq_diag = resid_prec.sqrt().inverse();
// //   VectorXd theta_star = randn_theta.array()/prior_prec.array().sqrt();
// //   theta_star += prior_mean;
// //   VectorXd e_star = randn_e.array() * R_sq_diag.array();
// //   MatrixXd QtW_theta_star = X * theta_star;
// //   VectorXd eta_resid = y - QtW_theta_star - e_star;
// //   MatrixXd RinvSqQtW = R_sq_diag.cwiseInverse().asDiagonal() * X;
// //   VectorXd eta_std = eta_resid.array()/R_sq_diag.array();
// //   VectorXd WtURinvy = RinvSqQtW.transpose() * eta_std;
// //
// //   VectorXd theta_tilda;
// //   if(b < n) {
// //     MatrixXd C = RinvSqQtW.transpose() * RinvSqQtW;
// //     C.diagonal() += prior_prec;
// //     theta_tilda = C.llt().solve(WtURinvy);
// //   } else{
// //     MatrixXd VAi = X * prior_prec.cwiseInverse().asDiagonal();
// //     MatrixXd inner = VAi*X.transpose();
// //     inner.diagonal() += R_sq_diag.cwiseProduct(R_sq_diag);
// //     VectorXd VAiWtURinvy = VAi * WtURinvy;
// //     VectorXd outerWtURinvy = VAi.transpose() * inner.ldlt().solve(VAiWtURinvy);
// //     theta_tilda = WtURinvy.array() / prior_prec.array();
// //     theta_tilda -= outerWtURinvy;
// //   }
// //
// //   VectorXd coefs = theta_tilda + theta_star;
// //   return coefs;
// // }
