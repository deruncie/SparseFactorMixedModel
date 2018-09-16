#include <math.h>
#include <iostream>
#include "BSFG_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace RcppParallel;

// [[Rcpp::export()]]
MatrixXd matrix_multiply_toDense(SEXP X_, SEXP Y_){
  if(Rf_isNull(X_)) return(as<Map<MatrixXd> >(Y_));
  if(Rf_isMatrix(X_)) {
    Map<MatrixXd> X = as<Map<MatrixXd> >(X_);
    if(Rf_isMatrix(Y_)) {
      Map<MatrixXd> Y = as<Map<MatrixXd> >(Y_);
      if(X.cols() != Y.rows()) stop("Wrong dimensions of matrices");
      return(X*Y);
    } else{
      MSpMat Y = as<MSpMat>(Y_);
      if(X.cols() != Y.rows()) stop("Wrong dimensions of matrices");
      return(X*Y);
    }
  }
  else {
    MSpMat X = as<MSpMat>(X_);
    if(Rf_isMatrix(Y_)) {
      Map<MatrixXd> Y = as<Map<MatrixXd> >(Y_);
      if(X.cols() != Y.rows()) stop("Wrong dimensions of matrices");
      return(X*Y);
    } else{
      MSpMat Y = as<MSpMat>(Y_);
      if(X.cols() != Y.rows()) stop("Wrong dimensions of matrices");
      return(X*Y);
    }
  }
}


// [[Rcpp::depends(RcppZiggurat)]]
#include <Ziggurat.h>
#include <ZigguratR.h>
using namespace Rcpp;
static Ziggurat::Ziggurat::Ziggurat zigg;

// [[Rcpp::export()]]
MatrixXd rstdnorm_mat2(int n,int p) {  // returns nxp matrix
  MatrixXd X_mat = as<VectorXd>(rnorm(n*p));
  // Map<MatrixXd> X_mat(X_vec.data(),n,p);
  X_mat.resize(n,p);
  return(X_mat);
}

// [[Rcpp::export()]]
MatrixXd rstdnorm_mat3(int n,int p) {  // returns nxp matrix
  VectorXd X_vec(n*p);
  for(int i = 0; i < n*p; i++){
    X_vec[i] = zigg.norm();
  }
  MatrixXd X_mat = Map<MatrixXd>(X_vec.data(),n,p);
  return(X_mat);
}

static Ziggurat::R::ZigguratR ziggr;
// [[Rcpp::export()]]
MatrixXd rstdnorm_mat4(int n,int p) {  // returns nxp matrix
  VectorXd X_vec(n*p);
  for(int i = 0; i < n*p; i++){
    X_vec[i] = ziggr.norm();
  }
  MatrixXd X_mat = Map<MatrixXd>(X_vec.data(),n,p);
  return(X_mat);
}




// [[Rcpp::export()]]
unsigned long int myzigg_getSeed(){
  return(zigg.getSeed());
}
// [[Rcpp::export()]]
void myzigg_setSeed(unsigned long int s){
  zigg.setSeed(s);
}

// [[Rcpp::export()]]
void test_rstdnorm2(int n,int p) {
  MatrixXd x = rstdnorm_mat2(n,p);
}
// [[Rcpp::export()]]
void test_rstdnorm3(int n,int p) {
  MatrixXd x = rstdnorm_mat3(n,p);
}

// // [[Rcpp::export]]
// uint32_t test(uint32_t jsr) {
//   uint32_t jz = jsr;
//   jsr^=(jsr<<13);
//   jsr^=(jsr>>17);
//   jsr^=(jsr<<5);
//   return(jz+jsr);
// }


// [[Rcpp::export()]]
unsigned long int ce(unsigned long int i,unsigned long int j){
  i ^= j;
  return(i);
}

// -------------------------------------------- //
// ---------- sample_MME_fixedEffects --------- //
// -------------------------------------------- //

VectorXd sample_MME_single_diagK2(  // returns b x 1 vector
    const Ref<const VectorXd>& y,           // nx1
    const MatrixXd& X,                      // nxb
    const Ref<const VectorXd>& prior_mean,  // bx1
    const Ref<const VectorXd>& prior_prec,  // bx1
    const SEXP& chol_R_,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
    double tot_Eta_prec,                    // double
    const Ref<const VectorXd>& randn_theta, // bx1
    const Ref<const VectorXd>& randn_e      // 0x1 or nx1. 0x1 if b<n
){
  if(randn_e.size() == 0){
    MatrixXd RinvSqX;
    VectorXd XtRinvy;
    if(Rf_isMatrix(chol_R_)) {
      Map<MatrixXd> chol_R = as<Map<MatrixXd> >(chol_R_);
      RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
      XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y * sqrt(tot_Eta_prec));
    } else{
      MSpMat chol_R = as<MSpMat>(chol_R_);
      RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
      XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y * sqrt(tot_Eta_prec));
    }

    VectorXd XtRinvy_std_mu = XtRinvy + prior_prec.asDiagonal()*prior_mean;
    MatrixXd C = RinvSqX.transpose() * RinvSqX;
    C.diagonal() += prior_prec;
    LLT<MatrixXd> C_llt;
    C_llt.compute(C);
    MatrixXd chol_C = C_llt.matrixU();

    VectorXd b = chol_C.transpose().triangularView<Lower>().solve(XtRinvy_std_mu);
    b += randn_theta;
    b = chol_C.triangularView<Upper>().solve(b);
    return(b);
  } else {
    // Using algorithm from Bhattacharya et al 2016 Biometrika. https://academic.oup.com/biomet/article/103/4/985/2447851
    MatrixXd Phi;
    VectorXd alpha;
    if(Rf_isMatrix(chol_R_)) {
      Map<MatrixXd> chol_R = as<Map<MatrixXd> >(chol_R_);
      Phi = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
      alpha = chol_R.transpose().triangularView<Lower>().solve(y * sqrt(tot_Eta_prec));
    } else{
      MSpMat chol_R = as<MSpMat>(chol_R_);
      Phi = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
      alpha = chol_R.transpose().triangularView<Lower>().solve(y * sqrt(tot_Eta_prec));
    }

    VectorXd u = randn_theta.array() / prior_prec.cwiseSqrt().array();
    u += prior_mean;
    VectorXd v = Phi * u + randn_e;
    VectorXd alpha_v = alpha-v;

    MatrixXd D_PhiT = prior_prec.cwiseInverse().asDiagonal() * Phi.transpose();
    MatrixXd cov = Phi * D_PhiT;
    cov.diagonal().array() += 1.0;

    VectorXd w = cov.ldlt().solve(alpha_v);

    VectorXd theta = u + D_PhiT * w;

    return(theta);
  }
}


// RcppParallel struct for sampling coefficients from a set of parallel regression problems:
// Y = X %*% B + E, where
// b[i,j] ~ N(prior_mean[i,j],1/prior_prec[i,j])
// E[,j] ~ N(0,t(chol_R) %*% chol(R) / tot_Eta_prec[j])
// where chol_R is selected from a list based on h2s_index[j]
// Y is complete, so everythign has the same dimensions
struct sample_MME_single_diagK_worker2 : public RcppParallel::Worker {
  const Map<MatrixXd>& Y;
  const Map<MatrixXd>& X;
  const Map<MatrixXd>& prior_mean, prior_prec;
  const MatrixXd& randn_theta, randn_e;
  // const std::vector<MSpMat> chol_R_list;
  const Rcpp::List &chol_R_list;
  const VectorXi h2s_index;
  const VectorXd tot_Eta_prec;
  MatrixXd &coefs;

  sample_MME_single_diagK_worker2(
    const Map<MatrixXd> &Y,           // nxp
    const Map<MatrixXd>& X,           // nxb
    const Map<MatrixXd>& prior_mean,  // bxp
    const Map<MatrixXd>& prior_prec,  // bxp    // const std::vector<MSpMat> &chol_R_list, // each element a MSpMat nxn upper-triangular cholesky decomposition of R
    const Rcpp::List &chol_R_list,
    const VectorXi h2s_index,   // px1, 1-based index
    const VectorXd tot_Eta_prec,// px1
    const MatrixXd& randn_theta, // bxp
    const MatrixXd& randn_e,     // 0xp or nxp. The former if b<n, the latter if b >= n
    MatrixXd &coefs       // bxp
  ):
    Y(Y), X(X), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
    chol_R_list(chol_R_list), h2s_index(h2s_index), tot_Eta_prec(tot_Eta_prec), coefs(coefs) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++){
      int h2_index = h2s_index[j] - 1;
      const Rcpp::List Sigma_Choleskys_i = chol_R_list[h2_index];
      const SEXP chol_R = Sigma_Choleskys_i["chol_Sigma"];
      coefs.col(j) = sample_MME_single_diagK2(Y.col(j), X, prior_mean.col(j), prior_prec.col(j), chol_R, tot_Eta_prec[j], randn_theta.col(j),randn_e.col(j));
    }
  }
};


// -------------------------------------------------- //
// -- Versions of the independent residuals regression model --- //
// -------------------------------------------------- //


// basic version - all regression have same dimensions
// [[Rcpp::export()]]
MatrixXd sample_MME_fixedEffects_c2(  // returns bxp matrix
    Map<MatrixXd> Y,              // nxp
    Map<MatrixXd> X,              // nxb
    Rcpp::List chol_R_list,   // list of Cholesky decompositions of residuals. nxn upper triangular, dgCMatrix
    VectorXi h2s_index,           // px1 index of Cholesky matrix for each column
    Map<VectorXd> tot_Eta_prec,   // px1
    Map<MatrixXd> prior_mean,     // bxp
    Map<MatrixXd> prior_prec,     // bxp
    int grainSize) {

  int b = X.cols();
  int p = Y.cols();
  int n = Y.rows();

  // Map<MatrixXd> randn_theta(as<VectorXd>(rnorm(b*p)).data(),b,p);
  // // MatrixXd randn_theta = rstdnorm_mat2(b,p);
  // Map<MatrixXd> randn_e(as<VectorXd>(rnorm(0)).data(),0,p);
  // if(b < n) {
  //   // Rcout << "b<n" << std::endl;
  //   randn_e = rstdnorm_mat2(0,p);
  // } else{
  //   // Rcout << "b>=n" << std::endl;
  //   new (&randn_e) Map<MatrixXd>(as<VectorXd>(rnorm(n*p)).data(),n,p);
  // }
  MatrixXd randn_theta = rstdnorm_mat2(b,p);
  MatrixXd randn_e;
  if(b < n) {
    randn_e = rstdnorm_mat2(0,p);
  } else{
    randn_e = rstdnorm_mat2(n,p);
  }
  if(randn_theta.rows() != b || randn_theta.cols() != p) stop("wrong dimension of randn_theta");
  if(randn_e.cols() != p) stop("wrong dimension of randn_e");

  // std::vector<MSpMat> chol_R_list;
  // for(int i = 0; i < h2s_index.maxCoeff(); i++){
  //   Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);
  //   chol_R_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
  // }

  MatrixXd coefs(b,p);

  sample_MME_single_diagK_worker2 sampler(Y,X,prior_mean,prior_prec,chol_R_list,h2s_index,tot_Eta_prec,randn_theta,randn_e, coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return(coefs);
}

struct sample_MME_single_diagK_worker3 : public RcppParallel::Worker {
  MatrixXd Y;
  MatrixXd X;
  MatrixXd prior_mean, prior_prec;
  MatrixXd randn_theta, randn_e;
  // const std::vector<MSpMat> chol_R_list;
  Rcpp::List &chol_R_list;
  VectorXi h2s_index;
  VectorXd tot_Eta_prec;
  MatrixXd &coefs;

  sample_MME_single_diagK_worker3(
    MatrixXd Y,           // nxp
    MatrixXd X,           // nxb
    MatrixXd prior_mean,  // bxp
    MatrixXd prior_prec,  // bxp    // const std::vector<MSpMat> &chol_R_list, // each element a MSpMat nxn upper-triangular cholesky decomposition of R
    Rcpp::List &chol_R_list,
    VectorXi h2s_index,   // px1, 1-based index
    VectorXd tot_Eta_prec,// px1
    MatrixXd randn_theta, // bxp
    MatrixXd randn_e,     // 0xp or nxp. The former if b<n, the latter if b >= n
    MatrixXd &coefs       // bxp
  ):
    Y(Y), X(X), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
    chol_R_list(chol_R_list), h2s_index(h2s_index), tot_Eta_prec(tot_Eta_prec), coefs(coefs) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++){
      int h2_index = h2s_index[j] - 1;
      Rcpp::List Sigma_Choleskys_i = chol_R_list[h2_index];
      SEXP chol_R = Sigma_Choleskys_i["chol_Sigma"];
      coefs.col(j) = sample_MME_single_diagK2(Y.col(j), X, prior_mean.col(j), prior_prec.col(j), chol_R, tot_Eta_prec[j], randn_theta.col(j),randn_e.col(j));
    }
  }
};


// -------------------------------------------------- //
// -- Versions of the independent residuals regression model --- //
// -------------------------------------------------- //


// basic version - all regression have same dimensions
// [[Rcpp::export()]]
MatrixXd sample_MME_fixedEffects_c3(  // returns bxp matrix
    Map<MatrixXd> Y,              // nxp
    Map<MatrixXd> X,              // nxb
    Rcpp::List chol_R_list,   // list of Cholesky decompositions of residuals. nxn upper triangular, dgCMatrix
    VectorXi h2s_index,           // px1 index of Cholesky matrix for each column
    Map<VectorXd> tot_Eta_prec,   // px1
    Map<MatrixXd> prior_mean,     // bxp
    Map<MatrixXd> prior_prec,     // bxp
    int grainSize) {

  int b = X.cols();
  int p = Y.cols();
  int n = Y.rows();

  // Map<MatrixXd> randn_theta(as<VectorXd>(rnorm(b*p)).data(),b,p);
  // // MatrixXd randn_theta = rstdnorm_mat2(b,p);
  // Map<MatrixXd> randn_e(as<VectorXd>(rnorm(0)).data(),0,p);
  // if(b < n) {
  //   // Rcout << "b<n" << std::endl;
  //   randn_e = rstdnorm_mat2(0,p);
  // } else{
  //   // Rcout << "b>=n" << std::endl;
  //   new (&randn_e) Map<MatrixXd>(as<VectorXd>(rnorm(n*p)).data(),n,p);
  // }
  MatrixXd randn_theta = rstdnorm_mat4(b,p);
  MatrixXd randn_e;
  if(b < n) {
    randn_e = rstdnorm_mat4(0,p);
  } else{
    randn_e = rstdnorm_mat4(n,p);
  }
  if(randn_theta.rows() != b || randn_theta.cols() != p) stop("wrong dimension of randn_theta");
  if(randn_e.cols() != p) stop("wrong dimension of randn_e");

  // std::vector<MSpMat> chol_R_list;
  // for(int i = 0; i < h2s_index.maxCoeff(); i++){
  //   Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);
  //   chol_R_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
  // }

  MatrixXd coefs(b,p);

  sample_MME_single_diagK_worker3 sampler(Y,X,prior_mean,prior_prec,chol_R_list,h2s_index,tot_Eta_prec,randn_theta,randn_e, coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return(coefs);
}




//// New strategy

VectorXd sample_MME_single_diagK4(  // returns b x 1 vector
    const Ref<const VectorXd>& y,           // nx1
    const MatrixXd& RinvSqX,                // nxb
    const MatrixXd& C_,                     // bxb
    const Ref<const VectorXd>& prior_mean,  // bx1
    const Ref<const VectorXd>& prior_prec,  // bx1
    const SEXP& chol_R_,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
    double tot_Eta_prec,                    // double
    const Ref<const VectorXd>& randn_theta, // bx1
    const Ref<const VectorXd>& randn_e      // 0x1 or nx1. 0x1 if b<n
){
  if(randn_e.size() == 0){
    // MatrixXd RinvSqX;
    VectorXd XtRinvy;
    // RinvSqX
    if(Rf_isMatrix(chol_R_)) {
      Map<MatrixXd> chol_R = as<Map<MatrixXd> >(chol_R_);
      // RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
      XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y * tot_Eta_prec);
    } else{
      MSpMat chol_R = as<MSpMat>(chol_R_);
      // RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
      XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y * tot_Eta_prec);
    }

    VectorXd XtRinvy_std_mu = XtRinvy + prior_prec.asDiagonal()*prior_mean;
    // MatrixXd C = RinvSqX.transpose() * RinvSqX * tot_Eta_prec;
    MatrixXd C = C_ * tot_Eta_prec;
    C.diagonal() += prior_prec;
    LLT<MatrixXd> C_llt;
    C_llt.compute(C);
    MatrixXd chol_C = C_llt.matrixU();

    VectorXd b = chol_C.transpose().triangularView<Lower>().solve(XtRinvy_std_mu);
    b += randn_theta;
    b = chol_C.triangularView<Upper>().solve(b);
    return(b);
  } else {
    // Using algorithm from Bhattacharya et al 2016 Biometrika. https://academic.oup.com/biomet/article/103/4/985/2447851
    MatrixXd Phi = RinvSqX * sqrt(tot_Eta_prec);
    VectorXd alpha;
    if(Rf_isMatrix(chol_R_)) {
      Map<MatrixXd> chol_R = as<Map<MatrixXd> >(chol_R_);
      // Phi = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
      alpha = chol_R.transpose().triangularView<Lower>().solve(y * sqrt(tot_Eta_prec));
    } else{
      MSpMat chol_R = as<MSpMat>(chol_R_);
      // Phi = chol_R.transpose().triangularView<Lower>().solve(X * sqrt(tot_Eta_prec));
      alpha = chol_R.transpose().triangularView<Lower>().solve(y * sqrt(tot_Eta_prec));
    }

    VectorXd u = randn_theta.array() / prior_prec.cwiseSqrt().array();
    u += prior_mean;
    VectorXd v = Phi * u + randn_e;
    VectorXd alpha_v = alpha-v;

    MatrixXd D_PhiT = prior_prec.cwiseInverse().asDiagonal() * Phi.transpose();
    MatrixXd cov = Phi * D_PhiT;
    cov.diagonal().array() += 1.0;

    VectorXd w = cov.ldlt().solve(alpha_v);

    VectorXd theta = u + D_PhiT * w;

    return(theta);
  }
}


struct sample_MME_single_diagK_worker4 : public RcppParallel::Worker {
  MatrixXd Y;
  MatrixXd X;
  MatrixXd prior_mean, prior_prec;
  MatrixXd randn_theta, randn_e;
  // const std::vector<MSpMat> chol_R_list;
  Rcpp::List &chol_R_list;
  VectorXi h2s_index;
  VectorXd tot_Eta_prec;
  MatrixXd &coefs;

  sample_MME_single_diagK_worker4(
    MatrixXd Y,           // nxp
    MatrixXd X,           // nxb
    MatrixXd prior_mean,  // bxp
    MatrixXd prior_prec,  // bxp    // const std::vector<MSpMat> &chol_R_list, // each element a MSpMat nxn upper-triangular cholesky decomposition of R
    Rcpp::List &chol_R_list,
    VectorXi h2s_index,   // px1, 1-based index
    VectorXd tot_Eta_prec,// px1
    MatrixXd randn_theta, // bxp
    MatrixXd randn_e,     // 0xp or nxp. The former if b<n, the latter if b >= n
    MatrixXd &coefs       // bxp
  ):
    Y(Y), X(X), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
    chol_R_list(chol_R_list), h2s_index(h2s_index), tot_Eta_prec(tot_Eta_prec), coefs(coefs) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++){
      int h2_index = i;
      int num_columns = 0;
      MatrixXd RinvSqX;
      MatrixXd C;
      SEXP chol_R_;
      for(int j = 0; j < h2s_index.size(); j++){
        if(h2s_index[j] == h2_index) {
          if(num_columns == 0) {
            num_columns = 1;
            Rcpp::List Sigma_Choleskys_i = chol_R_list[h2_index - 1];
            chol_R_ = Sigma_Choleskys_i["chol_Sigma"];
            if(Rf_isMatrix(chol_R_)) {
              Map<MatrixXd> chol_R = as<Map<MatrixXd> >(chol_R_);
              RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
            } else{
              MSpMat chol_R = as<MSpMat>(chol_R_);
              RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
            }
            C = RinvSqX.transpose() * RinvSqX;
          }
          coefs.col(j) = sample_MME_single_diagK4(Y.col(j), RinvSqX,C, prior_mean.col(j), prior_prec.col(j), chol_R_, tot_Eta_prec[j], randn_theta.col(j),randn_e.col(j));
        }
      }
    }
  }
};



// basic version - all regression have same dimensions
// [[Rcpp::export()]]
MatrixXd sample_MME_fixedEffects_c4(  // returns bxp matrix
    Map<MatrixXd> Y,              // nxp
    Map<MatrixXd> X,              // nxb
    Rcpp::List chol_R_list,   // list of Cholesky decompositions of residuals. nxn upper triangular, dgCMatrix
    VectorXi h2s_index,           // px1 index of Cholesky matrix for each column
    Map<VectorXd> tot_Eta_prec,   // px1
    Map<MatrixXd> prior_mean,     // bxp
    Map<MatrixXd> prior_prec,     // bxp
    int grainSize) {

  int b = X.cols();
  int p = Y.cols();
  int n = Y.rows();

  MatrixXd randn_theta = rstdnorm_mat4(b,p);
  MatrixXd randn_e;
  if(b < n) {
    randn_e = rstdnorm_mat4(0,p);
  } else{
    randn_e = rstdnorm_mat4(n,p);
  }
  if(randn_theta.rows() != b || randn_theta.cols() != p) stop("wrong dimension of randn_theta");
  if(randn_e.cols() != p) stop("wrong dimension of randn_e");

  MatrixXd coefs(b,p);

  sample_MME_single_diagK_worker4 sampler(Y,X,prior_mean,prior_prec,chol_R_list,h2s_index,tot_Eta_prec,randn_theta,randn_e, coefs);
  RcppParallel::parallelFor(h2s_index.minCoeff(),h2s_index.maxCoeff()+1,sampler,grainSize);
  return(coefs);
}

// new2

struct sample_MME_single_diagK_worker5 : public RcppParallel::Worker {
  const MatrixXd& Y;
  const MatrixXd& X;
  const MatrixXd& RinvSqX;
  const MatrixXd& C;
  SEXP chol_R_;
  const IntegerVector& cols;
  const MatrixXd& prior_mean, prior_prec;
  const MatrixXd& randn_theta, randn_e;
  // const std::vector<MSpMat> chol_R_list;
  // Rcpp::List &chol_R_list;
  // VectorXi h2s_index;
  const VectorXd& tot_Eta_prec;
  MatrixXd &coefs;

  sample_MME_single_diagK_worker5(
    const MatrixXd& Y,           // nxp
    const MatrixXd& X,           // nxb
    const MatrixXd& RinvSqX,
    const MatrixXd& C,
    SEXP chol_R_,
    const IntegerVector& cols,
    const MatrixXd& prior_mean,  // bxp
    const MatrixXd& prior_prec,  // bxp    // const std::vector<MSpMat> &chol_R_list, // each element a MSpMat nxn upper-triangular cholesky decomposition of R
    const VectorXd& tot_Eta_prec,// px1
    const MatrixXd& randn_theta, // bxp
    const MatrixXd& randn_e,     // 0xp or nxp. The former if b<n, the latter if b >= n
    MatrixXd &coefs       // bxp
  ):
    Y(Y), X(X), RinvSqX(RinvSqX), C(C), chol_R_(chol_R_), cols(cols),prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
    tot_Eta_prec(tot_Eta_prec), coefs(coefs) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t i = begin; i < end; i++){
      int j = cols[i];
      coefs.col(j) = sample_MME_single_diagK4(Y.col(j), RinvSqX,C, prior_mean.col(j), prior_prec.col(j), chol_R_, tot_Eta_prec[j], randn_theta.col(j),randn_e.col(j));
    }
  }
};


// [[Rcpp::export]]
Rcpp::IntegerVector which2(Rcpp::LogicalVector x) {
  Rcpp::IntegerVector v = Rcpp::seq(0, x.size()-1);
  return v[x];
}

// basic version - all regression have same dimensions
// [[Rcpp::export()]]
MatrixXd sample_MME_fixedEffects_c5(  // returns bxp matrix
    Map<MatrixXd> Y,              // nxp
    Map<MatrixXd> X,              // nxb
    Rcpp::List chol_R_list,   // list of Cholesky decompositions of residuals. nxn upper triangular, dgCMatrix
    IntegerVector h2s_index,           // px1 index of Cholesky matrix for each column
    Map<VectorXd> tot_Eta_prec,   // px1
    Map<MatrixXd> prior_mean,     // bxp
    Map<MatrixXd> prior_prec,     // bxp
    int grainSize) {

  int b = X.cols();
  int p = Y.cols();
  int n = Y.rows();

  MatrixXd randn_theta = rstdnorm_mat4(b,p);
  MatrixXd randn_e;
  if(b < n) {
    randn_e = rstdnorm_mat4(0,p);
  } else{
    randn_e = rstdnorm_mat4(n,p);
  }
  if(randn_theta.rows() != b || randn_theta.cols() != p) stop("wrong dimension of randn_theta");
  if(randn_e.cols() != p) stop("wrong dimension of randn_e");

  MatrixXd coefs(b,p);

  for(int i = min(h2s_index); i <= max(h2s_index); i++) {
    int h2_index = i;
    IntegerVector cols = which2(h2s_index == h2_index);
    if(cols.size() > 0){
      MatrixXd RinvSqX;
      MatrixXd C;
      Rcpp::List Sigma_Choleskys_i = chol_R_list[h2_index - 1];
      SEXP chol_R_ = Sigma_Choleskys_i["chol_Sigma"];
      if(Rf_isMatrix(chol_R_)) {
        Map<MatrixXd> chol_R = as<Map<MatrixXd> >(chol_R_);
        RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
      } else{
        MSpMat chol_R = as<MSpMat>(chol_R_);
        RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
      }
      C = RinvSqX.transpose() * RinvSqX;
      sample_MME_single_diagK_worker5 sampler(Y,X,RinvSqX,C,chol_R_,cols,prior_mean,prior_prec,tot_Eta_prec,randn_theta,randn_e, coefs);
      RcppParallel::parallelFor(0,cols.size(),sampler,grainSize);
    }
  }
  return(coefs);
}

