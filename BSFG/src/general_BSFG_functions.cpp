#include <math.h>
#include <iostream>
#include "BSFG_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace RcppParallel;


// -------------------------------------------- //
// ---------- helper functions --------- //
// -------------------------------------------- //
// functions to speed up sparse multiplication and conversion to dense matrices

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



// [[Rcpp::export()]]
MatrixXd rstdnorm_mat(int n,int p) {  // returns nxp matrix
  VectorXd X_vec(n*p);
  for(int i = 0; i < n*p; i++){
    X_vec[i] = ziggr.norm();
  }
  MatrixXd X_mat = Map<MatrixXd>(X_vec.data(),n,p);
  return(X_mat);
}



// [[Rcpp::export()]]
VectorXd find_candidate_states(
    MatrixXd h2s_matrix,
    double step_size,
    int old_state
) {
  VectorXd dists = (h2s_matrix.colwise() - h2s_matrix.col(old_state)).cwiseAbs().colwise().sum();
  VectorXd indices(dists.size());
  int count = 0;
  for(int i = 0; i < dists.size(); i++){
    if(dists[i] < step_size & dists[i] > 0) {
      indices[count] = i;
      count++;
    }
  }
  if(count == 0) {  // return all indices as candidates
    for(int i = 0; i < dists.size(); i++){
      indices[count] = i;
      count++;
    }
  }
  return indices.head(count);
}


// code to convert list of R matrices (sparse or dense) into a thread-safe object
struct R_matrix {
  Map<MatrixXd> dense;
  MSpMat sparse;
  bool isDense;
  R_matrix(Map<MatrixXd> dense_, MSpMat sparse_,bool isDense_) : dense(dense_), sparse(sparse_), isDense(isDense_) {}
};

void load_R_matrices_list(const Rcpp::List X_list, std::vector<R_matrix>& X_vector){
  // null_matrices
  MatrixXd null_d = MatrixXd::Zero(0,0);
  Map<MatrixXd> M_null_d(null_d.data(),0,0);
  SpMat null_s = null_d.sparseView();
  MSpMat M_null_s(0,0,0,null_s.outerIndexPtr(),null_s.innerIndexPtr(),null_s.valuePtr());

  int p = X_list.size();
  X_vector.reserve(p);
  for(int i = 0; i < p; i++){
    SEXP Xi_ = X_list[i];
    if(Rf_isMatrix(Xi_)){
      Map<MatrixXd> Xi = as<Map<MatrixXd> >(Xi_);
      R_matrix Xim(Xi,M_null_s,true);
      X_vector.push_back(Xim);
    } else{
      MSpMat Xi = as<MSpMat>(Xi_);
      R_matrix Xim(M_null_d,Xi,false);
      X_vector.push_back(Xim);
    }
  }
}

// -------------------------------------------- //
// ---------- regression_sampler --------- //
// -------------------------------------------- //

Rcpp::IntegerVector which(Rcpp::LogicalVector x) {
  Rcpp::IntegerVector v = Rcpp::seq(0, x.size()-1);
  return v[x];
}
MatrixXd get_RinvSqX(const R_matrix& chol_R, MatrixXd X){
  MatrixXd RinvSqX;
  if(chol_R.isDense) {
    RinvSqX = chol_R.dense.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
  } else{
    RinvSqX = chol_R.sparse.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
  }
  return(RinvSqX);
}

VectorXd regression_sampler_v1(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b < n
    const Ref<const VectorXd>& y,           // nx1
    const MatrixXd& W,           // nxa
    const MatrixXd& RinvSqX,                // nxb
    const MatrixXd& C,                     // bxb
    const VectorXd& prior_prec_alpha, // ax 1
    const Ref<const VectorXd>& prior_mean_beta,  // bx1
    const Ref<const VectorXd>& prior_prec_beta,  // bx1
    const R_matrix& chol_R,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
    double Y_prec,                    // double
    const VectorXd& randn_alpha,
    const Ref<const VectorXd>& randn_beta,
    const double rgamma_1,
    const double Y_prec_b0
){
  int n = y.size();
  int a = W.cols();
  int b = RinvSqX.cols();

  // Check inputs
  if(W.rows() != n) stop("Wrong dimension of W");
  if(RinvSqX.rows() != n) stop("Wrong dimension of X");
  if(C.rows() != b || C.cols() != b) stop("Wrong dimension of C");
  if(prior_prec_alpha.size() != a) stop("Wrong length of prior_prec_alpha");
  if(prior_mean_beta.size() != b) stop("Wrong length of prior_mean_beta");
  if(prior_prec_beta.size() != b) stop("Wrong length of prior_prec_beta");
  if(randn_alpha.size() != a) stop("Wrong length of randn_alpha");
  if(randn_beta.size() != b) stop("Wrong length of randn_beta");

  // Calculate cholesky of A_beta
  // C = Xt(RtR)^-1X
  MatrixXd C_beta = C;
  C_beta.diagonal() += prior_prec_beta;
  LLT<MatrixXd> A_beta_llt;
  A_beta_llt.compute(C_beta);
  MatrixXd chol_A_beta = A_beta_llt.matrixU();
  // Y_prec * chol_A_beta^\T * chol_A_beta = A_beta

  // Step 1
  VectorXd alpha(a);
  VectorXd y_tilde = y;
  if(a > 0) {
    // Sample alpha
    // Calculate A_alpha = Y_prec*W^T*Sigma_beta^{-1}*W + D_alpha^{-1}
    // We don't need to actually calculate Sigma_beta^{-1} directly.
    MatrixXd RinvSqW = get_RinvSqX(chol_R,W);  // n*n*a -> n x a
    MatrixXd WtRinvX = RinvSqW.transpose() * RinvSqX; // a*n*b -> a*b
    MatrixXd invSqAbXtRinvW = chol_A_beta.transpose().triangularView<Lower>().solve(WtRinvX.transpose()); // b*b*a -> b x a

    MatrixXd A_alpha = Y_prec * (RinvSqW.transpose() * RinvSqW - invSqAbXtRinvW.transpose() * invSqAbXtRinvW);
    A_alpha.diagonal() += prior_prec_alpha;

    LLT<MatrixXd> A_alpha_llt;
    A_alpha_llt.compute(A_alpha);
    MatrixXd chol_A_alpha = A_alpha_llt.matrixU();

    VectorXd Rinvsqy = get_RinvSqX(chol_R,y); // n*n*q -> n x 1;
    VectorXd XtRinvy = RinvSqX.transpose() * Rinvsqy; // b*n*1 >- b x 1
    VectorXd invSqAbXtRinvy = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy); // b*b*1 -> b*1

    VectorXd WtSbinvy = RinvSqW.transpose() * Rinvsqy - invSqAbXtRinvW.transpose() * invSqAbXtRinvy;

    alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(WtSbinvy) * Y_prec + randn_alpha;
    alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
    y_tilde = y - W * alpha;
  }

  // Step 2 - sample Y_prec
  // We don't need to actually calculate Sigma_beta^{-1} directly.
  VectorXd RinvSqy = get_RinvSqX(chol_R,y_tilde);
  VectorXd XtRinvy = RinvSqX.transpose() * RinvSqy;
  VectorXd prod1 = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy);
  double score = Y_prec_b0 + (RinvSqy.dot(RinvSqy) - prod1.dot(prod1))/2;
  Y_prec = rgamma_1/score;

  // Step 3 - sample beta
  VectorXd XtRinvy_std_mu = XtRinvy*Y_prec + prior_prec_beta.asDiagonal()*prior_mean_beta;
  VectorXd beta = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy_std_mu) / sqrt(Y_prec) + randn_beta;
  beta = chol_A_beta.triangularView<Upper>().solve(beta) / sqrt(Y_prec);

  VectorXd result(1+a+b);
  result << Y_prec,alpha,beta;

  return(result);
}


VectorXd regression_sampler_v2(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b > n
    const Ref<const VectorXd>& y,           // nx1
    const MatrixXd& W,           // nxa
    const MatrixXd& X,           // nxm or nxb
    const VectorXd& prior_prec_alpha, // ax 1
    const Ref<const VectorXd>& prior_mean_beta,  // bx1
    const Ref<const VectorXd>& prior_prec_beta,  // bx1
    const R_matrix& chol_R,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
    const MatrixXd& R,
    double Y_prec,                    // double
    const VectorXd& randn_alpha,
    const Ref<const VectorXd>& randn_beta,
    const Ref<const VectorXd>& randn_e,
    const double rgamma_1,
    const double Y_prec_b0
){
  int n = y.size();
  int a = W.cols();
  int b = X.cols();

  // Check inputs
  if(W.rows() != n) stop("Wrong dimension of W");
  if(X.rows() != n) stop("Wrong dimension of X");
  if(prior_prec_alpha.size() != a) stop("Wrong length of prior_prec_alpha");
  if(prior_mean_beta.size() != b) stop("Wrong length of prior_mean_beta");
  if(prior_prec_beta.size() != b) stop("Wrong length of prior_prec_beta");
  if(randn_alpha.size() != a) stop("Wrong length of randn_alpha");
  if(randn_beta.size() != b) stop("Wrong length of randn_beta");

  // Calculate inverse of Sigma_beta
  MatrixXd DXt = prior_prec_beta.cwiseInverse().asDiagonal() * X.transpose();
  MatrixXd Sigma_beta = X * DXt + R;
  LDLT<MatrixXd> Sigma_beta_ldlt;
  Sigma_beta_ldlt.compute(Sigma_beta);

  // Step 1
  VectorXd alpha(a);
  VectorXd y_tilde = y;
  if(a > 0) {
    // Sample alpha
    MatrixXd SbinvW = Sigma_beta_ldlt.solve(W);
    MatrixXd A_alpha = Y_prec * SbinvW.transpose() * W;
    A_alpha.diagonal() += prior_prec_alpha;

    LLT<MatrixXd> A_alpha_llt;
    A_alpha_llt.compute(A_alpha);
    MatrixXd chol_A_alpha = A_alpha_llt.matrixU();

    VectorXd WtSbinvy = SbinvW.transpose() * y;
    alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(WtSbinvy) * Y_prec + randn_alpha;
    alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
    y_tilde = y - W * alpha;
  }

  // Step 2 - sample Y_prec
  VectorXd e2 = y_tilde.transpose() * Sigma_beta_ldlt.solve(y_tilde);
  double score = Y_prec_b0 + e2[0]/2;
  Y_prec = rgamma_1/score;

  // Step 3 - sample beta
  // what about prior mean?
  VectorXd u = randn_beta.array() / (prior_prec_beta * Y_prec).cwiseSqrt().array() + prior_mean_beta.array();
  VectorXd v = sqrt(Y_prec) * X * u;
  if(chol_R.isDense) {
    v += chol_R.dense.transpose().triangularView<Lower>() * randn_e;
  } else{
    v += chol_R.sparse.transpose().triangularView<Lower>() * randn_e;
  }
  VectorXd w = Sigma_beta_ldlt.solve(y_tilde * sqrt(Y_prec) - v);
  VectorXd beta = u + DXt * w / sqrt(Y_prec);

  VectorXd result(1+a+b);
  result << Y_prec,alpha,beta;

  return(result);
}


VectorXd regression_sampler_v3(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b > n > m
    const Ref<const VectorXd>& y,           // nx1
    const MatrixXd& W,           // nxa
    const MatrixXd& U,           // nxm or nxb
    const MatrixXd& V,           // mxb
    const VectorXd& prior_prec_alpha, // ax 1
    const Ref<const VectorXd>& prior_mean_beta,  // bx1
    const Ref<const VectorXd>& prior_prec_beta,  // bx1
    const R_matrix& chol_R,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
    const MatrixXd& Rinv,
    const MatrixXd& RinvU,
    const MatrixXd& UtRinvU,
    double Y_prec,                    // double
    const VectorXd& randn_alpha,
    const Ref<const VectorXd>& randn_beta,
    const Ref<const VectorXd>& randn_e,
    const double rgamma_1,
    const double Y_prec_b0
){
  int n = y.size();
  int a = W.cols();
  if(V.rows() != U.cols()) stop("Wrong dimensions of V");
  MatrixXd X = U*V;
  int b = X.cols();

  // Check inputs
  if(W.rows() != n) stop("Wrong dimension of W");
  if(X.rows() != n) stop("Wrong dimension of X");
  if(prior_prec_alpha.size() != a) stop("Wrong length of prior_prec_alpha");
  if(prior_mean_beta.size() != b) stop("Wrong length of prior_mean_beta");
  if(prior_prec_beta.size() != b) stop("Wrong length of prior_prec_beta");
  if(randn_alpha.size() != a) stop("Wrong length of randn_alpha");
  if(randn_beta.size() != b) stop("Wrong length of randn_beta");
  if(randn_e.size() != n) stop("Wrong length of randn_e");

  // Calculate inverse of Sigma_beta
  // MatrixXd Sigma_beta_inv;
  // Using Ainv - Ainv * U * (I + BVAinvU)inv * BVAinv in case B = VDVt is singular
  MatrixXd DVt = prior_prec_beta.cwiseInverse().asDiagonal() * V.transpose();
  MatrixXd VDVt = V * DVt;
  if(RinvU.rows() != n) stop("Wrong dimensions of RinvU");
  MatrixXd inner = VDVt * UtRinvU;
  inner.diagonal().array() += 1.0;
  LDLT<MatrixXd> inner_ldlt;
  inner_ldlt.compute(inner);
  // Sigma_beta_inv = Rinv - RinvU * inner.ldlt().solve(VDVt * RinvU.transpose());  // Don't actually calculate this. Stay in mxm space

  // Step 1
  VectorXd alpha(a);
  VectorXd y_tilde = y;
  if(a > 0) {
    // Sample alpha
    // MatrixXd SbinvW = Sigma_beta_inv * W;
    // MatrixXd A_alpha = Y_prec * SbinvW.transpose() * W;
    MatrixXd RinvW = Rinv * W;
    MatrixXd UtRinvW = U.transpose() * RinvW;
    MatrixXd A_alpha = Y_prec * (W.transpose() * RinvW - UtRinvW.transpose() * inner_ldlt.solve(VDVt * UtRinvW));
    A_alpha.diagonal() += prior_prec_alpha;

    LLT<MatrixXd> A_alpha_llt;
    A_alpha_llt.compute(A_alpha);
    MatrixXd chol_A_alpha = A_alpha_llt.matrixU();

    // VectorXd WtSbinvy = SbinvW.transpose() * y;
    VectorXd UtRinvy = RinvU.transpose() * y;
    VectorXd WtSbinvy = RinvW.transpose() * y - UtRinvW.transpose() * inner_ldlt.solve(VDVt * UtRinvy);

    alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(WtSbinvy) * Y_prec + randn_alpha;
    alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
    y_tilde = y - W * alpha;
  }

  // Step 2 - sample Y_prec
  // VectorXd e2 = y_tilde.transpose() * Sigma_beta_inv * y_tilde;
  VectorXd Rinv_y = Rinv * y_tilde;
  VectorXd UtRinvy = RinvU.transpose() * y_tilde;
  VectorXd e2 = y_tilde.transpose() * Rinv_y - UtRinvy.transpose() * inner_ldlt.solve(VDVt * UtRinvy);

  double score = Y_prec_b0 + e2[0]/2;
  Y_prec = rgamma_1/score;

  // Step 3 - sample beta
  // what about prior mean?
  VectorXd u = randn_beta.array() / (prior_prec_beta * Y_prec).cwiseSqrt().array() + prior_mean_beta.array();
  VectorXd v = sqrt(Y_prec) * X * u;
  if(chol_R.isDense) {
    v += chol_R.dense.transpose().triangularView<Lower>() * randn_e;
  } else{
    v += chol_R.sparse.transpose().triangularView<Lower>() * randn_e;
  }
  // VectorXd w = Sigma_beta_inv * (y_tilde * sqrt(Y_prec) - v);
  VectorXd e = y_tilde * sqrt(Y_prec) - v;
  VectorXd UtRinve = RinvU.transpose() * e;
  VectorXd w = Rinv * e - RinvU * inner_ldlt.solve(VDVt * UtRinve);

  VectorXd beta = u + DVt * (U.transpose() * w) / sqrt(Y_prec); //b*b*1 + b*n*1

  VectorXd result(1+a+b);
  result << Y_prec,alpha,beta;

  return(result);
}



struct regression_sampler_worker : public RcppParallel::Worker {
  const int sampler;  // which sampler to use?
  const VectorXi& trait_set;
  const Map<MatrixXd> Y;           // nx1
  const Map<MatrixXd> W_base;           // nxa
  const std::vector<R_matrix>& W_list;           // nxa
  const Map<MatrixXd> X_U;   // could be X or U.
  const Map<MatrixXd> V;
  const MatrixXd RinvSqX;                // nxb
  const MatrixXd& C;                     // bxb
  const R_matrix& chol_R;                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
  const MatrixXd& R;
  const MatrixXd& Rinv;
  const MatrixXd& RinvU;
  const MatrixXd& UtRinvU;
  const Map<MatrixXd> prior_prec_alpha1; // a1 x p matrix of prior precisions for alpha1
  const VectorXd& prior_prec_alpha2;     // p-vector of precision of alpha2s for each trait
  const Map<MatrixXd> prior_mean_beta; // b x p matrix of prior means of beta
  const Map<MatrixXd> prior_prec_beta; // b x p matrix of prior precisions of beta
  double Y_prec_b0;

  const MatrixXd& randn_alpha1;
  const std::vector<VectorXd>& randn_alpha2;
  const MatrixXd& randn_beta;
  const MatrixXd& randn_e;
  const VectorXd& rgamma_1;

  MatrixXd& alpha1;
  std::vector<VectorXd>& alpha2;
  MatrixXd& beta;
  VectorXd& Y_prec;

  regression_sampler_worker(
    const int sampler,
    const VectorXi& trait_set,
    const Map<MatrixXd> Y,           // nx1
    const Map<MatrixXd> W_base,           // nxa
    const std::vector<R_matrix>& W_list,           // nxa
    const Map<MatrixXd> X_U,
    const Map<MatrixXd> V,
    const MatrixXd RinvSqX,                // nxb
    const MatrixXd& C,                     // bxb
    const R_matrix& chol_R,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
    const MatrixXd& R,
    const MatrixXd& Rinv,
    const MatrixXd& RinvU,
    const MatrixXd& UtRinvU,
    const Map<MatrixXd> prior_prec_alpha1, // a1 x p matrix of prior precisions for alpha1
    const VectorXd& prior_prec_alpha2,     // p-vector of precision of alpha2s for each trait
    const Map<MatrixXd> prior_mean_beta, // b x p matrix of prior means of beta
    const Map<MatrixXd> prior_prec_beta, // b x p matrix of prior precisions of beta
    double Y_prec_b0,
    const MatrixXd& randn_alpha1,
    const std::vector<VectorXd>& randn_alpha2,
    const MatrixXd& randn_beta,
    const MatrixXd& randn_e,
    const VectorXd& rgamma_1,
    MatrixXd& alpha1,
    std::vector<VectorXd>& alpha2,
    MatrixXd& beta,
    VectorXd& Y_prec
  ):
    sampler(sampler), trait_set(trait_set),
    Y(Y), W_base(W_base), W_list(W_list), X_U(X_U), V(V), RinvSqX(RinvSqX),
    C(C), chol_R(chol_R), R(R), Rinv(Rinv), RinvU(RinvU), UtRinvU(UtRinvU),
    prior_prec_alpha1(prior_prec_alpha1), prior_prec_alpha2(prior_prec_alpha2), prior_mean_beta(prior_mean_beta), prior_prec_beta(prior_prec_beta),
    Y_prec_b0(Y_prec_b0),
    randn_alpha1(randn_alpha1), randn_alpha2(randn_alpha2), randn_beta(randn_beta), randn_e(randn_e), rgamma_1(rgamma_1),
    alpha1(alpha1), alpha2(alpha2), beta(beta), Y_prec(Y_prec)
  {}

  void operator()(std::size_t begin, std::size_t end) {
    int n = Y.rows();
    int a1 = W_base.cols();

    for(std::size_t i = begin; i < end; i++){
      int j = trait_set[i];
      MatrixXd W;
      int a;
      int a2 = 0;
      int b;
      VectorXd prior_prec_alpha;
      VectorXd randn_alpha;
      if(W_list.size() == 0) {
        W = W_base;
        a = a1;
        prior_prec_alpha = prior_prec_alpha1.col(j);
        randn_alpha = randn_alpha1.col(j);
      } else{
        Map<MatrixXd> W2 = W_list[j].dense;
        a2 = W2.cols();
        a = a1+a2;
        W = MatrixXd(n,a);
        W << W_base,W2;
        prior_prec_alpha = VectorXd(a);
        prior_prec_alpha.head(a1) = prior_prec_alpha1.col(j);
        prior_prec_alpha.tail(a2).array() = prior_prec_alpha2[j];
        randn_alpha = VectorXd(a);
        randn_alpha.head(a1) = randn_alpha1.col(j);
        randn_alpha.tail(a2) = randn_alpha2[j];
      }

      VectorXd samples;
      if(sampler == 1) {
        b = RinvSqX.cols();
        samples = regression_sampler_v1(Y.col(j), W, RinvSqX, C, prior_prec_alpha, prior_mean_beta.col(j),
                                        prior_prec_beta.col(j), chol_R, Y_prec[j], randn_alpha,
                                        randn_beta.col(j), rgamma_1[j],Y_prec_b0);
      } else if(sampler == 2) {
        b = X_U.cols();
        samples = regression_sampler_v2(Y.col(j), W, X_U, prior_prec_alpha, prior_mean_beta.col(j),
                                        prior_prec_beta.col(j), chol_R, R, Y_prec[j], randn_alpha,
                                        randn_beta.col(j), randn_e.col(j),rgamma_1[j],Y_prec_b0);
      } else if(sampler == 3) {
        b = V.cols();
        samples = regression_sampler_v3(Y.col(j), W, X_U, V, prior_prec_alpha, prior_mean_beta.col(j),
                                        prior_prec_beta.col(j), chol_R, Rinv, RinvU, UtRinvU, Y_prec[j], randn_alpha,
                                        randn_beta.col(j), randn_e.col(j),rgamma_1[j],Y_prec_b0);

      } else {
        stop("sampler not implemented");
      }

      // extract samples
      Y_prec[j] = samples[0];
      if(a1 > 0) alpha1.col(j) = samples.segment(1,a1);
      if(a2 > 0) alpha2[j] = samples.segment(1+a1,a2);
      if(b > 0) beta.col(j) = samples.tail(b);
    }
  }
};


// [[Rcpp::export]]
Rcpp::List regression_sampler_parallel(
    Map<MatrixXd> Y,               // n x p matrix of observations
    Map<MatrixXd> W_base,          // n x a1 matrix of W covariates common to all p. Can be NULL
    Rcpp::List W_list_,             // p-list of n x a2 matrices of W covariates unique to each p. Can be NULL
    Map<MatrixXd> X,               // either X, a n x b matrix, or U, a n x m matrix. If U, then V must be non-NULL
    SEXP V_,                       // m x b matrix if X is U
    Rcpp::IntegerVector h2s_index, // p-vector of indices for appropriate V of each trait
    Rcpp::List chol_V_list_,        // list of cholesky decompositions of V: RtR (each nxn). Can be either dense or sparse
    VectorXd Y_prec,               // p-vector of Y current precisions
    double Y_prec_a0,
    double Y_prec_b0,
    Map<MatrixXd> prior_prec_alpha1, // a1 x p matrix of prior precisions for alpha1
    VectorXd prior_prec_alpha2,     // p-vector of precision of alpha2s for each trait
    Map<MatrixXd> prior_mean_beta, // b x p matrix of prior means of beta
    Map<MatrixXd> prior_prec_beta, // b x p matrix of prior precisions of beta
    int grainSize) {

  int n = Y.rows();
  int p = Y.cols();

  // W_base
  if(W_base.rows() != n) stop("Wrong dimension of W_base");
  int a1 = W_base.cols();

  // W_list
  std::vector<R_matrix> W_list;
  load_R_matrices_list(W_list_, W_list);
  if(W_list.size() > 0) {
    if(W_list.size() != p) stop("Wrong length of W_list");
  }

  // X or U and V
  Map<MatrixXd> U = X;
  MatrixXd z = MatrixXd::Zero(0,0);
  Map<MatrixXd> V(z.data(),0,0);
  int b = X.cols();
  if(X.rows() != n) stop("Wrong dimension of X");
  if(Rf_isMatrix(V_)) {
    // Map<MatrixXd> V__ = as<Map<MatrixXd> >(V_);
    // new (&v) Map<MatrixXd> (V__,V__.rows(),V__.cols());
    new (&V) Map<MatrixXd> (as<Map<MatrixXd> >(V_));
    if(U.cols() != V.rows()) stop("X and V_ have incompatible dimensions");
    b = V.cols();
  }

  // chol_V_list
  std::vector<R_matrix> chol_V_list;
  load_R_matrices_list(chol_V_list_, chol_V_list);
  if(max(h2s_index) > chol_V_list.size()) {
    stop("max(h2s_index) > length(chol_V_list)");
  }

  // priors
  if(Y_prec.size() != p) {
    stop("Wrong length of Y_prec");
  }
  if(prior_prec_alpha1.rows() != a1 || prior_prec_alpha1.cols() != p) stop("Wrong dimensions of prior_prec_alpha1");
  if(W_list.size() > 0 && prior_prec_alpha2.size() != p) {
    stop("Wrong length of prior_prec_alpha2");
  }
  if(prior_mean_beta.rows() != b || prior_mean_beta.cols() != p) stop("Wrong dimensions of prior_mean_beta");
  if(prior_prec_beta.rows() != b || prior_prec_beta.cols() != p) stop("Wrong dimensions of prior_prec_beta");

  // generate random numbers
  MatrixXd randn_alpha1 = rstdnorm_mat(a1,p);
  std::vector<VectorXd> randn_alpha2;
  if(W_list.size() > 0){
    for(int i = 0; i < p; i++){
      randn_alpha2.push_back(rstdnorm_mat(W_list[i].dense.cols(),1));
    }
  }
  MatrixXd randn_beta = rstdnorm_mat(b,p);
  MatrixXd randn_e;
  if(b > n) {
    randn_e = rstdnorm_mat(n,p);
  }
  VectorXd rgamma_1 = as<VectorXd>(rgamma(p,Y_prec_a0 + n/2.0,1.0));

  // Results structures
  MatrixXd alpha1(a1,p);
  std::vector<VectorXd> alpha2;
  alpha2.reserve(W_list.size());
  int alpha2_size = 0;
  if(W_list.size() > 0){
    for(int i = 0; i < W_list.size(); i++){
      int a2 = W_list[i].dense.cols();
      alpha2.push_back(VectorXd::Zero(a2));
      alpha2_size += a2;
    }
  }
  MatrixXd beta(b,p);

  // go through h2s indices and sample columns with same index as a set
  for(int i = min(h2s_index); i <= max(h2s_index); i++) {
    int h2_index = i;
    VectorXi trait_set = as<VectorXi>(which(h2s_index == h2_index));  // list of traits with same h2_index

    if(trait_set.size() > 0){
      // prepare matrices for sampler
      MatrixXd RinvSqX, C, R, Rinv, RinvU, UtRinvU;
      R_matrix chol_R = chol_V_list[h2_index - 1];
      int which_sampler;
      // Decide which sampler to use
      if(b <= n) {
        // use regression_sampler_v1
        which_sampler = 1;
        if(chol_R.isDense) {
          RinvSqX = chol_R.dense.transpose().triangularView<Lower>().solve(X);
        } else{
          RinvSqX = chol_R.sparse.transpose().triangularView<Lower>().solve(X);
        }
        C = RinvSqX.transpose() * RinvSqX;
      }
      else if(V.cols() == 0) {
        // use regression_sampler_v2
        which_sampler = 2;
        if(chol_R.isDense) {
          R = chol_R.dense.transpose().triangularView<Lower>() * chol_R.dense;
        } else{
          R = chol_R.sparse.transpose().triangularView<Lower>() * chol_R.sparse;
        }
      } else {
        // use regression_sampler_v3
        which_sampler = 3;
        if(chol_R.isDense) {
          Rinv = chol_R.dense.triangularView<Upper>().solve(chol_R.dense.transpose().triangularView<Lower>().solve(MatrixXd::Identity(n,n)));
          RinvU = chol_R.dense.triangularView<Upper>().solve(chol_R.dense.transpose().triangularView<Lower>().solve(U));
        } else{
          Rinv = chol_R.sparse.triangularView<Upper>().solve(chol_R.sparse.transpose().triangularView<Lower>().solve(MatrixXd::Identity(n,n)));
          RinvU = chol_R.sparse.triangularView<Upper>().solve(chol_R.sparse.transpose().triangularView<Lower>().solve(U));
          // RinvU = Rinv * U;
          // Rcout << i << std::endl;
          // Rcout << Rinv.diagonal().transpose() << std::endl;
        }
        UtRinvU = U.transpose() * RinvU;
      }
      regression_sampler_worker sampler(which_sampler, trait_set,Y,W_base,W_list,X,V,RinvSqX,C,
                                        chol_R,R,Rinv,RinvU,UtRinvU,
                                        prior_prec_alpha1,prior_prec_alpha2,prior_mean_beta,prior_prec_beta,Y_prec_b0,
                                        randn_alpha1,randn_alpha2,randn_beta,randn_e, rgamma_1,
                                        alpha1,alpha2,beta,Y_prec);
      RcppParallel::parallelFor(0,trait_set.size(),sampler,grainSize);
    }
  }

  // collect alpha2 into a vector
  VectorXd alpha2_vec(alpha2_size);
  if(W_list.size() > 0){
    int index = 0;
    for(int i = 0; i < W_list.size(); i++){
      alpha2_vec.segment(index,alpha2[i].size()) = alpha2[i];
      index += alpha2[i].size();
    }
  }

  return(Rcpp::List::create(
      Named("alpha1") = alpha1,
      Named("alpha2") = alpha2_vec,
      Named("beta") = beta,
      Named("Y_prec") = Y_prec
  ));
}


// -------------------------------------------- //
// ------------ sample_MME_ZKZts -------------- //
// -------------------------------------------- //

// Samples from model:
// y = Zu + e
// u ~ N(0,K); solve(K) = t(chol_K_inv) %*% chol_K_inv
// e[i] ~ N(0,1/tot_Eta_prec)
// C = ZtRinvZ + diag(Kinv)
// [[Rcpp::export()]]
VectorXd sample_MME_single_diagR(
    VectorXd y,           // nx1
    MSpMat Z,             // nxr dgCMatrix
    MSpMat chol_ZtZ_Kinv,       // rxr CsparseMatrix upper triangular: chol(ZtRinvZ + diag(Kinv))
    double tot_Eta_prec,   // double
    double pe,            // double
    VectorXd randn_theta  // rx1
){
  VectorXd b = Z.transpose() * y * pe;
  b = chol_ZtZ_Kinv.transpose().triangularView<Lower>().solve(b / sqrt(tot_Eta_prec));
  b += randn_theta;
  b = chol_ZtZ_Kinv.triangularView<Upper>().solve(b / sqrt(tot_Eta_prec));
  return(b);
}


struct sample_MME_single_diagR_worker : public RcppParallel::Worker {
  const Map<MatrixXd> Y;
  const MSpMat Z;
  const std::vector<R_matrix> chol_ZtZ_Kinv_list;
  const ArrayXd& pes;
  const Map<VectorXd> tot_Eta_prec;
  const VectorXi& h2s_index;
  const MatrixXd& randn_theta;
  MatrixXd &coefs;

  sample_MME_single_diagR_worker(
    const Map<MatrixXd> Y,                                // nxp
    const MSpMat Z,                                  // nxr
    const std::vector<R_matrix> chol_ZtZ_Kinv_list,     // std::vector of rxr SpMat upper-triangle
    const ArrayXd& pes,                               // px1
    const Map<VectorXd> tot_Eta_prec,                     // px1
    const VectorXi& h2s_index,                        // px1, 1-based index
    const MatrixXd& randn_theta,                      // rxp
    MatrixXd &coefs                            // rxp
  ):
    Y(Y), Z(Z),
    chol_ZtZ_Kinv_list(chol_ZtZ_Kinv_list),
    pes(pes), tot_Eta_prec(tot_Eta_prec),h2s_index(h2s_index),
    randn_theta(randn_theta),
    coefs(coefs) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++){
      int h2_index = h2s_index[j] - 1;
      // ZtZ_Kinv needs to be scaled by tot_Eta_prec[j].
      coefs.col(j) = sample_MME_single_diagR(Y.col(j), Z, chol_ZtZ_Kinv_list[h2_index].sparse, tot_Eta_prec[j], pes[j],randn_theta.col(j));
    }
  }
};


// samples random effects from model:
// Y = ZU + E
// U[,j] ~ N(0,1/tot_Eta_prec[j] * h2[j] * K)
// E[,j] ~ N(0,1/tot_Eta_prec[j] * (1-h2[j]) * I_n)
// For complete data, ie no missing obs.
// [[Rcpp::export()]]
MatrixXd sample_MME_ZKZts_c(
    Map<MatrixXd> Y,                    // nxp
    MSpMat Z,                           // nxr
    Map<VectorXd> tot_Eta_prec,         // px1
    Rcpp::List chol_ZtZ_Kinv_list_,      // List or R st RtR = ZtZ_Kinv
    Map<MatrixXd> h2s,                  // n_RE x p
    VectorXi h2s_index,                 // px1
    int grainSize) {

  int p = Y.cols();
  int r = Z.cols();

  MatrixXd randn_theta = rstdnorm_mat(r,p);

  std::vector<R_matrix> chol_ZtZ_Kinv_list;
  load_R_matrices_list(chol_ZtZ_Kinv_list_, chol_ZtZ_Kinv_list);

  MatrixXd U(r,p);
  ArrayXd h2_e = 1.0 - h2s.colwise().sum().array();
  ArrayXd pes = tot_Eta_prec.array() / h2_e.array();

  sample_MME_single_diagR_worker sampler(Y,Z,chol_ZtZ_Kinv_list,pes,tot_Eta_prec,h2s_index,randn_theta,U);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return(U);
}

// -------------------------------------------- //
// ---------------- sample h2s ---------------- //
// -------------------------------------------- //

struct log_ps_worker : public RcppParallel::Worker {
  const Map<MatrixXd> Y;
  const Map<VectorXd> tot_Eta_prec;
  const std::vector<R_matrix>& chol_V_list;
  const Map<VectorXd> discrete_priors;
  MatrixXd &log_ps;

  log_ps_worker(
    const Map<MatrixXd> Y,                           // nxp
    const Map<VectorXd> tot_Eta_prec,                 // px1
    const std::vector<R_matrix>& chol_V_list, // std::vector of nxn SpMat, upper-triangular
    const Map<VectorXd> discrete_priors,              // n_h2 x 1
    MatrixXd &log_ps                       // n_h2 x p
  ):
    Y(Y), tot_Eta_prec(tot_Eta_prec), chol_V_list(chol_V_list),
    discrete_priors(discrete_priors), log_ps(log_ps) {}

  void operator()(std::size_t begin, std::size_t end) {
    int n = Y.rows();
    for(std::size_t i = begin; i < end; i++){
      R_matrix chol_R = chol_V_list[i];
      MatrixXd y_std(n,end-begin);
      double log_det_V;
      if(chol_R.isDense){
        y_std = chol_R.dense.transpose().triangularView<Lower>().solve(Y);
        log_det_V = 2*chol_R.dense.diagonal().array().log().sum();
      } else{
        y_std = chol_R.sparse.transpose().triangularView<Lower>().solve(Y);
        log_det_V = 0;
        for(int j = 0; j < chol_R.sparse.rows(); j++) {
          log_det_V += 2*std::log(chol_R.sparse.coeffRef(j,j));
        }
      }
      VectorXd scores2 = (y_std.transpose() * y_std).diagonal().array() * tot_Eta_prec.array();
      log_ps.row(i) = (-n/2.0 * log(2*M_PI) - 0.5 * (log_det_V - n*tot_Eta_prec.array().log()) -
        0.5 * scores2.array() + log(discrete_priors[i]));
    }
  }
};

// [[Rcpp::export()]]
MatrixXd log_p_h2s(
    Map<MatrixXd> Y,              // nxp
    Map<VectorXd> tot_Eta_prec,   // px1
    Rcpp::List chol_V_list_,       // List. Each element contains: R st RtR = V. may be sparse or dense
    Map<VectorXd> discrete_priors,// n_h2 x 1
    int grainSize)
{
  int b = discrete_priors.size();
  int p = Y.cols();

  std::vector<R_matrix> chol_V_list;
  load_R_matrices_list(chol_V_list_, chol_V_list);

  MatrixXd log_ps(b,p);

  log_ps_worker sampler(Y,tot_Eta_prec,chol_V_list,discrete_priors,log_ps);
  RcppParallel::parallelFor(0,b,sampler,grainSize);
  return(log_ps);
}

// [[Rcpp::export()]]
VectorXi sample_h2s(
    Map<ArrayXXd> log_ps,
    int grainSize
)
{

  int p = log_ps.cols();
  VectorXd rs = as<VectorXd>(runif(p));

  VectorXi h2s_index(p);


  struct sampleColumn : public RcppParallel::Worker {
    ArrayXXd log_ps;
    VectorXd rs;
    VectorXi &h2s_index;

    sampleColumn(ArrayXXd log_ps,
                 VectorXd rs,
                 VectorXi &h2s_index):
      log_ps(log_ps), rs(rs), h2s_index(h2s_index) {}

    void operator()(std::size_t begin, std::size_t end) {
      int b = log_ps.rows();
      for(std::size_t j = begin; j < end; j++){
        // for(int j = 0; j < p; j++){
        double max_col = log_ps.col(j).maxCoeff();
        double norm_factor = max_col + log((log_ps.col(j) - max_col).exp().sum());
        VectorXd ps_j = (log_ps.col(j) - norm_factor).exp();
        h2s_index[j] = 1;
        double cumsum = 0;
        for(int i = 0; i < b; i++){
          cumsum += ps_j[i];
          if(rs[j] > cumsum) {
            h2s_index[j] ++;
          }
        }
      }
    }
  };
  sampleColumn sampler(log_ps,rs,h2s_index);
  RcppParallel::parallelFor(0,p,sampler,grainSize);

  return h2s_index; // 1-based index
}


// -------------------------------------------- //
// --------------- sample h2s MH -------------- //
// -------------------------------------------- //

double log_prob_h2_c(
    const Ref<const VectorXd>& y,           // nx1
    R_matrix chol_R,     // nxn upper-triangular. Dense or sparse
    int n,                // int
    double tot_Eta_prec,  // double
    double discrete_prior // double
){
  VectorXd y_std;
  double log_det_V;
  if(chol_R.isDense){
    y_std = chol_R.dense.transpose().triangularView<Lower>().solve(y);
    log_det_V = 2*chol_R.dense.diagonal().array().log().sum();
  } else{
    y_std = chol_R.sparse.transpose().triangularView<Lower>().solve(y);
    log_det_V = 0;
    for(int j = 0; j < chol_R.sparse.rows(); j++) {
      log_det_V += 2*std::log(chol_R.sparse.coeffRef(j,j));
    }
  }
  double score2 = tot_Eta_prec * y_std.dot(y_std);

  double log_p = -n/2.0 * log(2*M_PI) - 0.5*(log_det_V - n*log(tot_Eta_prec)) - 0.5 * score2 + log(discrete_prior);
  return log_p;
}

struct sample_h2s_discrete_MH_worker : public RcppParallel::Worker {
  const Map<MatrixXd> Y;
  const MatrixXd h2s_matrix;
  const std::vector<R_matrix>& chol_V_list;
  const VectorXd tot_Eta_prec;
  const VectorXd discrete_priors;
  const VectorXd r_draws;
  const VectorXd state_draws;
  const VectorXi h2_index;
  const double step_size;
  VectorXi &new_index;

  sample_h2s_discrete_MH_worker(const Map<MatrixXd> Y,                 // nxp
                                const MatrixXd h2s_matrix,        // n_RE x n_h2
                                const std::vector<R_matrix>& chol_V_list,     // List of R st RtR = V
                                const VectorXd tot_Eta_prec,      // p x 1
                                const VectorXd discrete_priors,   // n_h2 x 1
                                const VectorXd r_draws,           // px1
                                const VectorXd state_draws,       // px1
                                const VectorXi h2_index,          // px1 - 1-based index
                                const double step_size,           // double
                                VectorXi &new_index               // px1 - 1-based index
  ):

    Y(Y), h2s_matrix(h2s_matrix),
    chol_V_list(chol_V_list),
    tot_Eta_prec(tot_Eta_prec),  discrete_priors(discrete_priors),
    r_draws(r_draws),state_draws(state_draws),h2_index(h2_index),step_size(step_size),new_index(new_index) {}

  void operator()(std::size_t begin, std::size_t end) {
    int n = Y.rows();
    for(std::size_t j = begin; j < end; j++){
      int old_state = h2_index[j] - 1;
      VectorXd candidate_new_states = find_candidate_states(h2s_matrix,step_size,old_state);
      int r = state_draws[j] * (candidate_new_states.size());
      int proposed_state = candidate_new_states[r];

      if(discrete_priors[proposed_state] == 0.0) {
        new_index[j] = old_state;  // don't bother with calculations if prior == 0.0
      } else{
        double old_log_p = log_prob_h2_c(Y.col(j),chol_V_list[old_state],n,tot_Eta_prec[j],discrete_priors[old_state]);
        double new_log_p = log_prob_h2_c(Y.col(j),chol_V_list[proposed_state],n,tot_Eta_prec[j],discrete_priors[proposed_state]);

        VectorXd candidate_states_from_new_state = find_candidate_states(h2s_matrix,step_size,proposed_state);

        double forward_prob = 1.0 / candidate_new_states.size();
        double back_prob = 1.0 / candidate_states_from_new_state.size();

        double log_MH_ratio = new_log_p - old_log_p + log(forward_prob) - log(back_prob);

        if(log(r_draws[j]) < log_MH_ratio) {
          new_index[j] = proposed_state;
        } else {
          new_index[j] = old_state;
        }
      }
      new_index[j] += 1;  // convert to 1-based index
    }
  }
};

// [[Rcpp::export()]]
VectorXi sample_h2s_discrete_MH_c(
    Map<MatrixXd> Y,                // nxp
    Map<VectorXd> tot_Eta_prec,     // px1
    Map<VectorXd> discrete_priors,  // n_h2 x 1
    VectorXi h2s_index,              // px1
    Map<MatrixXd> h2s_matrix,       // n_RE x n_h2
    Rcpp::List chol_V_list_,         // List of R st RtR = V, can be dense or sparse
    double step_size,               // double
    int grainSize
){

  int p = Y.cols();

  VectorXd r_draws = as<VectorXd>(runif(p));
  VectorXd state_draws = as<VectorXd>(runif(p));

  VectorXi new_index(p);

  std::vector<R_matrix> chol_V_list;
  load_R_matrices_list(chol_V_list_, chol_V_list);

  sample_h2s_discrete_MH_worker sampler(Y,h2s_matrix,chol_V_list,tot_Eta_prec,
                                        discrete_priors,r_draws,state_draws,h2s_index,step_size,new_index);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return new_index;
}


// -------------------------------------------------- //
// -- Sample factor scores --- //
// -------------------------------------------------- //

// Sample factor scores given factor loadings (F_a), factor residual variances (F_e_prec) and
// phenotype residuals
// Y - ZU = F * Lambda + E
// F[,k] ~ N(XFBF[,k] + ZUF[,k],1/F_e_prec[,j])
// E[,j] ~ N(0,1/resid_Eta_prec[,j])
// Sampling is done separately for each block of rows with the same pattern of missing observations
// [[Rcpp::export()]]
MatrixXd sample_factors_scores_c( // returns nxk matrix
    Map<MatrixXd> Eta_tilde,      // nxp
    Map<MatrixXd> prior_mean,     // nxk
    Map<MatrixXd> Lambda,         // nxp
    Map<VectorXd> resid_Eta_prec, // px1
    Map<VectorXd> F_e_prec        // kx1
) {
  int n = Eta_tilde.rows();
  int k = Lambda.cols();
  MatrixXd randn_draws = rstdnorm_mat(n,k);

  MatrixXd Lmsg = resid_Eta_prec.asDiagonal() * Lambda;
  MatrixXd Sigma = Lambda.transpose() * Lmsg;
  Sigma.diagonal() += F_e_prec;
  Eigen::LLT<MatrixXd> chol_Sigma;
  chol_Sigma.compute(Sigma);
  MatrixXd R = chol_Sigma.matrixU();

  MatrixXd Meta = R.transpose().triangularView<Lower>().solve((Eta_tilde * Lmsg + prior_mean * F_e_prec.asDiagonal()).transpose());

  MatrixXd Ft = R.triangularView<Upper>().solve(Meta + randn_draws.transpose());

  return Ft.transpose();
}



// -------------------------------------------------- //
// -- Sample tau2 and delta scores --- //
// -------------------------------------------------- //

VectorXd cumprod(const VectorXd& x) {
  int n = x.size();
  VectorXd res(n);
  res[0] = x[0];
  if(n > 1) {
    for(int i = 1; i < n; i++){
      res[i] = res[i-1]*x[i];
    }
  }
  return(res);
}


// [[Rcpp::export()]]
Rcpp::List sample_tau2_delta_c_Eigen_v2(
    double tau2,
    double xi,
    VectorXd delta,
    Map<VectorXd> scores,
    double tau_0,
    double delta_shape,
    double delta_rate,
    int p,
    int times
) {

  int K = scores.size();
  if(delta.size() != K) stop("Wrong size of delta");
  double shape;
  double rate;
  VectorXd cumprod_delta = cumprod(delta);
  for(int i = 0; i < times; i++){
    // sample tau2
    shape = (p*K + 1)/2.0;
    rate = 1.0/xi + cumprod_delta.dot(scores);
    tau2 = 1.0/R::rgamma(shape,1.0/rate);

    // sample xi
    shape = 1.0;
    rate = 1.0/(tau_0*tau_0) + 1.0/tau2;
    xi = 1.0/R::rgamma(shape,1.0/rate);

    for(int h = 1; h < K; h++) {
      // delta_h
      shape = delta_shape + p*(K-h)/2.0;
      rate = delta_rate + cumprod_delta.tail(K-h).dot(scores.tail(K-h)) / (tau2 * delta(h));
      delta[h] = R::rgamma(shape,1.0/rate);
      cumprod_delta = cumprod(delta);
    }
  }

  return(Rcpp::List::create(Named("tau2") = tau2,
                            Named("xi") = xi,
                            Named("delta") = delta));
}



