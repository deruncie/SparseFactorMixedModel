
typedef Eigen::SparseMatrix<float> SpMatf;


// [[Rcpp::export()]]
MatrixXf rstdnorm_matb(int n,int p) {  // returns nxp matrix
  VectorXf X_vec(n*p);
  for(int i = 0; i < n*p; i++){
    X_vec[i] = ziggr.norm();
  }
  MatrixXf X_mat = MatrixXf(X_vec.data(),n,p);
  return(X_mat);
}

struct R_matrixb {
  MatrixXf dense;
  SpMatf sparse;
  bool isDense;
  R_matrixb(MatrixXf dense_, SpMatf sparse_,bool isDense_) : dense(dense_), sparse(sparse_), isDense(isDense_) {}
};
// void load_R_matrices_list(const Rcpp::List X_list, std::vector<R_matrixb>& X_vector);


void load_R_matrices_listb(const Rcpp::List X_list, std::vector<R_matrixb>& X_vector){
  // null_matrices
  MatrixXf null_d = MatrixXf::Zero(0,0);
  MatrixXf M_null_d(null_d.data(),0,0);
  SpMatf null_s = null_d.sparseView();
  SpMatf M_null_s(0,0,0,null_s.outerIndexPtr(),null_s.innerIndexPtr(),null_s.valuePtr());
  
  int p = X_list.size();
  X_vector.reserve(p);
  for(int i = 0; i < p; i++){
    SEXP Xi_ = X_list[i];
    if(Rf_isMatrix(Xi_)){
      MatrixXf Xi = as<MatrixXf >(Xi_);
      R_matrixb Xim(Xi,M_null_s,true);
      X_vector.push_back(Xim);
    } else{
      SpMatf Xi = as<SpMatf>(Xi_);
      R_matrixb Xim(M_null_d,Xi,false);
      X_vector.push_back(Xim);
    }
  }
}

MatrixXf get_RinvSqXb(const R_matrixb& chol_R, MatrixXf X){
  MatrixXf RinvSqX;
  if(chol_R.isDense) {
    RinvSqX = chol_R.dense.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
  } else{
    RinvSqX = chol_R.sparse.transpose().triangularView<Lower>().solve(X);// * sqrt(tot_Eta_prec));
  }
  return(RinvSqX);
}



VectorXf regression_sampler_v1b(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b < n
    const Ref<const VectorXf>& y,           // nx1
    const MatrixXf& W,           // nxa
    const MatrixXf& RinvSqX,                // nxb
    const MatrixXf& C,                     // bxb
    const VectorXf& prior_prec_alpha, // ax 1
    const Ref<const VectorXf>& prior_mean_beta,  // bx1
    const Ref<const VectorXf>& prior_prec_beta,  // bx1
    const R_matrixb& chol_R,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
    double Y_prec,                    // double
    const VectorXf& randn_alpha,
    const Ref<const VectorXf>& randn_beta,
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
  MatrixXf C_beta = C;
  C_beta.diagonal() += prior_prec_beta;
  LLT<MatrixXf> A_beta_llt;
  A_beta_llt.compute(C_beta);
  MatrixXf chol_A_beta = A_beta_llt.matrixU();
  // Y_prec * chol_A_beta^\T * chol_A_beta = A_beta
  
  // Step 1
  VectorXf alpha(a);
  VectorXf y_tilde = y;
  if(a > 0) {
    // Sample alpha
    // Calculate A_alpha = Y_prec*W^T*Sigma_beta^{-1}*W + D_alpha^{-1}
    // We don't need to actually calculate Sigma_beta^{-1} directly.
    MatrixXf RinvSqW = get_RinvSqXb(chol_R,W);  // n*n*a -> n x a
    MatrixXf WtRinvX = RinvSqW.transpose() * RinvSqX; // a*n*b -> a*b
    MatrixXf invSqAbXtRinvW = chol_A_beta.transpose().triangularView<Lower>().solve(WtRinvX.transpose()); // b*b*a -> b x a
    
    MatrixXf A_alpha = Y_prec * (RinvSqW.transpose() * RinvSqW - invSqAbXtRinvW.transpose() * invSqAbXtRinvW);
    A_alpha.diagonal() += prior_prec_alpha;
    
    VectorXf Rinvsqy = get_RinvSqXb(chol_R,y); // n*n*q -> n x 1;
    VectorXf XtRinvy = RinvSqX.transpose() * Rinvsqy; // b*n*1 >- b x 1
    VectorXf invSqAbXtRinvy = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy); // b*b*1 -> b*1
    
    VectorXf WtSbinvy = RinvSqW.transpose() * Rinvsqy - invSqAbXtRinvW.transpose() * invSqAbXtRinvy;
    
    LLT<MatrixXf> A_alpha_llt;
    A_alpha_llt.compute(A_alpha);
    MatrixXf chol_A_alpha = A_alpha_llt.matrixU();
    
    alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(WtSbinvy) * Y_prec + randn_alpha;
    alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
    
    y_tilde = y - W * alpha;
  }
  
  // Step 2 - sample Y_prec
  // We don't need to actually calculate Sigma_beta^{-1} directly.
  VectorXf RinvSqy = get_RinvSqXb(chol_R,y_tilde);
  VectorXf XtRinvy = RinvSqX.transpose() * RinvSqy;
  VectorXf prod1 = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy);
  double score = Y_prec_b0 + (RinvSqy.dot(RinvSqy) - prod1.dot(prod1))/2;
  Y_prec = rgamma_1/score;
  
  // Step 3 - sample beta
  VectorXf XtRinvy_std_mu = XtRinvy*Y_prec + prior_prec_beta.asDiagonal()*prior_mean_beta;
  VectorXf beta = chol_A_beta.transpose().triangularView<Lower>().solve(XtRinvy_std_mu) / sqrt(Y_prec) + randn_beta;
  beta = chol_A_beta.triangularView<Upper>().solve(beta) / sqrt(Y_prec);
  
  VectorXf result(1+a+b);
  result << Y_prec,alpha,beta;
  
  return(result);
}


VectorXf regression_sampler_v2b(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b > n
    const Ref<const VectorXf>& y,           // nx1
    const MatrixXf& W,           // nxa
    const MatrixXf& X,           // nxm or nxb
    const VectorXf& prior_prec_alpha, // ax 1
    const Ref<const VectorXf>& prior_mean_beta,  // bx1
    const Ref<const VectorXf>& prior_prec_beta,  // bx1
    const R_matrixb& chol_R,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
    const MatrixXf& R,
    double Y_prec,                    // double
    const VectorXf& randn_alpha,
    const Ref<const VectorXf>& randn_beta,
    const Ref<const VectorXf>& randn_e,
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
  MatrixXf DXt = prior_prec_beta.cwiseInverse().asDiagonal() * X.transpose();
  MatrixXf Sigma_beta = X * DXt + R;
  LDLT<MatrixXf> Sigma_beta_ldlt;
  Sigma_beta_ldlt.compute(Sigma_beta);
  
  // Step 1
  VectorXf alpha(a);
  VectorXf y_tilde = y;
  if(a > 0) {
    // Sample alpha
    MatrixXf SbinvW = Sigma_beta_ldlt.solve(W);
    MatrixXf A_alpha = Y_prec * SbinvW.transpose() * W;
    A_alpha.diagonal() += prior_prec_alpha;
    
    LLT<MatrixXf> A_alpha_llt;
    A_alpha_llt.compute(A_alpha);
    MatrixXf chol_A_alpha = A_alpha_llt.matrixU();
    
    VectorXf WtSbinvy = SbinvW.transpose() * y;
    alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(WtSbinvy) * Y_prec + randn_alpha;
    alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
    y_tilde = y - W * alpha;
  }
  
  // Step 2 - sample Y_prec
  VectorXf e2 = y_tilde.transpose() * Sigma_beta_ldlt.solve(y_tilde);
  double score = Y_prec_b0 + e2[0]/2;
  Y_prec = rgamma_1/score;
  
  // Step 3 - sample beta
  // what about prior mean?
  VectorXf u = randn_beta.array() / (prior_prec_beta * Y_prec).cwiseSqrt().array() + prior_mean_beta.array();
  VectorXf v = sqrt(Y_prec) * X * u;
  if(chol_R.isDense) {
    v += chol_R.dense.transpose().triangularView<Lower>() * randn_e;
  } else{
    v += chol_R.sparse.transpose().triangularView<Lower>() * randn_e;
  }
  VectorXf w = Sigma_beta_ldlt.solve(y_tilde * sqrt(Y_prec) - v);
  VectorXf beta = u + DXt * w / sqrt(Y_prec);
  
  VectorXf result(1+a+b);
  result << Y_prec,alpha,beta;
  
  return(result);
}


VectorXf regression_sampler_v3b(  // returns vector of length 1 + a + b for y_prec, alpha, beta, useful when b > n > m
    const Ref<const VectorXf>& y,           // nx1
    const MatrixXf& W,           // nxa
    const MatrixXf& U,           // nxm or nxb
    const MatrixXf& V,           // mxb
    const VectorXf& prior_prec_alpha, // ax 1
    const Ref<const VectorXf>& prior_mean_beta,  // bx1
    const Ref<const VectorXf>& prior_prec_beta,  // bx1
    const R_matrixb& chol_R,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
    const MatrixXf& Rinv,
    const MatrixXf& RinvU,
    const MatrixXf& UtRinvU,
    double Y_prec,                    // double
    const VectorXf& randn_alpha,
    const Ref<const VectorXf>& randn_beta,
    const Ref<const VectorXf>& randn_e,
    const double rgamma_1,
    const double Y_prec_b0
){
  int n = y.size();
  int a = W.cols();
  if(V.rows() != U.cols()) stop("Wrong dimensions of V");
  MatrixXf X = U*V;
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
  // MatrixXf Sigma_beta_inv;
  // Using Ainv - Ainv * U * (I + BVAinvU)inv * BVAinv in case B = VDVt is singular
  MatrixXf DVt = prior_prec_beta.cwiseInverse().asDiagonal() * V.transpose();
  MatrixXf VDVt = V * DVt;
  if(RinvU.rows() != n) stop("Wrong dimensions of RinvU");
  MatrixXf inner = VDVt * UtRinvU;
  inner.diagonal().array() += 1.0;
  LDLT<MatrixXf> inner_ldlt;
  inner_ldlt.compute(inner);
  // Sigma_beta_inv = Rinv - RinvU * inner.ldlt().solve(VDVt * RinvU.transpose());  // Don't actually calculate this. Stay in mxm space
  
  // Step 1
  VectorXf alpha(a);
  VectorXf y_tilde = y;
  if(a > 0) {
    // Sample alpha
    // MatrixXf SbinvW = Sigma_beta_inv * W;
    // MatrixXf A_alpha = Y_prec * SbinvW.transpose() * W;
    MatrixXf RinvW = Rinv * W;
    MatrixXf UtRinvW = U.transpose() * RinvW;
    MatrixXf A_alpha = Y_prec * (W.transpose() * RinvW - UtRinvW.transpose() * inner_ldlt.solve(VDVt * UtRinvW));
    A_alpha.diagonal() += prior_prec_alpha;
    
    // VectorXf WtSbinvy = SbinvW.transpose() * y;
    VectorXf UtRinvy = RinvU.transpose() * y;
    VectorXf WtSbinvy = RinvW.transpose() * y - UtRinvW.transpose() * inner_ldlt.solve(VDVt * UtRinvy);
    
    LLT<MatrixXf> A_alpha_llt;
    A_alpha_llt.compute(A_alpha);
    MatrixXf chol_A_alpha = A_alpha_llt.matrixU();
    
    alpha = chol_A_alpha.transpose().triangularView<Lower>().solve(WtSbinvy) * Y_prec + randn_alpha;
    alpha = chol_A_alpha.triangularView<Upper>().solve(alpha);
    
    y_tilde = y - W * alpha;
  }
  
  // Step 2 - sample Y_prec
  // VectorXf e2 = y_tilde.transpose() * Sigma_beta_inv * y_tilde;
  VectorXf Rinv_y = Rinv * y_tilde;
  VectorXf UtRinvy = RinvU.transpose() * y_tilde;
  VectorXf e2 = y_tilde.transpose() * Rinv_y - UtRinvy.transpose() * inner_ldlt.solve(VDVt * UtRinvy);
  
  double score = Y_prec_b0 + e2[0]/2;
  Y_prec = rgamma_1/score;
  
  // Step 3 - sample beta
  // what about prior mean?
  VectorXf u = randn_beta.array() / (prior_prec_beta * Y_prec).cwiseSqrt().array() + prior_mean_beta.array();
  VectorXf v = sqrt(Y_prec) * X * u;
  if(chol_R.isDense) {
    v += chol_R.dense.transpose().triangularView<Lower>() * randn_e;
  } else{
    v += chol_R.sparse.transpose().triangularView<Lower>() * randn_e;
  }
  // VectorXf w = Sigma_beta_inv * (y_tilde * sqrt(Y_prec) - v);
  VectorXf e = y_tilde * sqrt(Y_prec) - v;
  VectorXf UtRinve = RinvU.transpose() * e;
  VectorXf w = Rinv * e - RinvU * inner_ldlt.solve(VDVt * UtRinve);
  
  VectorXf beta = u + DVt * (U.transpose() * w) / sqrt(Y_prec); //b*b*1 + b*n*1
  
  VectorXf result(1+a+b);
  result << Y_prec,alpha,beta;
  
  return(result);
}



struct regression_sampler_workerb : public RcppParallel::Worker {
  const int sampler;  // whicha sampler to use?
  const VectorXi& trait_set;
  const MatrixXf Y;           // nx1
  const MatrixXf W_base;           // nxa
  const std::vector<R_matrixb>& W_list;           // nxa
  const MatrixXf X_U;   // could be X or U.
  const MatrixXf V;
  const MatrixXf RinvSqX;                // nxb
  const MatrixXf& C;                     // bxb
  const R_matrixb& chol_R;                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
  const MatrixXf& R;
  const MatrixXf& Rinv;
  const MatrixXf& RinvU;
  const MatrixXf& UtRinvU;
  const MatrixXf prior_prec_alpha1; // a1 x p matrix of prior precisions for alpha1
  const VectorXf& prior_prec_alpha2;     // p-vector of precision of alpha2s for each trait
  const MatrixXf prior_mean_beta; // b x p matrix of prior means of beta
  const MatrixXf prior_prec_beta; // b x p matrix of prior precisions of beta
  double Y_prec_b0;
  
  const MatrixXf& randn_alpha1;
  const std::vector<VectorXf>& randn_alpha2;
  const MatrixXf& randn_beta;
  const MatrixXf& randn_e;
  const VectorXf& rgamma_1;
  
  MatrixXf& alpha1;
  std::vector<VectorXf>& alpha2;
  MatrixXf& beta;
  VectorXf& Y_prec;
  
  regression_sampler_workerb(
    const int sampler,
    const VectorXi& trait_set,
    const MatrixXf Y,           // nx1
    const MatrixXf W_base,           // nxa
    const std::vector<R_matrixb>& W_list,           // nxa
    const MatrixXf X_U,
    const MatrixXf V,
    const MatrixXf RinvSqX,                // nxb
    const MatrixXf& C,                     // bxb
    const R_matrixb& chol_R,                    // either a upper-triangular matrix or upper-triangular CsparseMatrix
    const MatrixXf& R,
    const MatrixXf& Rinv,
    const MatrixXf& RinvU,
    const MatrixXf& UtRinvU,
    const MatrixXf prior_prec_alpha1, // a1 x p matrix of prior precisions for alpha1
    const VectorXf& prior_prec_alpha2,     // p-vector of precision of alpha2s for each trait
    const MatrixXf prior_mean_beta, // b x p matrix of prior means of beta
    const MatrixXf prior_prec_beta, // b x p matrix of prior precisions of beta
    double Y_prec_b0,
    const MatrixXf& randn_alpha1,
    const std::vector<VectorXf>& randn_alpha2,
    const MatrixXf& randn_beta,
    const MatrixXf& randn_e,
    const VectorXf& rgamma_1,
    MatrixXf& alpha1,
    std::vector<VectorXf>& alpha2,
    MatrixXf& beta,
    VectorXf& Y_prec
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
      MatrixXf W;
      int a;
      int a2 = 0;
      int b;
      VectorXf prior_prec_alpha;
      VectorXf randn_alpha;
      if(W_list.size() == 0) {
        W = W_base;
        a = a1;
        prior_prec_alpha = prior_prec_alpha1.col(j);
        randn_alpha = randn_alpha1.col(j);
      } else{
        MatrixXf W2 = W_list[j].dense;
        a2 = W2.cols();
        a = a1+a2;
        W = MatrixXf(n,a);
        W << W_base,W2;
        prior_prec_alpha = VectorXf(a);
        prior_prec_alpha.head(a1) = prior_prec_alpha1.col(j);
        prior_prec_alpha.tail(a2).array() = prior_prec_alpha2[j];
        randn_alpha = VectorXf(a);
        randn_alpha.head(a1) = randn_alpha1.col(j);
        randn_alpha.tail(a2) = randn_alpha2[j];
      }
      
      VectorXf samples;
      if(sampler == 1) {
        b = RinvSqX.cols();
        samples = regression_sampler_v1b(Y.col(j), W, RinvSqX, C, prior_prec_alpha, prior_mean_beta.col(j),
                                         prior_prec_beta.col(j), chol_R, Y_prec[j], randn_alpha,
                                         randn_beta.col(j), rgamma_1[j],Y_prec_b0);
      } else if(sampler == 2) {
        b = X_U.cols();
        samples = regression_sampler_v2b(Y.col(j), W, X_U, prior_prec_alpha, prior_mean_beta.col(j),
                                         prior_prec_beta.col(j), chol_R, R, Y_prec[j], randn_alpha,
                                         randn_beta.col(j), randn_e.col(j),rgamma_1[j],Y_prec_b0);
      } else if(sampler == 3) {
        b = V.cols();
        samples = regression_sampler_v3b(Y.col(j), W, X_U, V, prior_prec_alpha, prior_mean_beta.col(j),
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
Rcpp::List regression_sampler_parallelb(
    MatrixXf Y,               // n x p matrix of observations
    MatrixXf W_base,          // n x a1 matrix of W covariates common to all p. Can be NULL
    Rcpp::List W_list_,             // p-list of n x a2 matrices of W covariates unique to each p. Can be NULL
    MatrixXf X,               // either X, a n x b matrix, or U, a n x m matrix. If U, then V must be non-NULL
    SEXP V_,                       // m x b matrix if X is U
    Rcpp::IntegerVector h2s_index, // p-vector of indices for appropriate V of each trait
    Rcpp::List chol_V_list_,        // list of cholesky decompositions of V: RtR (each nxn). Can be either dense or sparse
    VectorXf Y_prec,               // p-vector of Y current precisions
    double Y_prec_a0,
    double Y_prec_b0,
    MatrixXf prior_prec_alpha1, // a1 x p matrix of prior precisions for alpha1
    VectorXf prior_prec_alpha2,     // p-vector of precision of alpha2s for each trait
    MatrixXf prior_mean_beta, // b x p matrix of prior means of beta
    MatrixXf prior_prec_beta, // b x p matrix of prior precisions of beta
    int grainSize) {
  
  int n = Y.rows();
  int p = Y.cols();
  
  // W_base
  if(W_base.rows() != n) stop("Wrong dimension of W_base");
  int a1 = W_base.cols();
  
  // W_list
  std::vector<R_matrixb> W_list;
  load_R_matrices_listb(W_list_, W_list);
  if(W_list.size() > 0) {
    if(W_list.size() != p) stop("Wrong length of W_list");
  }
  
  // X or U and V
  MatrixXf U = X;
  MatrixXf z = MatrixXf::Zero(0,0);
  MatrixXf V(z.data(),0,0);
  int b = X.cols();
  if(X.rows() != n) stop("Wrong dimension of X");
  if(Rf_isMatrix(V_)) {
    // MatrixXf V__ = as<MatrixXf >(V_);
    // new (&v) MatrixXf (V__,V__.rows(),V__.cols());
    new (&V) MatrixXf (as<MatrixXf >(V_));
    if(U.cols() != V.rows()) stop("X and V_ have incompatible dimensions");
    b = V.cols();
  }
  
  // chol_V_list
  std::vector<R_matrixb> chol_V_list;
  load_R_matrices_listb(chol_V_list_, chol_V_list);
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
  MatrixXf randn_alpha1 = rstdnorm_matb(a1,p);
  std::vector<VectorXf> randn_alpha2;
  if(W_list.size() > 0){
    for(int i = 0; i < p; i++){
      randn_alpha2.push_back(rstdnorm_matb(W_list[i].dense.cols(),1));
    }
  }
  MatrixXf randn_beta = rstdnorm_matb(b,p);
  MatrixXf randn_e;
  if(b > n) {
    randn_e = rstdnorm_matb(n,p);
  }
  VectorXf rgamma_1 = as<VectorXf>(rgamma(p,Y_prec_a0 + n/2.0,1.0));
  
  // Results structures
  MatrixXf alpha1(a1,p);
  std::vector<VectorXf> alpha2;
  alpha2.reserve(W_list.size());
  int alpha2_size = 0;
  if(W_list.size() > 0){
    for(int i = 0; i < W_list.size(); i++){
      int a2 = W_list[i].dense.cols();
      alpha2.push_back(VectorXf::Zero(a2));
      alpha2_size += a2;
    }
  }
  MatrixXf beta(b,p);
  
  // go through h2s indices and sample columns with same index as a set
  for(int i = min(h2s_index); i <= max(h2s_index); i++) {
    int h2_index = i;
    VectorXi trait_set = as<VectorXi>(whicha(h2s_index == h2_index));  // list of traits with same h2_index
    
    if(trait_set.size() > 0){
      // prepare matrices for sampler
      MatrixXf RinvSqX, C, R, Rinv, RinvU, UtRinvU;
      R_matrixb chol_R = chol_V_list[h2_index - 1];
      int which_sampler;
      // Decide whicha sampler to use
      if(b <= n) {
        // use regression_sampler_v1b
        which_sampler = 1;
        if(chol_R.isDense) {
          RinvSqX = chol_R.dense.transpose().triangularView<Lower>().solve(X);
        } else{
          RinvSqX = chol_R.sparse.transpose().triangularView<Lower>().solve(X);
        }
        C = RinvSqX.transpose() * RinvSqX;
      }
      else if(V.cols() == 0) {
        // use regression_sampler_v2b
        which_sampler = 2;
        if(chol_R.isDense) {
          R = chol_R.dense.transpose().triangularView<Lower>() * chol_R.dense;
        } else{
          R = chol_R.sparse.transpose().triangularView<Lower>() * chol_R.sparse;
        }
      } else {
        // use regression_sampler_v3b
        which_sampler = 3;
        if(chol_R.isDense) {
          Rinv = chol_R.dense.triangularView<Upper>().solve(chol_R.dense.transpose().triangularView<Lower>().solve(MatrixXf::Identity(n,n)));
          RinvU = chol_R.dense.triangularView<Upper>().solve(chol_R.dense.transpose().triangularView<Lower>().solve(U));
        } else{
          Rinv = chol_R.sparse.triangularView<Upper>().solve(chol_R.sparse.transpose().triangularView<Lower>().solve(MatrixXf::Identity(n,n)));
          RinvU = chol_R.sparse.triangularView<Upper>().solve(chol_R.sparse.transpose().triangularView<Lower>().solve(U));
          // RinvU = Rinv * U;
          // Rcout << i << std::endl;
          // Rcout << Rinv.diagonal().transpose() << std::endl;
        }
        UtRinvU = U.transpose() * RinvU;
      }
      regression_sampler_workerb sampler(which_sampler, trait_set,Y,W_base,W_list,X,V,RinvSqX,C,
                                         chol_R,R,Rinv,RinvU,UtRinvU,
                                         prior_prec_alpha1,prior_prec_alpha2,prior_mean_beta,prior_prec_beta,Y_prec_b0,
                                         randn_alpha1,randn_alpha2,randn_beta,randn_e, rgamma_1,
                                         alpha1,alpha2,beta,Y_prec);
      RcppParallel::parallelFor(0,trait_set.size(),sampler,grainSize);
    }
  }
  
  // collect alpha2 into a vector
  VectorXf alpha2_vec(alpha2_size);
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


