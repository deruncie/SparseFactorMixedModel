// draws a sample of the vector b from the model:
// y = X %*% beta + e
// with beta[i] ~ N(prior_mean[i],1/prior_prec[i]), i=1:b
// with e ~ N(0,t(chol_R) %*% chol(R)), i=1:n
// Uses sampling method from MCMCglmm, which requires draws b+n draws from N(0,1), which are passed as randn_theta and randn_e
// If b >= n, inverts the C matrix using Binomial Inverse Theorem
// [[Rcpp::export()]]
VectorXd sample_MME_single_diagK(  // returns b x 1 vector
    VectorXd y,           // nx1
    MatrixXd X,           // nxb
    VectorXd prior_mean,  // bx1
    VectorXd prior_prec,  // bx1
    SpMat chol_R,        // nxn upper triangular Cholesky decomposition of R. Format: dgCMatrix
    VectorXd randn_theta, // bx1
    VectorXd randn_e      // nx1
){
  int n = y.size();
  int b = X.cols();

  VectorXd theta_star = randn_theta.array() / prior_prec.cwiseSqrt().array();
  theta_star += prior_mean;
  VectorXd e_star = chol_R * randn_e;
  MatrixXd X_theta_star = X * theta_star;
  VectorXd y_resid = y - X_theta_star - e_star;

  MatrixXd RinvSqX = chol_R.transpose().triangularView<Lower>().solve(X);
  VectorXd XtRinvy = RinvSqX.transpose() * chol_R.transpose().triangularView<Lower>().solve(y_resid);

  VectorXd theta_tilda;

  if(b < n) {
    MatrixXd C = RinvSqX.transpose() * RinvSqX;
    C.diagonal() += prior_prec;
    theta_tilda = C.llt().solve(XtRinvy);
  } else{
    MatrixXd VAi = X * prior_prec.cwiseInverse().asDiagonal();
    MatrixXd inner = VAi*X.transpose() + chol_R.transpose() * chol_R;
    VectorXd VAiXtURinvy = VAi * XtRinvy;
    VectorXd outerXtURinvy = VAi.transpose() * inner.ldlt().solve(VAiXtURinvy);
    theta_tilda = XtRinvy.array() / prior_prec.array();
    theta_tilda -= outerXtURinvy;
  }

  VectorXd theta = theta_star + theta_tilda;

  return theta;
}

// RcppParallel struct for sampling coefficients from a set of parallel regression problems:
// Y = X %*% B + E, where
// b[i,j] ~ N(prior_mean[i,j],1/prior_prec[i,j])
// E[,j] ~ N(0,t(chol_R) %*% chol(R) / tot_Eta_prec[j])
// where chol_R is selected from a list based on h2s_index[j]
// Y is complete, so everythign has the same dimensions
struct sample_MME_single_diagK_worker : public RcppParallel::Worker {
  MatrixXd Y;
  MatrixXd X;
  MatrixXd prior_mean, prior_prec, randn_theta, randn_e;
  const std::vector<MSpMat> chol_R_list;
  VectorXi h2s_index;
  VectorXd tot_Eta_prec;
  MatrixXd &coefs;

  sample_MME_single_diagK_worker(
    MatrixXd Y,           // nxp
    MatrixXd X,           // nxb
    MatrixXd prior_mean,  // bxp
    MatrixXd prior_prec,  // bxp
    const std::vector<MSpMat> chol_R_list, // each element a MSpMat nxn upper-triangular cholesky decomposition of R
    VectorXi h2s_index,   // px1, 1-based index
    VectorXd tot_Eta_prec,// px1
    MatrixXd randn_theta, // bxp
    MatrixXd randn_e,     // nxp
    MatrixXd &coefs       // bxp
  ):
    Y(Y), X(X), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
    chol_R_list(chol_R_list), h2s_index(h2s_index), tot_Eta_prec(tot_Eta_prec), coefs(coefs) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++){
      int h2_index = h2s_index[j] - 1;
      SpMat chol_R = chol_R_list[h2_index];
      chol_R *= 1/sqrt(tot_Eta_prec[j]);
      coefs.col(j) = sample_MME_single_diagK(Y.col(j), X, prior_mean.col(j), prior_prec.col(j), chol_R, randn_theta.col(j),randn_e.col(j));
    }
  }
};


// -------------------------------------------------- //
// -- Versions of the independent residuals regression model --- //
// -------------------------------------------------- //


// basic version - all regression have same dimensions
// [[Rcpp::export()]]
MatrixXd sample_MME_fixedEffects_c(  // returns bxp matrix
    Map<MatrixXd> Y,              // nxp
    Map<MatrixXd> X,              // nxb
    Rcpp::List Sigma_Choleskys,   // list of Cholesky decompositions of residuals. nxn upper triangular, dgCMatrix
    VectorXi h2s_index,           // px1 index of Cholesky matrix for each column
    Map<VectorXd> tot_Eta_prec,   // px1
    Map<MatrixXd> prior_mean,     // bxp
    Map<MatrixXd> prior_prec,     // bxp
    int grainSize) {

  int n = X.rows();
  int b = X.cols();
  int p = Y.cols();

  MatrixXd randn_theta = rstdnorm_mat(b,p);
  MatrixXd randn_e = rstdnorm_mat(n,p);

  std::vector<MSpMat> chol_R_list;
  for(int i = 0; i < h2s_index.maxCoeff(); i++){
    Rcpp::List Sigma_Choleskys_i = Rcpp::as<Rcpp::List>(Sigma_Choleskys[i]);
    chol_R_list.push_back(Rcpp::as<MSpMat>(Sigma_Choleskys_i["chol_Sigma"]));
  }

  MatrixXd coefs(b,p);

  sample_MME_single_diagK_worker sampler(Y,X,prior_mean,prior_prec,chol_R_list,h2s_index,tot_Eta_prec,randn_theta,randn_e,coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return(coefs);
}
