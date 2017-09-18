// [[Rcpp::export]]
VectorXd sample_MME_single_diagR(
    VectorXd y,           // nx1
    SpMat Z,              // nxr dgCMatrix
    SpMat chol_C,         // rxr dgCMatrix upper triangular: chol(ZtRinvZ + diag(Kinv))
    double pe,            // double
    SpMat chol_K_inv,     // rxr dgCMatrix upper triangular: chol(solve(K))
    VectorXd randn_theta, // rx1
    VectorXd randn_e      // nx1
){
  VectorXd theta_star = chol_K_inv.triangularView<Upper>().solve(randn_theta);
  VectorXd e_star = randn_e / sqrt(pe);
  MatrixXd Z_theta_star = Z * theta_star;

  VectorXd y_resid = y - Z_theta_star - e_star;
  VectorXd ZtRiy_resid = Z.transpose() * (y_resid * pe);

  VectorXd theta_tilda = chol_C.triangularView<Upper>().solve(chol_C.transpose().triangularView<Lower>().solve(ZtRiy_resid));

  VectorXd theta = theta_tilda + theta_star;

  return theta;
}

struct sample_MME_single_diagR_worker : public RcppParallel::Worker {
  MatrixXd Y;
  SpMat Z;
  const std::vector<MSpMat> chol_C_list,chol_K_inv_list;
  ArrayXd pes;
  VectorXd tot_Eta_prec;
  VectorXi h2s_index;
  MatrixXd randn_theta, randn_e;
  MatrixXd &coefs;

  sample_MME_single_diagR_worker(MatrixXd Y,              // nxp
                                 SpMat Z,                                   // nxr
                                 const std::vector<MSpMat> chol_C_list,     // std::vector of rxr SpMat upper-triangle
                                 const std::vector<MSpMat> chol_K_inv_list, // std::vector of rxr SpMat upper-triangle
                                 ArrayXd pes,                               // px1
                                 VectorXd tot_Eta_prec,                     // px1
                                 VectorXi h2s_index,                        // px1, 1-based index
                                 MatrixXd randn_theta,                      // rxp
                                 MatrixXd randn_e,                          // nxp
                                 MatrixXd &coefs                            // rxp
  ):
    Y(Y), Z(Z),
    chol_C_list(chol_C_list), chol_K_inv_list(chol_K_inv_list),
    pes(pes), tot_Eta_prec(tot_Eta_prec),h2s_index(h2s_index),
    randn_theta(randn_theta), randn_e(randn_e),
    coefs(coefs) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++){
      int h2_index = h2s_index[j] - 1;
      SpMat chol_C = chol_C_list[h2_index] * sqrt(tot_Eta_prec[j]);  // scale C by tot_Eta_prec[j]. C = Zt*Rinv*Z + Kinv, where Rinv and Kinv are scaled by h2.
      SpMat chol_K_inv = chol_K_inv_list[h2_index];
      chol_K_inv *= sqrt(tot_Eta_prec[j]);
      coefs.col(j) = sample_MME_single_diagR(Y.col(j), Z, chol_C, pes[j],chol_K_inv, randn_theta.col(j),randn_e.col(j));
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
    Rcpp::List randomEffect_C_Choleskys, // List. Each element contains chol_C and chol_K_inv (both rxr dgCMatrix upper-triangle)
    Map<MatrixXd> h2s,                  // n_RE x p
    VectorXi h2s_index,                 // px1
    int grainSize) {

  int n = Y.rows();
  int p = Y.cols();
  int r = Z.cols();

  MatrixXd randn_theta = rstdnorm_mat(r,p);
  MatrixXd randn_e = rstdnorm_mat(n,p);

  std::vector<MSpMat> chol_C_list,chol_K_inv_list;
  for(int i = 0; i < h2s_index.maxCoeff(); i++){
    Rcpp::List randomEffect_C_Cholesky_i = Rcpp::as<Rcpp::List>(randomEffect_C_Choleskys[i]);
    chol_C_list.push_back(Rcpp::as<MSpMat>(randomEffect_C_Cholesky_i["chol_C"]));
    chol_K_inv_list.push_back(Rcpp::as<MSpMat>(randomEffect_C_Cholesky_i["chol_K_inv"]));
  }

  MatrixXd U(r,p);
  ArrayXd h2_e = 1.0 - h2s.colwise().sum().array();
  ArrayXd pes = tot_Eta_prec.array() / h2_e.array();

  sample_MME_single_diagR_worker sampler(Y,Z,chol_C_list,chol_K_inv_list,pes,tot_Eta_prec,h2s_index,randn_theta,randn_e,U);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return(U);
}


// [[Rcpp::export()]]
VectorXd sample_MME_single_diagR2(
    VectorXd y,           // nx1
    SpMat ZQt,              // nxr dgCMatrix
    SpMat Q,              // nxr dgCMatrix
    SpMat chol_S,         // rxr dgCMatrix upper triangular: chol(ZtRinvZ + diag(Kinv))
    double pe,            // double
    VectorXd randn_theta // rx1
){
  VectorXd b = ZQt * y * pe;
  b = chol_S.transpose().triangularView<Lower>().solve(b);
  b += randn_theta;
  b = Q * chol_S.triangularView<Upper>().solve(b);
  return(b);
}


struct sample_MME_single_diagR_worker2 : public RcppParallel::Worker {
  MatrixXd Y;
  SpMat ZQt,Q;
  const std::vector<MSpMat> chol_C_list;
  ArrayXd pes;
  VectorXd tot_Eta_prec;
  VectorXi h2s_index;
  MatrixXd randn_theta;
  MatrixXd &coefs;

  sample_MME_single_diagR_worker2(MatrixXd Y,              // nxp
                                  SpMat ZQt,                                   // nxr
                                  SpMat Q,                                   // nxr
                                  const std::vector<MSpMat> chol_C_list,     // std::vector of rxr SpMat upper-triangle
                                  ArrayXd pes,                               // px1
                                  VectorXd tot_Eta_prec,                     // px1
                                  VectorXi h2s_index,                        // px1, 1-based index
                                  MatrixXd randn_theta,                      // rxp
                                  MatrixXd &coefs                            // rxp
  ):
    Y(Y), ZQt(ZQt),Q(Q),
    chol_C_list(chol_C_list),
    pes(pes), tot_Eta_prec(tot_Eta_prec),h2s_index(h2s_index),
    randn_theta(randn_theta),
    coefs(coefs) {}

  void operator()(std::size_t begin, std::size_t end) {
    for(std::size_t j = begin; j < end; j++){
      int h2_index = h2s_index[j] - 1;
      SpMat chol_C = chol_C_list[h2_index] * sqrt(tot_Eta_prec[j]);  // scale C by tot_Eta_prec[j]. C = Zt*Rinv*Z + Kinv, where Rinv and Kinv are scaled by h2.
      coefs.col(j) = sample_MME_single_diagR2(Y.col(j), ZQt,Q, chol_C, pes[j],randn_theta.col(j));
    }
  }
};


// samples random effects from model:
// Y = ZU + E
// U[,j] ~ N(0,1/tot_Eta_prec[j] * h2[j] * K)
// E[,j] ~ N(0,1/tot_Eta_prec[j] * (1-h2[j]) * I_n)
// For complete data, ie no missing obs.
// [[Rcpp::export()]]
MatrixXd sample_MME_ZKZts_c2(
    Map<MatrixXd> Y,                    // nxp
    MSpMat ZQt,                           // nxr
    MSpMat Q,                           // nxr
    Map<VectorXd> tot_Eta_prec,         // px1
    Rcpp::List randomEffect_C_Choleskys, // List. Each element contains chol_C and chol_K_inv (both rxr dgCMatrix upper-triangle)
    Map<MatrixXd> h2s,                  // n_RE x p
    VectorXi h2s_index,                 // px1
    int grainSize) {

  int p = Y.cols();
  int r = ZQt.rows();

  MatrixXd randn_theta = rstdnorm_mat2(r,p);

  std::vector<MSpMat> chol_C_list;
  for(int i = 0; i < h2s_index.maxCoeff(); i++){
    Rcpp::List randomEffect_C_Cholesky_i = Rcpp::as<Rcpp::List>(randomEffect_C_Choleskys[i]);
    chol_C_list.push_back(Rcpp::as<MSpMat>(randomEffect_C_Cholesky_i["chol_C"]));
  }

  MatrixXd U(r,p);
  ArrayXd h2_e = 1.0 - h2s.colwise().sum().array();
  ArrayXd pes = tot_Eta_prec.array() / h2_e.array();

  sample_MME_single_diagR_worker2 sampler(Y,ZQt,Q,chol_C_list,pes,tot_Eta_prec,h2s_index,randn_theta,U);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return(U);
}
