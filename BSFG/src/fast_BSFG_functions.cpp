#include <math.h>
#include <iostream>
#include "BSFG_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace RcppParallel;

// -------------------------- //
// -------------------------- //
// -------------------------- //

VectorXd sample_coefs_single(
    VectorXd UtEta,
    MatrixXd UtW,
    VectorXd prior_mean,
    VectorXd prior_prec,
    double h2,
    double tot_Eta_prec,
    VectorXd randn_theta,
    VectorXd randn_e,
    VectorXd s,
    int b,
    int n
) {

  VectorXd R_sq_diag = ((h2 * s.array() + (1.0-h2))/tot_Eta_prec).sqrt();
  VectorXd theta_star = randn_theta.array()/prior_prec.array().sqrt();
  theta_star += prior_mean;
  VectorXd e_star = randn_e.array() * R_sq_diag.array();
  MatrixXd UtW_theta_star = UtW * theta_star;
  VectorXd eta_resid = UtEta - UtW_theta_star - e_star;
  MatrixXd RinvSqUtW = R_sq_diag.cwiseInverse().asDiagonal() * UtW;
  VectorXd eta_std = eta_resid.array()/R_sq_diag.array();
  VectorXd WtURinvy = RinvSqUtW.transpose() * eta_std;

  VectorXd theta_tilda;
  if(b < n) {
    MatrixXd C = RinvSqUtW.transpose() * RinvSqUtW;
    C.diagonal() = C.diagonal() + prior_prec;
    theta_tilda = C.householderQr().solve(WtURinvy);
  } else{
    MatrixXd VAi = UtW * prior_prec.cwiseInverse().asDiagonal();
    MatrixXd inner = VAi*UtW.transpose();
    for(int i = 0; i < n; i++) {
      inner(i,i) += (h2 * s(i) + (1.0-h2)) / tot_Eta_prec;
    }
    VectorXd VAiWtURinvy = VAi * WtURinvy;
    VectorXd outerWtURinvy = VAi.transpose() * inner.householderQr().solve(VAiWtURinvy);
    theta_tilda = WtURinvy.array() / prior_prec.array();
    theta_tilda -= outerWtURinvy;
  }

  VectorXd coefs = theta_tilda + theta_star;
  return coefs;
}

// [[Rcpp::export()]]
MatrixXd sample_coefs_parallel_sparse_c_Eigen(
    MSpMat Ut,
    Map<MatrixXd> Eta,
    Map<MatrixXd> W,
    Map<VectorXd> h2,
    Map<VectorXd> tot_Eta_prec,
    Map<VectorXd> s,
    Map<MatrixXd> prior_mean,
    Map<MatrixXd> prior_prec,
    Map<MatrixXd> randn_theta,
    Map<MatrixXd> randn_e,
    int grainSize){

  // Sample regression coefficients
  // columns of matrices are independent
  // each column conditional posterior is a MVN due to conjugacy


  struct sampleColumn : public Worker {
    MatrixXd UtW, UtEta;
    MatrixXd prior_mean, prior_prec, randn_theta, randn_e;
    VectorXd h2, tot_Eta_prec, s;
    int b,n;
    MatrixXd &coefs;

    sampleColumn(MatrixXd UtW, MatrixXd UtEta, MatrixXd prior_mean, MatrixXd prior_prec, VectorXd h2, VectorXd tot_Eta_prec,
                 MatrixXd randn_theta, MatrixXd randn_e,
                 VectorXd s, int b, int n,
                 MatrixXd &coefs) :
      UtW(UtW), UtEta(UtEta), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
      h2(h2), tot_Eta_prec(tot_Eta_prec),
      s(s), b(b), n(n),
      coefs(coefs) {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        coefs.col(j) = sample_coefs_single(UtEta.col(j), UtW, prior_mean.col(j), prior_prec.col(j), h2(j), tot_Eta_prec(j), randn_theta.col(j),randn_e.col(j),s,b,n);
      }
    }
  };

  MatrixXd UtEta = Ut*Eta;
  MatrixXd UtW = Ut*W;

  int p = UtEta.cols();
  int b = UtW.cols();
  int n = UtW.rows();

  MatrixXd coefs(b,p);

  sampleColumn sampler(UtW, UtEta, prior_mean, prior_prec, h2, tot_Eta_prec, randn_theta,randn_e, s, b, n, coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);

  return(coefs);
}

// [[Rcpp::export()]]
MatrixXd sample_coefs_set_c(
    Rcpp::List model_matrices,
    Rcpp::List randn_draws,
    Rcpp::List s_vectors,
    Map<VectorXd> h2s,
    Map<VectorXd> tot_Eta_prec,
    Map<MatrixXd> prior_mean,
    Map<MatrixXd> prior_prec,
    int n,
    int grainSize){

  struct sampleColumn : public Worker {
    std::vector<VectorXd> UtEta_list;
    std::vector<MatrixXd> UtW_list;
    std::vector<VectorXd> randn_theta_list;
    std::vector<VectorXd> randn_e_list;
    std::vector<VectorXd> s_list;
    VectorXd h2s,tot_Eta_prec;
    MatrixXd prior_mean,prior_prec;
    MatrixXd &coefs;

    sampleColumn(
      std::vector<VectorXd> UtEta_list,
      std::vector<MatrixXd> UtW_list,
      std::vector<VectorXd> randn_theta_list,
      std::vector<VectorXd> randn_e_list,
      std::vector<VectorXd> s_list,
      VectorXd h2s,
      VectorXd tot_Eta_prec,
      MatrixXd prior_mean,
      MatrixXd prior_prec,
      MatrixXd &coefs) :
      UtEta_list(UtEta_list), UtW_list(UtW_list), randn_theta_list(randn_theta_list), randn_e_list(randn_e_list),s_list(s_list),
      h2s(h2s), tot_Eta_prec(tot_Eta_prec),prior_mean(prior_mean), prior_prec(prior_prec),
      coefs(coefs)
    {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        int b = randn_theta_list[j].size();
        int n = randn_e_list[j].size();
        coefs.col(j) = sample_coefs_single(UtEta_list[j], UtW_list[j], prior_mean.col(j),
                  prior_prec.col(j), h2s(j), tot_Eta_prec(j), randn_theta_list[j],
                                                                              randn_e_list[j],s_list[j],b,n);
      }
    }
  };

  std::vector<VectorXd> UtEta_list;
  std::vector<MatrixXd> UtW_list;
  std::vector<VectorXd> randn_theta_list;
  std::vector<VectorXd> randn_e_list;
  std::vector<VectorXd> s_list;
  for(int i = 0; i < n; i++){
    Rcpp::List model_matrix_i = Rcpp::as<Rcpp::List>(model_matrices[i]);
    Rcpp::List randn_draws_i = Rcpp::as<Rcpp::List>(randn_draws[i]);
    UtEta_list.push_back(Rcpp::as<VectorXd>(model_matrix_i["y"]));
    UtW_list.push_back(Rcpp::as<MatrixXd>(model_matrix_i["X"]));
    randn_theta_list.push_back(Rcpp::as<VectorXd>(randn_draws_i["randn_theta"]));
    randn_e_list.push_back(Rcpp::as<VectorXd>(randn_draws_i["randn_e"]));
    s_list.push_back(Rcpp::as<VectorXd>(s_vectors[i]));
  }

  int b = randn_theta_list[0].size();
  int p = UtEta_list.size();


  MatrixXd coefs(b,p);

  sampleColumn sampler(UtEta_list, UtW_list,randn_theta_list,randn_e_list,s_list,h2s,tot_Eta_prec,prior_mean,prior_prec,coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);

  return(coefs);
}

// -------------------------- //
// -------------------------- //
// -------------------------- //

// [[Rcpp::export()]]
VectorXd tot_prec_scores_c (
    Map<MatrixXd> UtEta,
    Map<VectorXd> h2,
    Map<VectorXd> s
) {

  int p = UtEta.cols();

  VectorXd scores(p);

  for(int i = 0; i < p; i++){
    ArrayXd Sigma_sqrt = sqrt(h2(i) * s.array() + (1.0 - h2(i)));
    VectorXd SiUtEta_i = UtEta.col(i).array() / Sigma_sqrt;
    scores(i) = SiUtEta_i.dot(SiUtEta_i);
  }
  return scores;
}

// [[Rcpp::export()]]
VectorXd tot_prec_scores_withX_c (
    Map<MatrixXd> UtEta,
    Map<MatrixXd> B_F,
    Map<VectorXd> h2,
    Map<VectorXd> s,
    Map<MatrixXd> prec_B_F
) {

  int p = UtEta.cols();

  VectorXd scores(p);

  for(int i = 0; i < p; i++){
    ArrayXd Sigma_sqrt = sqrt(h2(i) * s.array() + (1.0 - h2(i)));
    VectorXd SiUtEta_i = UtEta.col(i).array() / Sigma_sqrt;
    VectorXd b_std = B_F.col(i).cwiseProduct(prec_B_F.col(i).cwiseSqrt());
    scores(i) = SiUtEta_i.dot(SiUtEta_i) + b_std.dot(b_std);
  }
  return scores;
}

// -------------------------- //
// -------------------------- //
// -------------------------- //

// [[Rcpp::export()]]
MatrixXd log_p_h2s_fast(
    Map<MatrixXd> UtEta,
    Map<VectorXd> tot_Eta_prec,
    Map<VectorXd> discrete_priors,
    Map<VectorXd> s,
    int grainSize)
{

  struct sampleColumn : public Worker {
    MatrixXd std_scores_b2;
    VectorXd discrete_priors;
    VectorXd s;
    MatrixXd &log_ps;

    sampleColumn(MatrixXd std_scores_b2,
                 VectorXd discrete_priors,
                 VectorXd s,
                 MatrixXd &log_ps):
      std_scores_b2(std_scores_b2), discrete_priors(discrete_priors), s(s), log_ps(log_ps) {}

    void operator()(std::size_t begin, std::size_t end) {
      int n = std_scores_b2.rows();
      int b = discrete_priors.size();
      for(std::size_t i = begin; i < end; i++){
        double h2 = double(i)/b;
        VectorXd s2s = h2*s.array() + (1.0-h2);
        MatrixXd std_scores2 = -0.5 * std_scores_b2 * s2s.cwiseInverse().asDiagonal();
        double det = -n/2 * log(2.0*M_PI) - 0.5*s2s.array().log().sum();
        log_ps.row(i) = std_scores2.rowwise().sum().array() + det + log(discrete_priors[i]);
      }
    }
  };

  int b = discrete_priors.size();
  int p = UtEta.cols();

  MatrixXd log_ps(b,p);

  MatrixXd std_scores_b = tot_Eta_prec.cwiseSqrt().asDiagonal() * UtEta.transpose();
  MatrixXd std_scores_b2 = std_scores_b.cwiseProduct(std_scores_b);

  sampleColumn sampler(std_scores_b2,discrete_priors,s, log_ps);
  RcppParallel::parallelFor(0,b,sampler,grainSize);
  return(log_ps);
}

// -------------------------- //
// -------------------------- //
// -------------------------- //



// [[Rcpp::export()]]
MatrixXd sample_randomEffects_parallel_sparse_c_Eigen (
    Map<MatrixXd> Eta,
    MSpMat Z,
    Map<ArrayXd> tot_prec,
    Map<ArrayXd> h2,
    List invert_aZZt_Kinv,
    Map<ArrayXXd> randn_draws,
    int grainSize) {
  //samples genetic effects on factors (F_a) conditional on the factor scores F:
  // F_i = F_{a_i} + E_i, E_i~N(0,s2*(1-h2)*I) for each latent trait i
  // U_i = zeros(r,1) if h2_i = 0
  // it is assumed that s2 = 1 because this scaling factor is absorbed in
  // Lambda
  // invert_aZZt_Kinv has parameters to diagonalize a*Z*Z' + b*K for fast
  // inversion:

  ArrayXd a_prec = tot_prec / h2;
  ArrayXd e_prec = tot_prec / (1.0-h2);

  SpMat U = as<MSpMat>(invert_aZZt_Kinv["U"]);
  VectorXd s1 = as<VectorXd>(invert_aZZt_Kinv["s1"]);
  VectorXd s2 = as<VectorXd>(invert_aZZt_Kinv["s2"]);

  int p = Eta.cols();
  int r = Z.cols();
  // MatrixXd b = U.transpose() * Z.transpose() * (Eta * e_prec.matrix().asDiagonal());
  MatrixXd b = (Z*U).transpose() * (Eta * e_prec.matrix().asDiagonal());

  // MatrixXd z = randn(r,p);

  MatrixXd effects(r,p);

  struct sampleColumn : public Worker {
    VectorXd s1, s2;
    ArrayXd a_prec, e_prec;
    SpMat U;
    MatrixXd b;
    ArrayXXd randn_draws;
    MatrixXd &effects;

    sampleColumn(VectorXd s1, VectorXd s2, ArrayXd a_prec, ArrayXd e_prec, SpMat U, MatrixXd b, ArrayXXd randn_draws, MatrixXd &effects)
      : s1(s1), s2(s2), a_prec(a_prec), e_prec(e_prec), U(U), b(b), randn_draws(randn_draws), effects(effects) {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        ArrayXd d = s2*a_prec(j) + s1*e_prec(j);
        ArrayXd mlam = b.col(j).array() / d;
        effects.col(j) = U * (mlam + randn_draws.col(j) / sqrt(d)).matrix();
      }
    }
  };

  sampleColumn sampler(s1, s2, a_prec, e_prec, U, b, randn_draws, effects);
  RcppParallel::parallelFor(0,p,sampler,grainSize);

  return(effects);
}

// -------------------------- //
// -------------------------- //
// -------------------------- //

// [[Rcpp::export()]]
MatrixXd sample_factors_scores_sparse_c_Eigen(
    Map<MatrixXd> Eta_tilde,
    Map<MatrixXd> prior_mean,
    Map<MatrixXd> Lambda,
    Map<VectorXd> resid_Eta_prec,
    Map<VectorXd> F_e_prec,
    Map<MatrixXd> randn_draws
) {
  //Sample factor scores given factor loadings (F_a), factor residual variances (F_e_prec) and
  //phenotype residuals
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

// -------------------------- //
// -------------------------- //
// -------------------------- //

// [[Rcpp::export()]]
VectorXd sample_delta_c_Eigen(
    VectorXd delta,
    VectorXd tauh,
    Map<MatrixXd> Lambda_prec,
    double delta_1_shape,
    double delta_1_rate,
    double delta_2_shape,
    double delta_2_rate,
    Map<MatrixXd> randg_draws,  // all done with rate = 1;
    Map<MatrixXd> Lambda2
) {
  int times = randg_draws.rows();
  int k = tauh.size();

  MatrixXd scores_mat = Lambda_prec.cwiseProduct(Lambda2);
  VectorXd scores = scores_mat.colwise().sum();

  double rate,delta_old;
  for(int i = 0; i < times; i++){
    delta_old = delta(0);
    rate = delta_1_rate + 0.5 * (1/delta(0)) * tauh.dot(scores);
    delta(0) = randg_draws(i,0) / rate;
    // tauh = cumprod(delta);
    tauh *= delta(0)/delta_old;   // replaces re-calculating cumprod

    for(int h = 1; h < k-1; h++) {
      delta_old = delta(h);
      rate = delta_2_rate + 0.5*(1/delta(h))*tauh.tail(k-h).dot(scores.tail(k-h));
      delta(h) = randg_draws(i,h) / rate;
      // tauh = cumprod(delta);
      tauh.tail(k-h) *= delta(h)/delta_old; // replaces re-calculating cumprod
      // Rcout << (tauh - cumprod(delta)).sum() << std::endl;
    }
  }
  return(delta);
}

