#include <math.h>
#include <iostream>
#include "BSFG_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace RcppParallel;

NumericVector tot_prec_scores_c2 (
    MatrixXd UtEta,
    VectorXd h2,
    ArrayXd s
) {

  int p = UtEta.cols();

  VectorXd scores(p);

  for(int i = 0; i < p; i++){
    ArrayXd Sigma_sqrt = sqrt(h2(i) * s + (1.0 - h2(i)));
    VectorXd SiUtEta_i = UtEta.col(i).array() / Sigma_sqrt;
    scores(i) = SiUtEta_i.dot(SiUtEta_i);
  }
  return wrap(scores);
}

MatrixXd log_p_h2s_fast2(
    MatrixXd UtEta,
    VectorXd tot_Eta_prec,
    VectorXd discrete_priors,
    VectorXd s,
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

Rcpp::IntegerVector sample_h2s2(
    ArrayXXd log_ps,
    VectorXd rs,
    int grainSize
)
{

  struct sampleColumn : public Worker {
    ArrayXXd log_ps;
    VectorXd rs;
    RVector<int> h2s_index;

    sampleColumn(ArrayXXd log_ps,
                 VectorXd rs,
                 Rcpp::IntegerVector h2s_index):
      log_ps(log_ps), rs(rs), h2s_index(h2s_index) {}

    void operator()(std::size_t begin, std::size_t end) {
      int b = log_ps.rows();
      for(std::size_t j = begin; j < end; j++){
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
  int p = log_ps.cols();

  Rcpp::IntegerVector h2s_index(p);

  sampleColumn sampler(log_ps,rs,h2s_index);
  RcppParallel::parallelFor(0,p,sampler,grainSize);

  return h2s_index;
}


MatrixXd sample_randomEffects_parallel_sparse_c_Eigen2 (
    MatrixXd Eta,
    SpMat Z,
    ArrayXd tot_prec,
    ArrayXd h2,
    List invert_aZZt_Kinv,
    ArrayXXd randn_draws,
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

VectorXd sample_coefs_single2(
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

MatrixXd sample_coefs_parallel_sparse_c_Eigen2(
    SpMat Ut,
    MatrixXd Eta,
    SpMat W,
    VectorXd h2,
    VectorXd tot_Eta_prec,
    VectorXd s,
    MatrixXd prior_mean,
    MatrixXd prior_prec,
    MatrixXd randn_theta,
    MatrixXd randn_e,
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
        coefs.col(j) = sample_coefs_single2(UtEta.col(j), UtW, prior_mean.col(j), prior_prec.col(j), h2(j), tot_Eta_prec(j), randn_theta.col(j),randn_e.col(j),s,b,n);
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


MatrixXd sample_factors_scores_sparse_c_Eigen2(
    MatrixXd Eta_tilde,
    MatrixXd prior_mean,
    MatrixXd Lambda,
    VectorXd resid_Eta_prec,
    VectorXd F_e_prec,
    MatrixXd randn_draws
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

// [[Rcpp::export()]]
void sample_latent_traits_fast(
  const MSpMat X,
  const MSpMat X_F,
  const MSpMat Z,
  const Map<MatrixXd> h2s_matrix,
  const LogicalVector X_F_zero_variance,
  const Map<MatrixXd> Eta,
  const Map<MatrixXd> Lambda,
  Map<MatrixXd> F,
  Map<MatrixXd> B,
  Map<MatrixXd> B_F,
  Map<MatrixXd> U_R,
  Map<MatrixXd> U_F,
  Map<MatrixXd> resid_h2,
  Map<MatrixXd> F_h2,
  Map<MatrixXd> tot_Eta_prec,
  Map<MatrixXd> tot_F_prec,
  Map<MatrixXd> B_F_prec,
  double tot_Eta_prec_shape,
  const NumericVector tot_Eta_prec_rate,
  const double tot_F_prec_shape,
  const double tot_F_prec_rate,
  const Map<VectorXd> h2_priors_resids,
  const Map<VectorXd> h2_priors_factors,
  const MSpMat Ut,
  const Map<ArrayXd> s,
  List invert_aZZt_Kinv,
  int nitt
) {
  int k = Lambda.cols();
  int p = Lambda.rows();
  int n = F.rows();
  int r = Z.cols();
  int b_F = X_F.cols();


  MatrixXd XB = X * B;
  MatrixXd Eta_tilde = Eta - XB - F * Lambda.transpose();
  MatrixXd UtEta_tilde = Ut * Eta_tilde;

  NumericVector scores = tot_prec_scores_c2(UtEta_tilde,resid_h2.row(0),s);
  for(int j = 0; j < p; j++){
    tot_Eta_prec.coeffRef(0,j) = rgamma(1,tot_Eta_prec_shape + n/2,1.0/(tot_Eta_prec_rate[j] + 0.5*scores[j]))[0];
  }

  MatrixXd log_ps = log_p_h2s_fast2(UtEta_tilde,tot_Eta_prec.row(0),h2_priors_resids,s,1);
  VectorXd rs = as<VectorXd>(runif(p,0,1));
  IntegerVector resid_h2_index = sample_h2s2(log_ps,rs,1);
  for(int j = 0; j < p; j++){
    resid_h2.col(j) = h2s_matrix.col(resid_h2_index[j]-1);
  }

  VectorXd randn_v = as<VectorXd>(rnorm(r*p));
  Map<MatrixXd> randn(randn_v.data(),r,p);
  U_R = sample_randomEffects_parallel_sparse_c_Eigen2(Eta_tilde, Z, tot_Eta_prec.row(0), resid_h2.row(0),
                                                      invert_aZZt_Kinv, randn,1);

  ArrayXd resid_Eta_prec = tot_Eta_prec.row(0).array() / (1.0 - resid_h2.row(0).array());

 // B_F
  MatrixXd prior_mean(b_F,k);
  prior_mean.setZero();
  MatrixXd prior_prec = B_F_prec * tot_F_prec.row(0).asDiagonal();

  randn_v = as<VectorXd>(rnorm(b_F*k));
  Map<MatrixXd> randn_theta(randn_v.data(),b_F,k);
  randn_v = as<VectorXd>(rnorm(n*k));
  Map<MatrixXd> randn_e(randn_v.data(),n,k);
  B_F = sample_coefs_parallel_sparse_c_Eigen2(Ut,F,X_F,F_h2.row(0),
                                             tot_F_prec.row(0),s,
                                             prior_mean,prior_prec,
                                             randn_theta,randn_e,1);

  MatrixXd XFBF = X_F * B_F;
  MatrixXd F_tilde = F - XFBF;
  MatrixXd UtF_tilde = Ut * F_tilde;

  scores = tot_prec_scores_c2(UtF_tilde,F_h2.row(0),s);
  for(int j = 0; j < b_F; j++){// [[Rcpp::export]]
    if(!X_F_zero_variance[j]){
      NumericVector scores_B_F = wrap(B_F.row(j).array().pow(2) * B_F_prec.row(j).array());
      scores = scores + scores_B_F;
    }
  }

  for(int j = 0; j < k; j++){
    tot_F_prec.coeffRef(0,j) = rgamma(1,tot_F_prec_shape + n/2 + sum(!X_F_zero_variance)/2,
                                      1.0/(tot_F_prec_rate + 0.5*scores[j]))[0];
  }

  log_ps = log_p_h2s_fast2(UtF_tilde,tot_F_prec.row(0),h2_priors_factors,s,1);
  rs = as<VectorXd>(runif(k,0,1));
  IntegerVector F_h2_index = sample_h2s2(log_ps,rs,1);
  for(int j = 0; j < k; j++){
    F_h2.col(j) = h2s_matrix.col(F_h2_index[j]-1);
  }

  randn_v = as<VectorXd>(rnorm(r*k));
  Map<MatrixXd> randn_uf(randn_v.data(),r,k);
  U_F = sample_randomEffects_parallel_sparse_c_Eigen2(F_tilde, Z, tot_F_prec.row(0), F_h2.row(0),
                                                      invert_aZZt_Kinv, randn_uf,1);

  Eta_tilde = Eta - XB - Z * U_R;

  ArrayXd F_e_prec = tot_F_prec.row(0).array() / (1.0 - F_h2.row(0).array());
  prior_mean = Z * U_F;
  if(b_F > 0) {
    prior_mean += XFBF;
  }
  randn_v = as<VectorXd>(rnorm(n*k));
  Map<MatrixXd> randn_f(randn_v.data(),n,k);
  F = sample_factors_scores_sparse_c_Eigen2( Eta_tilde, prior_mean,Lambda,resid_Eta_prec,F_e_prec,randn_f );

}


