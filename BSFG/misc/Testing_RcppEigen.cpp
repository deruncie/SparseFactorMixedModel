#include <RcppEigen.h>
#include <RcppParallel.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace RcppParallel;
using Eigen::Map;               	// 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
// using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers
using Eigen::Upper;
using Eigen::Lower;
typedef Eigen::MappedSparseMatrix<double> MSpMat;
typedef Eigen::SparseMatrix<double> SpMat;

// // [[Rcpp::export]]
// VectorXd getEigenValues(Map<MatrixXd> M) {
//   SelfAdjointEigenSolver<MatrixXd> es(M);
//   return es.eigenvalues();
// }
//
// // [[Rcpp::export]]
// VectorXd chol_solve(Map<MatrixXd> R, Map<VectorXd> b) {
//   VectorXd x;
//   x = R.triangularView<Upper>().solve(b);
//   return x;
// }
//
// // [[Rcpp::export]]
// VectorXd chol_sp_solve(MSpMat R, Map<VectorXd> b) {
//   // const SpMat L(Eigen::as<MSpMat>(mL));
//   VectorXd x;
//   x = R.triangularView<Upper>().solve(b);
//   return x;
// }
//
// // [[Rcpp::export]]
// VectorXd chol_sp_solve2(MSpMat L, Map<VectorXd> b) {
//   // const SpMat L(Eigen::as<MSpMat>(mL));
//   VectorXd x;
//   x = L.triangularView<Eigen::Lower>().solve(b);
//   return x;
// }


// // [[Rcpp::export]]
// VectorXd sample_MME_single_diagK_c(
//     const VectorXd prior_prec
// ) {//,
//   return prior_prec.array().sqrt();
// }

VectorXd sample_MME_single_diagK(
      VectorXd y,
      MatrixXd W,
      VectorXd prior_mean,
      VectorXd prior_prec,
      SpMat chol_R,
      double tot_Eta_prec,
      VectorXd randn_theta,
      VectorXd randn_e
  ){

  VectorXd theta_star = randn_theta.array() / prior_prec.cwiseSqrt().array();
  theta_star += prior_mean;
  VectorXd e_star = chol_R * randn_e / sqrt(tot_Eta_prec);
  MatrixXd W_theta_star = W * theta_star;
  VectorXd y_resid = y - W_theta_star - e_star;

  MatrixXd RinvSqW = chol_R.triangularView<Lower>().solve(W);
  VectorXd WtRinvy = RinvSqW.transpose() * chol_R.triangularView<Lower>().solve(y_resid) * tot_Eta_prec;

  VectorXd theta_tilda;

  if(W.cols() < W.rows()) {
    MatrixXd C = RinvSqW.transpose() * RinvSqW * tot_Eta_prec;
    C.diagonal() = C.diagonal() + prior_prec;
    theta_tilda = C.householderQr().solve(WtRinvy);
  } else{
    MatrixXd R = chol_R * chol_R.transpose() / tot_Eta_prec;
    MatrixXd AiU = (W * prior_prec.cwiseInverse().asDiagonal()).transpose();
    MatrixXd R_VAiU = R + W * AiU;
    MatrixXd inner = AiU * R_VAiU.householderQr().solve(AiU.transpose());
    theta_tilda = WtRinvy.array()/prior_prec.array();
    theta_tilda -= inner * WtRinvy;
  }

  VectorXd theta = theta_star + theta_tilda;

  return theta;
}

// [[Rcpp::export]]
VectorXd sample_MME_single_diagK_c(
    Map<VectorXd> y,
    Map<MatrixXd> W,
    Map<VectorXd> prior_mean,
    Map<VectorXd> prior_prec,
    MSpMat chol_R,
    double tot_Eta_prec,
    Map<VectorXd> randn_theta,
    Map<VectorXd> randn_e
){
  return sample_MME_single_diagK(y,W,prior_mean,prior_prec,chol_R,tot_Eta_prec,randn_theta,randn_e);
}

// [[Rcpp::export()]]
MatrixXd sample_MME_fixedEffects_c(
    Map<MatrixXd> Y,
    Map<MatrixXd> W,
    Rcpp::List chol_Rs,
    Rcpp::IntegerVector h2s_index,
    Map<VectorXd> tot_Eta_prec,
    Map<MatrixXd> prior_mean,
    Map<MatrixXd> prior_prec,
    Map<MatrixXd> randn_theta,
    Map<MatrixXd> randn_e,
    int grainSize) {

  struct sampleColumn : public Worker {
    MatrixXd Y, W, prior_mean, prior_prec, randn_theta, randn_e;
    Rcpp::List chol_Rs;
    Rcpp::IntegerVector h2s_index;
    VectorXd tot_Eta_prec;
    MatrixXd &coefs;

    sampleColumn(MatrixXd Y, MatrixXd W, MatrixXd prior_mean, MatrixXd prior_prec,
                 Rcpp::List chol_Rs,
                 Rcpp::IntegerVector h2s_index, VectorXd tot_Eta_prec,
                 MatrixXd randn_theta, MatrixXd randn_e,
                 MatrixXd &coefs):
      Y(Y), W(W), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
      chol_Rs(chol_Rs), h2s_index(h2s_index), tot_Eta_prec(tot_Eta_prec), coefs(coefs) {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        int h2_index = h2s_index[j] - 1;
        SpMat chol_R(Rcpp::as<MSpMat>(chol_Rs[h2_index]));
        coefs.col(j) = sample_MME_single_diagK(Y.col(j), W, prior_mean.col(j), prior_prec.col(j), chol_R, tot_Eta_prec[j], randn_theta.col(j),randn_e.col(j));
      }
    }
  };

  int b = randn_theta.rows();
  int p = randn_theta.cols();

  MatrixXd coefs(b,p);

  sampleColumn sampler(Y,W,prior_mean,prior_prec,chol_Rs,h2s_index,tot_Eta_prec,randn_theta,randn_e, coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return(coefs);
}

