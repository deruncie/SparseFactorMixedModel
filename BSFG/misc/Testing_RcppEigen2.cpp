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


MatrixXd sample_MME_multiple_diagR(
    MatrixXd Y,
    SpMat W,
    SpMat chol_C,
    VectorXd pe,
    MSpMat chol_K_inv,
    VectorXd tot_Eta_prec,
    MatrixXd randn_theta,
    MatrixXd randn_e
){
  MatrixXd theta_star = chol_K_inv.triangularView<Upper>().solve(randn_theta);
  MatrixXd e_star = randn_e * pe.cwiseSqrt().cwiseInverse().asDiagonal();
  MatrixXd W_theta_star = W * theta_star;

  MatrixXd Y_resid = Y - W_theta_star - e_star;
  MatrixXd WtRiy = W.transpose() * (Y_resid * pe.asDiagonal());

  MatrixXd theta_tilda = chol_C.transpose().triangularView<Upper>().solve(chol_C.triangularView<Lower>().solve(WtRiy));

  MatrixXd theta = theta_tilda * tot_Eta_prec.cwiseInverse().asDiagonal() + theta_star;

  return theta;
}

// [[Rcpp::export]]
VectorXd sample_MME_multiple_diagR_c(
    Map<MatrixXd> Y,
    MSpMat W,
    MSpMat chol_C,
    Map<VectorXd> pe,
    MSpMat chol_K_inv,
    Map<VectorXd> tot_Eta_prec,
    Map<MatrixXd> randn_theta,
    Map<MatrixXd> randn_e
){
  return sample_MME_multiple_diagR(Y,W,chol_C,pe,chol_K_inv,tot_Eta_prec,randn_theta,randn_e);
}


VectorXd sample_MME_multiple_diagR2(
    VectorXd Y,
    SpMat W,
    SpMat chol_C,
    double pe,
    SpMat chol_K_inv,
    double tot_Eta_prec,
    VectorXd randn_theta,
    VectorXd randn_e
){
  VectorXd theta_star = chol_K_inv.triangularView<Upper>().solve(randn_theta);
  // VectorXd e_star = randn_e * pe.cwiseSqrt().cwiseInverse().asDiagonal();
  VectorXd e_star = randn_e / sqrt(pe);
  MatrixXd W_theta_star = W * theta_star;

  VectorXd Y_resid = Y - W_theta_star - e_star;
  VectorXd WtRiy = W.transpose() * (Y_resid * pe);

  VectorXd theta_tilda = chol_C.transpose().triangularView<Upper>().solve(chol_C.triangularView<Lower>().solve(WtRiy));

  VectorXd theta = theta_tilda / tot_Eta_prec + theta_star;

  return theta;
}

// [[Rcpp::export]]
VectorXd sample_MME_multiple_diagR_c2(
    Map<VectorXd> Y,
    MSpMat W,
    MSpMat chol_C,
    double pe,
    MSpMat chol_K_inv,
    double tot_Eta_prec,
    Map<VectorXd> randn_theta,
    Map<VectorXd> randn_e
){
  return sample_MME_multiple_diagR2(Y,W,chol_C,pe,chol_K_inv,tot_Eta_prec,randn_theta,randn_e);
}


// // [[Rcpp::export]]
// VectorXd sample_MME_multiple_diagR_c(
//       Map<MatrixXd> Y,
//       MSpMat W,
//       MSpMat chol_C,
//       Map<VectorXd> pe,
//       MSpMat chol_K_inv,
//       Map<VectorXd> tot_Eta_prec,
//       Map<MatrixXd> randn_theta,
//       Map<MatrixXd> randn_e
//   ){
//     return sample_MME_multiple_diagR(Y,W,chol_C,pe,chol_K_inv,tot_Eta_prec,randn_theta,randn_e);
//   }


// [[Rcpp::export()]]
MatrixXd sample_MME_ZKZts_c(
    Map<MatrixXd> Y,
    MSpMat W,
    Map<VectorXd> tot_Eta_prec,
    Rcpp::List chol_Cs,
    Rcpp::List chol_K_invs,
    Map<VectorXd> h2s,
    Rcpp::IntegerVector h2s_index,
    Map<MatrixXd> randn_theta,
    Map<MatrixXd> randn_e,
    int grainSize) {

  struct sampleColumn : public Worker {
    MatrixXd Y, randn_theta, randn_e;
    SpMat W;
    Rcpp::List chol_Cs,chol_K_invs;
    Rcpp::IntegerVector h2s_index;
    VectorXd h2s,pes,tot_Eta_prec;
    MatrixXd &coefs;

    sampleColumn(MatrixXd Y,
                 SpMat W,
                 Rcpp::List chol_Cs,
                 VectorXd pes,
                 VectorXd tot_Eta_prec,
                 Rcpp::List chol_K_invs,
                 Rcpp::IntegerVector h2s_index,
                 MatrixXd randn_theta, MatrixXd randn_e,
                 MatrixXd &coefs):
      Y(Y), W(W), chol_Cs(chol_Cs), pes(pes), tot_Eta_prec(tot_Eta_prec), chol_K_invs(chol_K_invs),
      randn_theta(randn_theta), randn_e(randn_e),
      h2s_index(h2s_index), coefs(coefs) {}

    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        int h2_index = h2s_index[j] - 1;
        SpMat chol_C(Rcpp::as<MSpMat>(chol_Cs[h2_index]));
        SpMat chol_K_inv(Rcpp::as<MSpMat>(chol_K_invs[h2_index]));
        coefs.col(j) = sample_MME_multiple_diagR2(Y.col(j), W, chol_C, pes[j],chol_K_inv, tot_Eta_prec[j], randn_theta.col(j),randn_e.col(j));
      }
    }
  };

  int b = randn_theta.rows();
  int p = randn_theta.cols();

  MatrixXd coefs(b,p);
  VectorXd pes = tot_Eta_prec.array() / (1.0-h2s.colwise().sum().array());

  sampleColumn sampler(Y,W,chol_Cs,pes,tot_Eta_prec,chol_K_invs,h2s_index,randn_theta,randn_e,coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return(coefs);
}

