#include <RcppEigen.h>
#include <RcppParallel.h>

// [[Rcpp::depends(RcppEigen)]]

using namespace Rcpp;
using namespace RcppParallel;
using Eigen::Map;               	// 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
// using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers
using Eigen::Upper;
using Eigen::Lower;
typedef Eigen::MappedSparseMatrix<double> MSpMat;
typedef Eigen::SparseMatrix<double> SpMat;

SpMat make_chol_C(SpMat chol_K_inv,VectorXd h2s, SpMat ZtZ){
  SpMat Ki = chol_K_inv.transpose() * chol_K_inv;
  SpMat C = ZtZ;
  C /= (1.0 - h2s.sum());
  C += Ki;
  Eigen::LLT<MatrixXd> chol_C(C);
  MatrixXd chol_C_R = chol_C.matrixL().transpose();
  return chol_C_R.sparseView(0,1e-10);
}

SpMat make_Chol_K(std::vector<SpMat> chol_Ki_mats, VectorXd h2s){
  int h = h2s.size();
  VectorXd sizes(h);
  int total_size = 0;
  for(int i = 0; i < h; i++){
    sizes(i) = chol_Ki_mats[i].rows();
    total_size += sizes(i);
  }
  MatrixXd chol_K_dense(total_size,total_size);
  chol_K_dense.setZero();
  int curr_row = 0;
  int curr_col = 0;
  for(int i = 0; i < h; i++){
    if(h2s[i] == 0) {
      chol_K_dense.block(curr_row,curr_col,sizes(i),sizes(i)).diagonal().setOnes();
      chol_K_dense.block(curr_row,curr_col,sizes(i),sizes(i)).diagonal() /= 0;
    } else{
      chol_K_dense.block(curr_row,curr_col,sizes(i),sizes(i)) = chol_Ki_mats[i];
      chol_K_dense.block(curr_row,curr_col,sizes(i),sizes(i)) /= sqrt(h2s[i]);
    }
    curr_row += sizes(i);
    curr_col += sizes(i);
  }
  SpMat chol_K = chol_K_dense.sparseView(0,1e-10);
  return chol_K;
}

// [[Rcpp::export()]]
SpMat make_Chol_K_R(List chol_Ki_mats_list, VectorXd h2s){
  std::vector<SpMat> chol_Ki_mats;
  int h = h2s.size();
  for(int i = 0; i < h; i++){
    chol_Ki_mats.push_back(Rcpp::as<MSpMat>(chol_Ki_mats_list[i]));
  }
  return make_Chol_K(chol_Ki_mats,h2s);
}


// [[Rcpp::export()]]
SpMat make_chol_C_R(List chol_Ki_mats_list, VectorXd h2s, MSpMat ZtZ){
  std::vector<SpMat> chol_Ki_mats;
  int h = h2s.size();
  for(int i = 0; i < h; i++){
    chol_Ki_mats.push_back(Rcpp::as<MSpMat>(chol_Ki_mats_list[i]));
  }
  return make_chol_C(make_Chol_K(chol_Ki_mats,h2s),h2s,ZtZ);
}

struct randomEffect_C_Cholesky{
  SpMat chol_Ci,chol_K_inv;
};

// [[Rcpp::export()]]
void make_randomEffect_C_Choleskys(List chol_Ki_mats_list, MatrixXd h2s_matrix, MSpMat ZtZ) {
  std::vector<SpMat> chol_Ki_mats;
  int h = h2s_matrix.rows();
  int n = h2s_matrix.cols();
  for(int i = 0; i < h; i++){
    chol_Ki_mats.push_back(Rcpp::as<MSpMat>(chol_Ki_mats_list[i]));
  }
  std::vector<randomEffect_C_Cholesky> randomEffect_C_Cholesky_list;
  for(int j = 0; j < n; j++){
    Rcout << j << " ";
    VectorXd h2s = h2s_matrix.col(j);
    randomEffect_C_Cholesky randomEffect_C_Cholesky_j;
    randomEffect_C_Cholesky_j.chol_K_inv = make_Chol_K(chol_Ki_mats,h2s);
    randomEffect_C_Cholesky_j.chol_Ci = make_chol_C(randomEffect_C_Cholesky_j.chol_K_inv,h2s,ZtZ);
  }
}

// class randomEffect_C_Cholesky {
// public:
//   randomEffect_C_Cholesky(VectorXd h2s_): h2s(h2s_){}
//   void make_mats();
//   SpMat chol_Ci, chol_K_inv;
// private:
//   VectorXd h2s;
// };
// randomEffect_C_Cholesky::make_mats(){
//
// }
//
