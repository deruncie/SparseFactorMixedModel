#include "BSFG_types.h"

using namespace Rcpp;
using namespace Eigen;
using namespace RcppParallel;

#include <math.h>
#include <iostream>
#include "BSFG_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace RcppParallel;


// // [[Rcpp::export()]]
// VectorXd find_candidate_states(
//     MatrixXd h2s_matrix,
//     double step_size,
//     int old_state
// ) {
//   VectorXd dists = (h2s_matrix.colwise() - h2s_matrix.col(old_state)).cwiseAbs().colwise().sum();
//   VectorXd indices(dists.size());
//   int count = 0;
//   for(int i = 0; i < dists.size(); i++){
//     if(dists[i] < step_size & dists[i] > 0) {
//       indices[count] = i;
//       count++;
//     }
//   }
//   if(count == 0) {  // return all indices as candidates
//     for(int i = 0; i < dists.size(); i++){
//       indices[count] = i;
//       count++;
//     }
//   }
//   return indices.head(count);
// }


// [[Rcpp::export()]]
List LDLt_sparse(MSpMat A) {
  Eigen::SimplicialLDLT<SpMat> chol_A;
  chol_A.compute(A);
  MatrixXd I = MatrixXd::Identity(chol_A.rows(), chol_A.rows());
  MatrixXd P = chol_A.permutationP() * I;
  return(List::create(
      Named("P") = P.sparseView(),
      Named("L") =  chol_A.matrixL(),
      Named("d") = chol_A.vectorD()
  ));
}

// [[Rcpp::export()]]
List LDLt_notSparse(Map<MatrixXd> A) {
  Eigen::LDLT<MatrixXd> chol_A;
  chol_A.compute(A);
  MatrixXd I = MatrixXd::Identity(chol_A.rows(), chol_A.rows());
  MatrixXd P = chol_A.transpositionsP() * I;
  VectorXd d = chol_A.vectorD();
  MatrixXd L = chol_A.matrixL();
  return(List::create(
      Named("P") = P.sparseView(),
      Named("L") =L.sparseView(),
      Named("d") = d
  ));
}

SpMat make_C(SpMat chol_K_inv,VectorXd h2s, SpMat ZtZ){
  SpMat Ki = chol_K_inv.transpose() * chol_K_inv;
  SpMat C = ZtZ;
  C /= (1.0 - h2s.sum());
  C += Ki;
  return C;
}

SpMat make_Chol_K(std::vector<SpMat> chol_Ki_mats, VectorXd h2s,double tol){
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
  SpMat chol_K = chol_K_dense.sparseView(0,tol);
  return chol_K;
}

struct randomEffect_C_Cholesky{
  SpMat chol_Ci,chol_K_inv;
  MatrixXd P;
};

class randomEffect_C_Cholesky_database {
public:
  randomEffect_C_Cholesky_database(List chol_Ki_mats_list, MatrixXd h2s_matrix, MSpMat ZtZ,double tol, int grainSize);
  SpMat get_chol_Ci(int i) {return randomEffect_C_Cholesky_list[i-1].chol_Ci;}
  SpMat get_chol_K_inv_i(int i) {return randomEffect_C_Cholesky_list[i-1].chol_K_inv;}
private:
  VectorXd h2s;
  std::vector<randomEffect_C_Cholesky> randomEffect_C_Cholesky_list;
};

randomEffect_C_Cholesky_database::randomEffect_C_Cholesky_database(List chol_Ki_mats_list, MatrixXd h2s_matrix, MSpMat ZtZ,double tol,int grainSize) {
  std::vector<SpMat> chol_Ki_mats;
  int h = h2s_matrix.rows();
  int n = h2s_matrix.cols();
  for(int i = 0; i < h; i++){
    chol_Ki_mats.push_back(Rcpp::as<MSpMat>(chol_Ki_mats_list[i]));
  }

  struct calc_matrices : public Worker {
    std::vector<SpMat> chol_Ki_mats;
    MatrixXd h2s_matrix;
    SpMat ZtZ;
    double tol;
    std::vector<randomEffect_C_Cholesky> &randomEffect_C_Cholesky_list_parallel;

    calc_matrices(std::vector<SpMat> chol_Ki_mats,
                  MatrixXd h2s_matrix,
                  SpMat ZtZ,
                  double tol,
                  std::vector<randomEffect_C_Cholesky> &randomEffect_C_Cholesky_list_parallel
                ):
      chol_Ki_mats(chol_Ki_mats),
      h2s_matrix(h2s_matrix), ZtZ(ZtZ),tol(tol),
      randomEffect_C_Cholesky_list_parallel(randomEffect_C_Cholesky_list_parallel)
      {}
    void operator()(std::size_t begin, std::size_t end) {
        for(std::size_t j = begin; j < end; j++){
          VectorXd h2s = h2s_matrix.col(j);
          randomEffect_C_Cholesky randomEffect_C_Cholesky_j;
          randomEffect_C_Cholesky_j.chol_K_inv = make_Chol_K(chol_Ki_mats,h2s,tol);
          SpMat C = make_C(randomEffect_C_Cholesky_j.chol_K_inv,h2s,ZtZ);
          // Eigen::SimplicialLLT<SpMat> chol_C;   // This might work, but would need to re-build all uses of this matrix in terms of the permutation matrix.
          // chol_C.compute(C);
          // SpMat chol_C_R = chol_C.matrixU();
          // randomEffect_C_Cholesky_j.chol_Ci = chol_C_R;
          // PermutationMatrix<Eigen::Dynamic> perm = chol_C.permutationP();
          // MatrixXd P = perm.toDenseMatrix().cast<double>();
          // randomEffect_C_Cholesky_j.P = P;
          // chol_C_L = (P*chol_C_L).sparseView(0,1e-10);
          // chol_C_R = chol_C_R;//*chol_C.permutationPinv();
          // randomEffect_C_Cholesky_j.chol_Ci = chol_C_L.transpose();
          Eigen::LLT<MatrixXd> chol_C;
          chol_C.compute(C);
          MatrixXd chol_C_R = chol_C.matrixU();
          randomEffect_C_Cholesky_j.chol_Ci = chol_C_R.sparseView(0,tol);
          randomEffect_C_Cholesky_list_parallel[j] = randomEffect_C_Cholesky_j;
        }
    }
  };

  randomEffect_C_Cholesky_list.resize(n);

  calc_matrices calculator(chol_Ki_mats,h2s_matrix,ZtZ,tol,randomEffect_C_Cholesky_list);
  RcppParallel::parallelFor(0,n,calculator,grainSize);
}

SpMat make_Sigma(std::vector<SpMat> ZKZts, VectorXd h2s,double tol){
  int n = ZKZts[0].rows();
  int h = h2s.size();
  MatrixXd R(n,n);
  R.setZero();
  for(int i = 0; i < h; i++){
    R += h2s[i] * ZKZts[i];
  }
  R.diagonal().array() += (1.0-h2s.sum());
  return R.sparseView(0,tol);
}

struct Sigma_Cholesky{
  SpMat chol_Sigma;
  double log_det;
};


class Sigma_Cholesky_database {
public:
  Sigma_Cholesky_database(List chol_Ki_mats_list, MatrixXd h2s_matrix, double tol, int grainSize);
  SpMat get_chol_Sigma(int i) {return Sigma_Cholesky_list[i-1].chol_Sigma;}
  double get_log_det(int i) {return Sigma_Cholesky_list[i-1].log_det;}
private:
  VectorXd h2s;
  std::vector<Sigma_Cholesky> Sigma_Cholesky_list;
};
Sigma_Cholesky_database::Sigma_Cholesky_database(List ZtZs_list, MatrixXd h2s_matrix,double tol, int grainSize) {
  std::vector<SpMat> ZKZts;
  int h = h2s_matrix.rows();
  int n = h2s_matrix.cols();
  for(int i = 0; i < h; i++){
    ZKZts.push_back(Rcpp::as<MSpMat>(ZtZs_list[i]));
  }

  struct calc_matrices : public Worker {
    std::vector<SpMat> ZKZts;
    MatrixXd h2s_matrix;
    double tol;
    std::vector<Sigma_Cholesky> &Sigma_Cholesky_list_parallel;

    calc_matrices(std::vector<SpMat> ZKZts,
                  MatrixXd h2s_matrix,
                  double tol,
                  std::vector<Sigma_Cholesky> &Sigma_Cholesky_list_parallel
    ):
      ZKZts(ZKZts),
      h2s_matrix(h2s_matrix),
      tol(tol),
      Sigma_Cholesky_list_parallel(Sigma_Cholesky_list_parallel)
    {}
    void operator()(std::size_t begin, std::size_t end) {
      for(std::size_t j = begin; j < end; j++){
        VectorXd h2s = h2s_matrix.col(j);

        SpMat Sigma = make_Sigma(ZKZts, h2s,tol);
        Eigen::LLT<MatrixXd> chol_Sigma(Sigma);
        MatrixXd chol_SigmaU = chol_Sigma.matrixU();
        Sigma_Cholesky_list_parallel[j].chol_Sigma = chol_SigmaU.sparseView(0,tol);
        Sigma_Cholesky_list_parallel[j].log_det = 2*chol_SigmaU.diagonal().array().log().sum();
      }
    }
  };

  Sigma_Cholesky_list.resize(n);

  calc_matrices calculator(ZKZts,h2s_matrix,tol,Sigma_Cholesky_list);
  RcppParallel::parallelFor(0,n,calculator,grainSize);
}


RCPP_MODULE(Pre_calculations) {
  class_<randomEffect_C_Cholesky_database>("randomEffect_C_Cholesky_database")
  .constructor<List, MatrixXd, MSpMat,double,int>()
  .method("get_chol_Ci",&randomEffect_C_Cholesky_database::get_chol_Ci)
  .method("get_chol_K_inv_i",&randomEffect_C_Cholesky_database::get_chol_K_inv_i)
  ;
  class_<Sigma_Cholesky_database>("Sigma_Cholesky_database")
    .constructor<List,MatrixXd,double,int>()
    .method("get_chol_Sigma",&Sigma_Cholesky_database::get_chol_Sigma)
    .method("get_log_det",&Sigma_Cholesky_database::get_log_det)
  ;
}
