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

using namespace Eigen;
using namespace Rcpp;

// [[Rcpp::export()]]
List test(MSpMat X,int n){
  std::vector<SpMat> mats;
  for(int i = 0; i < n; i++){
    mats.push_back(X);
  }
  return List::create(mats = mats);
}

// [[Rcpp::export()]]
double test2(List l,int n){
  std::vector<MSpMat> mats;
  for(int i = 0; i < n; i++){
    mats.push_back(as<MSpMat>(l[i]));
  }
  return mats[9].coeffRef(1,1);
}

// [[Rcpp::export()]]
double test3(List l,int n){
  std::vector<SpMat> mats;
  for(int i = 0; i < n; i++){
    mats.push_back(as<MSpMat>(l[i]));
  }
  return mats[9].coeffRef(1,1);
}
