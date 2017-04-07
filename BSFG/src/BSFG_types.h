#include <RcppArmadillo.h>
#include <RcppEigen.h>
#include <RcppParallel.h>

using Eigen::Map;               	      // 'Eigen::Maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::ArrayXXd;                  // variable size array, double precision
using Eigen::Upper;
using Eigen::Lower;
typedef Eigen::MappedSparseMatrix<double> MSpMat;
typedef Eigen::SparseMatrix<double> SpMat;

using arma::mat;
using arma::vec;
using arma::sp_mat;
using arma::rowvec;
