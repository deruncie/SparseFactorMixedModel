#include <RcppEigen.h>
#include <RcppParallel.h>
#include <GIGrvg.h>


using Eigen::Map;               	      // 'Eigen::Maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::VectorXi;                  // variable size vector, int
using Eigen::ArrayXXd;                  // variable size array, double precision
using Eigen::ArrayXd;                  // variable size array, double precision
using Eigen::Upper;
using Eigen::Lower;
typedef Eigen::SparseMatrix<double> SpMat;
typedef Eigen::Map<SpMat> MSpMat;

using namespace Rcpp;
using namespace RcppParallel;
using namespace Eigen;

VectorXd find_candidate_states(MatrixXd, double, int);
