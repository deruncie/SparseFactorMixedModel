
#include <Rcpp.h>
#include <RcppEigen.h>
using namespace Rcpp;
#include <math.h>
#include <iostream>

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;

// [[Rcpp::export()]]
void addSample(Map<VectorXd> X_, Map<VectorXd> Y,IntegerVector dims,int i){
  // dimensions are in order 2,3,1 or original R array
  if(dims.size() != 3) stop("X should have 3 dimensions");
  if(Y.size() != dims[1]*dims[2]) stop("wrong dimension of Y");
  if(i < 0 || i > dims[1]) stop("index out of bounds");
  Map<MatrixXd> X = Map<MatrixXd>(X_.data(),dims[1]*dims[2],dims[0]);
  X.col(i-1) = Map<VectorXd>(Y.data(),Y.size());
}
