#include <math.h>
#include <iostream>
#include "BSFG_types.h"

// [[Rcpp::depends(RcppEigen)]]
using namespace Eigen;
using namespace RcppParallel;


// [[Rcpp::export()]]
VectorXd sample_delta_c_Eigen(
    VectorXd delta,
    VectorXd tauh,
    Map<VectorXd> scores,
    double delta_1_rate,
    double delta_2_rate,
    Map<MatrixXd> randg_draws  // all done with rate = 1;
) {
  int times = randg_draws.rows();
  int k = tauh.size();

  double rate,delta_old;
  for(int i = 0; i < times; i++){
    delta_old = delta(0);
    rate = delta_1_rate + (1/delta(0)) * tauh.dot(scores);
    delta(0) = randg_draws(i,0) / rate;
    // tauh = cumprod(delta);
    tauh *= delta(0)/delta_old;   // replaces re-calculating cumprod

    for(int h = 1; h < k; h++) {
      delta_old = delta(h);
      rate = delta_2_rate + (1/delta(h))*tauh.tail(k-h).dot(scores.tail(k-h));
      delta(h) = randg_draws(i,h) / rate;
      // tauh = cumprod(delta);
      tauh.tail(k-h) *= delta(h)/delta_old; // replaces re-calculating cumprod
      // Rcout << (tauh - cumprod(delta)).sum() << std::endl;
    }
  }
  return(delta);
}


// -------------------------- //
// -------------------------- //
// -------------------------- //

// [[Rcpp::export()]]
NumericVector rgig_multiple(
  int n,
  NumericVector lambda,
  NumericVector chi,
  NumericVector psi
){
  NumericVector res(n);
  SEXP (*fun)(SEXP,SEXP,SEXP,SEXP) = NULL;
  if (!fun)
    fun = (SEXP(*)(SEXP,SEXP,SEXP,SEXP)) R_GetCCallable("GIGrvg", "rgig");
  for(int i = 0; i < n; i++){
    res[i] = as<double>(fun(wrap(1.0),wrap(lambda[i]),wrap(chi[i]),wrap(psi[i])));
  }
  return(res);
}
