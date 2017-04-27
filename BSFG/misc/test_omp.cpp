#include <omp.h>
#include <Rcpp.h>

// [[Rcpp::export]]
double mean_serial(Rcpp::NumericVector x)
{
  double sum = 0.0;
  for (int i=0; i<x.size(); i++)
    sum += x[i];

  return sum / x.size();
}

// [[Rcpp::export]]
double mean_parallel(Rcpp::NumericVector x)
{
  double sum = 0.0;
#pragma omp parallel for simd shared(x) reduction(+:sum)
  for (int i=0; i<x.size(); i++)
    sum += x[i];

  return sum / x.size();
}
