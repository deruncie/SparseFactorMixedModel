// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <iostream>
using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;

// Enable C++11 via this plugin (Rcpp 0.10.3 or later)
// [[Rcpp::plugins(cpp11)]]

// Note: functions contain commented code to use R's random number generator for testing to ensure identical results to the R functions

sp_mat convertSparse(S4 mat) {         // slight improvement with two non-nested loops
	// from http://gallery.rcpp.org/articles/armadillo-sparse-matrix/

    IntegerVector dims = mat.slot("Dim");
    arma::urowvec i = Rcpp::as<arma::urowvec>(mat.slot("i"));
    arma::urowvec p = Rcpp::as<arma::urowvec>(mat.slot("p"));     
    arma::vec x     = Rcpp::as<arma::vec>(mat.slot("x"));

    int nrow = dims[0], ncol = dims[1];
    arma::sp_mat res(nrow, ncol);

    // create space for values, and copy
    arma::access::rw(res.values) = arma::memory::acquire_chunked<double>(x.size() + 1);
    arma::arrayops::copy(arma::access::rwp(res.values), x.begin(), x.size() + 1);

    // create space for row_indices, and copy
    arma::access::rw(res.row_indices) = arma::memory::acquire_chunked<arma::uword>(i.size() + 1);
    arma::arrayops::copy(arma::access::rwp(res.row_indices), i.begin(), i.size() + 1);
    
    // create space for col_ptrs, and copy 
    arma::access::rw(res.col_ptrs) = arma::memory::acquire<arma::uword>(p.size() + 2);
    arma::arrayops::copy(arma::access::rwp(res.col_ptrs), p.begin(), p.size() + 1);

    // important: set the sentinel as well
    arma::access::rwp(res.col_ptrs)[p.size()+1] = std::numeric_limits<arma::uword>::max();
    
    // set the number of non-zero elements
    arma::access::rw(res.n_nonzero) = x.size();

    return(res);

    // Rcout << "SpMat res:\n" << res << std::endl;
}

// [[Rcpp::export()]]
vec solve_sp_c(S4 S4_mat,vec x){
	sp_mat cA = convertSparse(S4_mat);
	vec result1 = spsolve(cA.t(),x);
	vec result2 = spsolve(cA,x);
	return(result2);
}
// [[Rcpp::export()]]
vec solve_c(mat cA,vec x){
	vec result1 = solve(cA.t(),x);
	vec result2 = solve(cA,result1);
	return(result2);
}