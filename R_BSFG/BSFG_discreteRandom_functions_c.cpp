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
mat test_spTimes(mat Y, S4 U){
	sp_mat mU;
	mat res;

	mU = convertSparse(U);

	res = mU.t() * Y;
	return(res);
}

mat sweep_times(mat x, int MARGIN, vec STATS){
	int m = x.n_rows;
	int n = x.n_cols;
	mat sweep_mat;
	if(MARGIN == 1) sweep_mat = repmat(STATS,1,n);
	if(MARGIN == 2) sweep_mat = repmat(STATS.t(),m,1);

	x = x % sweep_mat;

	return(x);
}


// [[Rcpp::export()]]
mat sample_coefs_DR_parallel_sparse_c(
					mat Y,
					mat W,
					vec tot_Y_prec,
					mat prior_mean,
					mat prior_prec,
					List inverse_Sigmas,   // List of inverses of Sigma for each trait
					vec inverseList,
					int grainSize) {
	
	// Sample regression coefficients
	// columns of matrices are independent
	// each column conditional posterior is a MVN due to conjugacy


	struct sampleColumn : public Worker {
		mat W, Y, prior_prec, Z;
		vec inverseList, tot_Y_prec;
		List inverse_Sigmas;
		mat &coefs;

		sampleColumn(mat W, mat Y, mat prior_prec, mat Z, vec inverseList, vec tot_Y_prec, List inverse_Sigmas, mat &coefs) :
			W(W), Y(Y), prior_prec(prior_prec), Z(Z), inverseList(inverseList), tot_Y_prec(tot_Y_prec), inverse_Sigmas(inverse_Sigmas), coefs(coefs) {}

      	void operator()(std::size_t begin, std::size_t end) {
			mat cholQ;
			sp_mat Sigma_inv, Q;
			vec means,v,m,y,z;
			for(std::size_t j = begin; j < end; j++){
				S4 inverse_Sigmas_j = as<S4>(inverse_Sigmas[inverseList(j)]["Sigma_inv"]);
				// chol_Sigma_inv = tot_Y_prec(j)*convertSparse(inverse_Chols_j);
				// Sigma_inv = chol_Sigma_inv * chol_Sigma_inv.t()
				Sigma_inv = tot_Y_prec(j)*convertSparse(inverse_Sigmas_j);
				means = W.t()*Sigma_inv * Y.col(j);
				Q = Sigma_inv;
				Q.diag() += prior_prec.col(j);
				// for(int i = 0; i < prior_prec.n_cols; i++) {
				// 	Q(i,i) += prior_prec(i,j);
				// }

				cholQ = chol(Q);
				v = solve(cholQ.t(),means);
				m = solve(cholQ,v);
				z = Z.col(j);
				y = solve(cholQ,z);

				coefs.col(j) = y + m;
			}
		}
	};

	int p = tot_Y_prec.n_elem;

	mat Z = randn(b,p);
	mat coefs = zeros(b,p);

	sampleColumn sampler(W,Y, Z, inverseList, tot_Y_prec, inverse_Chols, coefs);
	parallelFor(0,p,sampler,grainSize);
	return(coefs);
}

// [[Rcpp::export()]]
vec sample_tot_prec_sparse_c (mat Y,
					   double tot_Y_prec_shape,
					   double tot_Y_prec_rate,
					   List inverse_Sigmas,   // List of inverses of Sigma for each trait
					   vec inverseList,
					  ) {

	int n = Y.n_rows;
	int p = Y.n_cols;

	vec tot_Y_prec = zeros(p);

	for(int i = 0; i < p; i++){
		S4 inverse_Sigmas_i = as<S4>(inverse_Sigmas[inverseList(i)]["Sigma_inv"]);
		sp_mat Sigma_inv = convertSparse(inverse_Sigmas_i);
		vec Yi = Y.col(i);
		vec prec = randg(1,distr_param(tot_Y_prec_shape + n/2, 1.0/(tot_Y_prec_rate + 0.5 * Yi.t() * Sigma_inv * Yi)));
		tot_Y_prec(i) = prec(0);
	}
	return(tot_Y_prec);
}

// [[Rcpp::export()]]
vec sample_discrete_h2s_given_p_sparse_c (mat Y,
						vec tot_prec,
						vec discrete_priors,
					   List inverse_Sigmas   // List of inverses of Sigma for each trait
					  ){

	int p = Y.n_cols;
	int n = Y.n_rows;

	int discrete_bins = inverse_Sigmas.n_elem;

	vec discrete_h2_index = zeros(p);

	for(int i = 0; i < discrete_bins; i++){
		vec scores_2 = zeros(p);
		S4 inverse_Sigmas_i = as<S4>(inverse_Sigmas[i]["Sigma_inv"]);
		sp_mat Sigma_inv_i = convertSparse(inverse_Sigmas_i);

		double det = as<double>(inverse_Sigmas[i]["sqrt_det"]);
		for(int j = 0; j < p; j++) {
			vec Yj = Y.col(j);
			scores_2(j) = tot_prec(j) * Yj.t() * Sigma_inv_i * Y.t();

		}
		log_ps.col(i) = -n/2 * log(2.0*M_PI) + det - 1/2 * sum(scores_2,1) + log(discrete_priors(i));
	}
	for(int j =0; j < p; j++){
		double norm_factor = max(log_ps.row(j))+log(sum(exp(log_ps.row(j)-max(log_ps.row(j)))));
		mat ps_j = exp(log_ps.row(j) - norm_factor);
		log_ps.row(j) = ps_j;
		vec r = randu(1);
		uvec selected = find(repmat(r,1,h2_divisions)>cumsum(ps_j,1));
		discrete_h2_index(j) = double(selected.n_elem);
	}

	return(discrete_h2_index);
}

// [[Rcpp::export()]]
mat sample_randomEffects_parallel_sparse_c (mat Y,
				S4 Z_sparse,
				vec tot_prec,
			    List inverse_Sigmas,   // List of inverses of Sigma for each trait
			    vec inverseList
				int grainSize ) {
	//samples genetic effects on factors (F_a) conditional on the factor scores F:
	// F_i = F_{a_i} + E_i, E_i~N(0,s2*(1-h2)*I) for each latent trait i
	// U_i = zeros(r,1) if h2_i = 0
	// it is assumed that s2 = 1 because this scaling factor is absorbed in
	// Lambda
	// invert_aZZt_Ainv has parameters to diagonalize a*Z*Z' + b*I for fast
	// inversion:

	vec a_prec = tot_prec / h2;
	vec e_prec = tot_prec / (1-h2);

	sp_mat Z = convertSparse(Z_sparse);

	mat U = as<mat>(invert_aZZt_Ainv["U"]);
	vec s1 = as<vec>(invert_aZZt_Ainv["s1"]);
	vec s2 = as<vec>(invert_aZZt_Ainv["s2"]);

	int p = Y.n_cols;
	int r = Z.n_cols;
	mat b = U.t() * Z.t() * sweep_times(Y,2,e_prec);

	mat z = randn(r,p);
	// Environment stats("package:stats");
	// Function rnorm = stats["rnorm"];
	// vec z_v = as<vec>(rnorm(r*k));
	// mat z = reshape(z_v,r,k);

	mat effects = zeros(r,p);

	struct sampleColumn : public Worker {
		vec s1, s2, a_prec, e_prec;
		mat U, b, z;
		mat &effects;

		sampleColumn(vec s1, vec s2, vec a_prec, vec e_prec, mat U, mat b, mat z, mat &effects) 
			: s1(s1), s2(s2), a_prec(a_prec), e_prec(e_prec), U(U), b(b), z(z), effects(effects) {}

      	void operator()(std::size_t begin, std::size_t end) {
			vec d, mlam;
			for(std::size_t j = begin; j < end; j++){
				vec d = s2*a_prec(j) + s1*e_prec(j);
				vec mlam = b.col(j) / d;
				effects.col(j) = U * (mlam + z.col(j)/sqrt(d));
			}
		}
	};

	sampleColumn sampler(s1, s2, a_prec, e_prec, U, b, z, effects);
	parallelFor(0,p,sampler,grainSize);

	return(effects);
}

// [[Rcpp::export()]]
mat sample_factors_scores_ipx_sparse_c(mat Y_tilde,
							S4 Z_sparse,
							mat Lambda,
							vec resid_Y_prec,
							mat F_a,
							vec F_e_prec
							 ) {
//Sample factor scores given factor loadings (F_a), factor residual variances (F_e_prec) and
//phenotype residuals

	sp_mat Z = convertSparse(Z_sparse);

	mat Lmsg = sweep_times(Lambda,1,resid_Y_prec);
	mat S = chol(Lambda.t() * Lmsg + diagmat(F_e_prec));
	mat tS = S.t();

	mat Meta = trans(solve(tS,trans(Y_tilde * Lmsg + sweep_times(Z * F_a,2,F_e_prec))));

	mat Zlams = randn(Meta.n_rows,Meta.n_cols);	
	// Environment stats("package:stats");
	// Function rnorm = stats["rnorm"];
	// vec z = as<vec>(rnorm(Meta.n_rows*Meta.n_cols));
	// mat Zlams = reshape(z,Meta.n_rows,Meta.n_cols);

	mat F = trans(solve(S,trans(Meta + Zlams)));

	return(F);
}
