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
// [[Rcpp::export()]]
mat test_spTimes2(mat Y, mat U){
	mat res;
	res = U.t() * Y;
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

mat dnorm_std_log(mat x){
	return(log(1/sqrt(2*M_PI))-0.5* (x % x));
}

// [[Rcpp::export()]]
List GSVD_2_c(mat A, mat B){
	int n = B.n_cols;
	mat U,V;
	vec d;
	svd(U,d,V,A * inv(B));

	vec norm_factor = sqrt(ones(n) + d%d);

	mat C = diagmat(d / norm_factor);
	mat S = diagmat(1/norm_factor);
	mat X = sweep_times(B.t() * V,2,norm_factor);

	return(List::create(_["U"]=U,_["V"]=V,_["X"]=X,_["C"]=C,_["S"]=S));
}

// [[Rcpp::export()]]
mat sample_means_c(mat Y_tilde,
				   vec resid_Y_prec,
				   vec E_a_prec,
				   List invert_aPXA_bDesignDesignT ) {
	// when used to sample [B;E_a]:
	//  W - F*Lambda' = X*B + Z_1*E_a + E, vec(E)~N(0,kron(Psi_E,In)). 
	//  Note: conditioning on F, Lambda and W.
	// The vector [b_j;E_{a_j}] is sampled simultaneously. Each trait is sampled separately because their
	// conditional posteriors factor into independent MVNs.
	// note:invert_aPXA_bDesignDesignT has parameters to diagonalize mixed model equations for fast inversion: 
	// inv(a*blkdiag(fixed_effects_prec*eye(b),Ainv) + b*[X Z_1]'[X Z_1]) = U*diag(1./(a.*s1+b.*s2))*U'
	// Design_U = [X Z_1]*U, which doesn't change each iteration. 
	
	mat U = as<mat>(invert_aPXA_bDesignDesignT["U"]);
	vec s1 = as<vec>(invert_aPXA_bDesignDesignT["s1"]);
	vec s2 = as<vec>(invert_aPXA_bDesignDesignT["s2"]);
	mat Design_U = as<mat>(invert_aPXA_bDesignDesignT["Design_U"]);

	// int n = Y_tilde.n_rows;
	int p = Y_tilde.n_cols;
	int br = Design_U.n_cols;

	mat means = sweep_times(Design_U.t() * Y_tilde,2,resid_Y_prec);
	mat location_sample = zeros(br,p);

	mat Zlams = randn(br,p);	
	// Environment stats("package:stats");
	// Function rnorm = stats["rnorm"];
	// vec z = as<vec>(rnorm(br*p));
	// mat Zlams = reshape(z,br,p);

	vec d, mlam;
	for(int j =0; j<p; j++) {
		d = s1*E_a_prec(j) + s2*resid_Y_prec(j);
		mlam = means.col(j) /d;
		location_sample.col(j) = U * (mlam + Zlams.col(j)/sqrt(d));
	}

	return(location_sample);
}

// [[Rcpp::export()]]
mat sample_means_parallel_c(mat Y_tilde,
				   vec resid_Y_prec,
				   vec E_a_prec,
				   List invert_aPXA_bDesignDesignT,
				   int grainSize ) {
	// when used to sample [B;E_a]:
	//  W - F*Lambda' = X*B + Z_1*E_a + E, vec(E)~N(0,kron(Psi_E,In)). 
	//  Note: conditioning on F, Lambda and W.
	// The vector [b_j;E_{a_j}] is sampled simultaneously. Each trait is sampled separately because their
	// conditional posteriors factor into independent MVNs.
	// note:invert_aPXA_bDesignDesignT has parameters to diagonalize mixed model equations for fast inversion: 
	// inv(a*blkdiag(fixed_effects_prec*eye(b),Ainv) + b*[X Z_1]'[X Z_1]) = U*diag(1./(a.*s1+b.*s2))*U'
	// Design_U = [X Z_1]*U, which doesn't change each iteration. 
	
	struct sampleColumn : public Worker {
		vec E_a_prec, resid_Y_prec, s1, s2;
		mat means, Zlams, U;

		mat &location_sample;

		sampleColumn(vec E_a_prec, vec resid_Y_prec, vec s1, vec s2, mat means, mat Zlams, mat U, mat &location_sample)
			: E_a_prec(E_a_prec), resid_Y_prec(resid_Y_prec), s1(s1), s2(s2), means(means), Zlams(Zlams),U(U), location_sample(location_sample) {}

      	void operator()(std::size_t begin, std::size_t end) {
			vec d, mlam;
			for(std::size_t j = begin; j < end; j++){
				d = s1*E_a_prec(j) + s2*resid_Y_prec(j);
				mlam = means.col(j) /d;
				location_sample.col(j) = U * (mlam + Zlams.col(j)/sqrt(d));
			}
		}
	};

	mat U = as<mat>(invert_aPXA_bDesignDesignT["U"]);
	vec s1 = as<vec>(invert_aPXA_bDesignDesignT["s1"]);
	vec s2 = as<vec>(invert_aPXA_bDesignDesignT["s2"]);
	mat Design_U = as<mat>(invert_aPXA_bDesignDesignT["Design_U"]);

	// int n = Y_tilde.n_rows;
	int p = Y_tilde.n_cols;
	int br = Design_U.n_cols;

	mat means = sweep_times(Design_U.t() * Y_tilde,2,resid_Y_prec);
	mat location_sample = zeros(br,p);

	mat Zlams = randn(br,p);	
	// Environment stats("package:stats");
	// Function rnorm = stats["rnorm"];
	// vec z = as<vec>(rnorm(br*p));
	// mat Zlams = reshape(z,br,p);

	// vec d, mlam;
	// for(int j =0; j<p; j++) {
	// 	d = s1*E_a_prec(j) + s2*resid_Y_prec(j);
	// 	mlam = means.col(j) /d;
	// 	location_sample.col(j) = U * (mlam + Zlams.col(j)/sqrt(d));
	// }
	sampleColumn sampler(E_a_prec,resid_Y_prec,s1,s2,means,Zlams,U,location_sample);
	parallelFor(0,p,sampler,grainSize);

	return(location_sample);
}

// [[Rcpp::export()]]
mat sample_coefs_c(
					mat Y,
					mat W,
					vec h2,
					vec tot_Y_prec,
					mat prior_mean,
					mat prior_prec,
					List invert_aI_bZAZ) {
	
	// Sample regression coefficients
	// columns of matrices are independent
	// each column conditional posterior is a MVN due to conjugacy
	int p = tot_Y_prec.n_elem;
	int b = W.n_cols;

	mat U = as<mat>(invert_aI_bZAZ["U"]);
	vec s = as<vec>(invert_aI_bZAZ["s"]);

	mat WtU = W.t() * U;
	mat UtY = U.t() * Y;

	mat Z = randn(b,p);
	mat coefs = zeros(b,p);

	mat WtUDi, Q, cholQ;
	vec means,v,m,y,z;
	for(int j =0; j<p; j++){
		WtUDi = tot_Y_prec(j) * sweep_times(WtU,2,1.0/(h2(j) * s + (1-h2(j))));
		means = WtUDi * UtY.col(j);
		Q = WtUDi * WtU.t();
		for(int i = 0; i < b; i++) {
			Q(i,i) += prior_prec(i,j);
		}

		cholQ = chol(Q);
		v = solve(cholQ.t(),means);
		m = solve(cholQ,v);
		z = Z.col(j);
		y = solve(cholQ,z);

		coefs.col(j) = y + m;
	}
	return(coefs);
}

// [[Rcpp::export()]]
mat sample_coefs_parallel_c(
					mat Y,
					mat W,
					vec h2,
					vec tot_Y_prec,
					mat prior_mean,
					mat prior_prec,
					List invert_aI_bZAZ,
					int grainSize) {
	
	// Sample regression coefficients
	// columns of matrices are independent
	// each column conditional posterior is a MVN due to conjugacy


	struct sampleColumn : public Worker {
		mat WtU, UtY, prior_prec, Z;
		vec h2, tot_Y_prec, s;
		int b;
		mat &coefs;

		sampleColumn(mat WtU, mat UtY, mat prior_prec, mat Z, vec h2, vec tot_Y_prec, vec s, int b, mat &coefs) :
			WtU(WtU), UtY(UtY), prior_prec(prior_prec), Z(Z), h2(h2), tot_Y_prec(tot_Y_prec), s(s), b(b), coefs(coefs) {}

      	void operator()(std::size_t begin, std::size_t end) {
			mat WtUDi, Q, cholQ;
			vec means,v,m,y,z;
			for(std::size_t j = begin; j < end; j++){
				WtUDi = tot_Y_prec(j) * sweep_times(WtU,2,1.0/(h2(j) * s + (1-h2(j))));
				means = WtUDi * UtY.col(j);
				Q = WtUDi * WtU.t();
				for(int i = 0; i < b; i++) {
					Q(i,i) += prior_prec(i,j);
				}

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
	int b = W.n_cols;

	mat U = as<mat>(invert_aI_bZAZ["U"]);
	vec s = as<vec>(invert_aI_bZAZ["s"]);

	mat WtU = W.t() * U;
	mat UtY = U.t() * Y;

	mat Z = randn(b,p);
	mat coefs = zeros(b,p);

	sampleColumn sampler(WtU, UtY, prior_prec, Z, h2, tot_Y_prec, s, b, coefs);
	parallelFor(0,p,sampler,grainSize);
	return(coefs);
}

// [[Rcpp::export()]]
mat sample_coefs_parallel_sparse_c(
					mat Y,
					mat W,
					vec h2,
					vec tot_Y_prec,
					mat prior_mean,
					mat prior_prec,
					List invert_aI_bZAZ,
					int grainSize) {
	
	// Sample regression coefficients
	// columns of matrices are independent
	// each column conditional posterior is a MVN due to conjugacy


	struct sampleColumn : public Worker {
		mat WtU, UtY, prior_prec, Z;
		vec h2, tot_Y_prec, s;
		int b;
		mat &coefs;

		sampleColumn(mat WtU, mat UtY, mat prior_prec, mat Z, vec h2, vec tot_Y_prec, vec s, int b, mat &coefs) :
			WtU(WtU), UtY(UtY), prior_prec(prior_prec), Z(Z), h2(h2), tot_Y_prec(tot_Y_prec), s(s), b(b), coefs(coefs) {}

      	void operator()(std::size_t begin, std::size_t end) {
			mat WtUDi, Q, cholQ;
			vec means,v,m,y,z;
			for(std::size_t j = begin; j < end; j++){
				WtUDi = tot_Y_prec(j) * sweep_times(WtU,2,1.0/(h2(j) * s + (1-h2(j))));
				means = WtUDi * UtY.col(j);
				Q = WtUDi * WtU.t();
				for(int i = 0; i < b; i++) {
					Q(i,i) += prior_prec(i,j);
				}

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
	int b = W.n_cols;

	S4 U_sparse = as<S4>(invert_aI_bZAZ["U_sparse"]);
	vec s = as<vec>(invert_aI_bZAZ["s"]);

	sp_mat U = convertSparse(U_sparse);

	mat WtU = W.t() * U;
	mat UtY = U.t() * Y;

	mat Z = randn(b,p);
	mat coefs = zeros(b,p);

	sampleColumn sampler(WtU, UtY, prior_prec, Z, h2, tot_Y_prec, s, b, coefs);
	parallelFor(0,p,sampler,grainSize);
	return(coefs);
}

// [[Rcpp::export()]]
mat sample_Lambda_c(mat Y_tilde,
					mat F,
					vec resid_Y_prec,
					vec E_a_prec,
					mat Plam,
					List invert_aI_bZAZ ){
	
	// Sample factor loadings Lambda while marginalizing over residual
	// genetic effects: Y - Z_2W = F*Lambda' + E, vec(E)~N(0,kron(Psi_E,In) + kron(Psi_U, ZAZ^T))
	// note: conditioning on F, but marginalizing over E_a.
	// sampling is done separately by trait because each column of Lambda is
	// independent in the conditional posterior
	// note: invert_aI_bZAZ has parameters that diagonalize aI + bZAZ for fast
	// inversion: inv(aI + bZAZ) = 1/b*U*diag(1./(s+a/b))*U'
	
	int p = resid_Y_prec.n_elem;
	int k = F.n_cols;

	mat U = as<mat>(invert_aI_bZAZ["U"]);
	vec s = as<vec>(invert_aI_bZAZ["s"]);

	mat FtU = F.t() * U;
	mat UtY = U.t() * Y_tilde;

	mat Zlams = randn(k,p);	
	// Environment stats("package:stats");
	// Function rnorm = stats["rnorm"];
	// vec z = as<vec>(rnorm(k*p));
	// mat Zlams = reshape(z,k,p);
	mat Lambda = zeros(p,k);

	mat FUDi, Qlam, Llam_t,Plam_row;
	vec means,vlam,mlam,ylam,Zcol;
	for(int j =0; j<p; j++){
		FUDi = E_a_prec(j) * sweep_times(FtU, 2,1/(s + E_a_prec(j)/resid_Y_prec(j)));
		means = FUDi * UtY.col(j);
		Qlam = FUDi * FtU.t() + diagmat(Plam.row(j));

		Llam_t = chol(Qlam);
		vlam = solve(Llam_t.t(),means);
		mlam = solve(Llam_t,vlam);
		Zcol = Zlams.col(j);
		ylam = solve(Llam_t,Zcol);

		Lambda.row(j) = ylam.t() + mlam.t();

	}

	return(Lambda);
}

// [[Rcpp::export()]]
mat sample_Lambda_parallel_c(mat Y_tilde,
					mat F,
					vec resid_Y_prec,
					vec E_a_prec,
					mat Plam,
					List invert_aI_bZAZ,
					int grainSize ){
	
	// Sample factor loadings Lambda while marginalizing over residual
	// genetic effects: Y - Z_2W = F*Lambda' + E, vec(E)~N(0,kron(Psi_E,In) + kron(Psi_U, ZAZ^T))
	// note: conditioning on F, but marginalizing over E_a.
	// sampling is done separately by trait because each column of Lambda is
	// independent in the conditional posterior
	// note: invert_aI_bZAZ has parameters that diagonalize aI + bZAZ for fast
	// inversion: inv(aI + bZAZ) = 1/b*U*diag(1./(s+a/b))*U'
	
	struct sampleColumn : public Worker {

		vec E_a_prec, s, resid_Y_prec;
		mat FtU, UtY, Plam, Zlams;
		mat &Lambda;

		sampleColumn(vec E_a_prec, 
					 vec s,
					 vec resid_Y_prec,
					 mat FtU,
					 mat UtY,
					 mat Plam,
					 mat Zlams,
					 mat &Lambda)
			: E_a_prec(E_a_prec), s(s), resid_Y_prec(resid_Y_prec),FtU(FtU),UtY(UtY),Plam(Plam),
				Zlams(Zlams),
				Lambda(Lambda) {
				Lambda.zeros();
			}

      	void operator()(std::size_t begin, std::size_t end) {
			mat FUDi, Qlam, Llam_t,Plam_row;
			vec means,vlam,mlam,ylam,Zcol;
			for(std::size_t j = begin; j < end; j++){
				FUDi = E_a_prec(j) * sweep_times(FtU, 2,1/(s + E_a_prec(j)/resid_Y_prec(j)));
				means = FUDi * UtY.col(j);
				Qlam = FUDi * FtU.t() + diagmat(Plam.row(j));

				Llam_t = chol(Qlam);
				vlam = solve(Llam_t.t(),means);
				mlam = solve(Llam_t,vlam);
				// Zcol = randn(Llam_t.n_cols);
				Zcol = Zlams.col(j);
				ylam = solve(Llam_t,Zcol);

				Lambda.row(j) = ylam.t() + mlam.t();
			}
		}
	};	

	int p = resid_Y_prec.n_elem;
	int k = F.n_cols;

	mat U = as<mat>(invert_aI_bZAZ["U"]);
	vec s = as<vec>(invert_aI_bZAZ["s"]);

	mat FtU = F.t() * U;
	mat UtY = U.t() * Y_tilde;

	mat Zlams = randn(k,p);	
	// Environment stats("package:stats");
	// Function rnorm = stats["rnorm"];
	// vec z = as<vec>(rnorm(k*p));
	// mat Zlams = reshape(z,k,p);
	mat Lambda = zeros(p,k);

	mat FUDi, Qlam, Llam_t,Plam_row;
	vec means,vlam,mlam,ylam,Zcol;

	sampleColumn sampler(E_a_prec,s,resid_Y_prec,FtU,UtY,Plam,
		Zlams,
		Lambda);
	parallelFor(0,p,sampler,grainSize);

	return(Lambda);
}

// [[Rcpp::export()]]
mat sample_Lambda_parallel_sparse_c(mat Y_tilde,
					mat F,
					vec resid_Y_prec,
					vec E_a_prec,
					mat Plam,
					List invert_aI_bZAZ,
					int grainSize ){
	
	// Sample factor loadings Lambda while marginalizing over residual
	// genetic effects: Y - Z_2W = F*Lambda' + E, vec(E)~N(0,kron(Psi_E,In) + kron(Psi_U, ZAZ^T))
	// note: conditioning on F, but marginalizing over E_a.
	// sampling is done separately by trait because each column of Lambda is
	// independent in the conditional posterior
	// note: invert_aI_bZAZ has parameters that diagonalize aI + bZAZ for fast
	// inversion: inv(aI + bZAZ) = 1/b*U*diag(1./(s+a/b))*U'
	
	struct sampleColumn : public Worker {

		vec E_a_prec, s, resid_Y_prec;
		mat FtU, UtY, Plam, Zlams;
		mat &Lambda;

		sampleColumn(vec E_a_prec, 
					 vec s,
					 vec resid_Y_prec,
					 mat FtU,
					 mat UtY,
					 mat Plam,
					 mat Zlams,
					 mat &Lambda)
			: E_a_prec(E_a_prec), s(s), resid_Y_prec(resid_Y_prec),FtU(FtU),UtY(UtY),Plam(Plam),
				Zlams(Zlams),
				Lambda(Lambda) {
				Lambda.zeros();
			}

      	void operator()(std::size_t begin, std::size_t end) {
			mat FUDi, Qlam, Llam_t,Plam_row;
			vec means,vlam,mlam,ylam,Zcol;
			for(std::size_t j = begin; j < end; j++){
				FUDi = E_a_prec(j) * sweep_times(FtU, 2,1/(s + E_a_prec(j)/resid_Y_prec(j)));
				means = FUDi * UtY.col(j);
				Qlam = FUDi * FtU.t() + diagmat(Plam.row(j));

				Llam_t = chol(Qlam);
				vlam = solve(Llam_t.t(),means);
				mlam = solve(Llam_t,vlam);
				// Zcol = randn(Llam_t.n_cols);
				Zcol = Zlams.col(j);
				ylam = solve(Llam_t,Zcol);

				Lambda.row(j) = ylam.t() + mlam.t();
			}
		}
	};	

	int p = resid_Y_prec.n_elem;
	int k = F.n_cols;

	S4 U_sparse = as<S4>(invert_aI_bZAZ["U_sparse"]);
	vec s = as<vec>(invert_aI_bZAZ["s"]);

	sp_mat U = convertSparse(U_sparse);

	mat FtU = F.t() * U;
	mat UtY = U.t() * Y_tilde;

	mat Zlams = randn(k,p);	
	// Environment stats("package:stats");
	// Function rnorm = stats["rnorm"];
	// vec z = as<vec>(rnorm(k*p));
	// mat Zlams = reshape(z,k,p);
	mat Lambda = zeros(p,k);

	mat FUDi, Qlam, Llam_t,Plam_row;
	vec means,vlam,mlam,ylam,Zcol;

	sampleColumn sampler(E_a_prec,s,resid_Y_prec,FtU,UtY,Plam,
		Zlams,
		Lambda);
	parallelFor(0,p,sampler,grainSize);

	return(Lambda);
}

// [[Rcpp::export()]]
mat sample_factors_scores_ipx_c(mat Y_tilde,
							mat Z,
							mat Lambda,
							vec resid_Y_prec,
							mat F_a,
							vec F_e_prec
							 ) {
//Sample factor scores given factor loadings (F_a), factor residual variances (F_e_prec) and
//phenotype residuals

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

// [[Rcpp::export()]]
mat sample_px_factors_scores_c(mat Y_tilde,
							mat Z_1,
							mat Lambda_px,
							vec resid_Y_prec,
							mat F_a_px,
							vec tau_e,
							vec px_factor ) {
//Sample factor scores given factor loadings (F_a), factor heritabilities (F_h2) and
//phenotype residuals
// as before, but px_factor gives the variance of the px_factor	

	mat Lmsg = sweep_times(Lambda_px,1,resid_Y_prec);
	// vec tau_e = 1.0 / (px_factor .* (1.0 - F_h2));
	mat S = chol(Lambda_px.t() * Lmsg + diagmat(tau_e));
	mat tS = S.t();

	mat Meta = trans(solve(tS,trans(Y_tilde * Lmsg + sweep_times(Z_1 * F_a_px,2,tau_e))));

	mat Zlams = randn(Meta.n_rows,Meta.n_cols);	
	// Environment stats("package:stats");
	// Function rnorm = stats["rnorm"];
	// vec z = as<vec>(rnorm(Meta.n_rows*Meta.n_cols));
	// mat Zlams = reshape(z,Meta.n_rows,Meta.n_cols);

	mat F_px = trans(solve(S,trans(Meta + Zlams)));

	return(F_px);
}


// [[Rcpp::export()]]
vec sample_tot_prec_c (mat Y,
					   vec h2,
					   double tot_Y_prec_shape,
					   double tot_Y_prec_rate,
					   List invert_aI_bZAZ
					  ) {

	mat U = as<mat>(invert_aI_bZAZ["U"]);
	vec s = as<vec>(invert_aI_bZAZ["s"]);

	mat UtY = U.t() * Y;

	int n = Y.n_rows;
	int p = Y.n_cols;

	vec tot_Y_prec = zeros(p);

	for(int i = 0; i < p; i++){
		vec Sigma_sqrt = sqrt(h2(i) * s + (1.0 - h2(i)));
		vec SiUtY_i = UtY.col(i) / Sigma_sqrt;
		vec prec = randg(1,distr_param(tot_Y_prec_shape + n/2, 1.0/(tot_Y_prec_rate + 0.5 * dot(SiUtY_i,SiUtY_i))));
		tot_Y_prec(i) = prec(0);
	}
	return(tot_Y_prec);
}

// [[Rcpp::export()]]
vec sample_tot_prec_sparse_c (mat Y,
					   vec h2,
					   double tot_Y_prec_shape,
					   double tot_Y_prec_rate,
					   List invert_aI_bZAZ
					  ) {

	S4 U_sparse = as<S4>(invert_aI_bZAZ["U_sparse"]);
	vec s = as<vec>(invert_aI_bZAZ["s"]);

	sp_mat U = convertSparse(U_sparse);

	mat UtY = U.t() * Y;

	int n = Y.n_rows;
	int p = Y.n_cols;

	vec tot_Y_prec = zeros(p);

	for(int i = 0; i < p; i++){
		vec Sigma_sqrt = sqrt(h2(i) * s + (1.0 - h2(i)));
		vec SiUtY_i = UtY.col(i) / Sigma_sqrt;
		vec prec = randg(1,distr_param(tot_Y_prec_shape + n/2, 1.0/(tot_Y_prec_rate + 0.5 * dot(SiUtY_i,SiUtY_i))));
		tot_Y_prec(i) = prec(0);
	}
	return(tot_Y_prec);
}

// [[Rcpp::export()]]
vec sample_h2s_discrete_given_p_c (mat Y,
						int h2_divisions,
						vec h2_priors,
						vec Tot_prec,
						List invert_aI_bZAZ){

	mat U = as<mat>(invert_aI_bZAZ["U"]);
	vec s = as<vec>(invert_aI_bZAZ["s"]);

	int p = Y.n_cols;
	vec h2 = zeros(p);

	mat log_ps = zeros(p,h2_divisions);
	mat std_scores_b = sweep_times(Y.t() * U,1,sqrt(Tot_prec));

	mat det_mat,std_scores;
	mat det;
	for(double i =0; i < h2_divisions; i+=1){
		double h2 = (i)/(h2_divisions);
		if(h2 > 0) {
			std_scores = sweep_times(std_scores_b,2,1/sqrt(h2*s + (1-h2)));
			det = sum(sum(log((h2*s + (1-h2)))/2)) * ones(1,p);
		} else {
			std_scores = Y.t();
			det = zeros(1,p);
		}
		log_ps.col(i) = sum(dnorm_std_log(std_scores),1) - det.t() + log(h2_priors(i));
	}
	for(int j =0; j < p; j++){
		double norm_factor = max(log_ps.row(j))+log(sum(exp(log_ps.row(j)-max(log_ps.row(j)))));
		mat ps_j = exp(log_ps.row(j) - norm_factor);
		log_ps.row(j) = ps_j;
		vec r = randu(1);
		uvec selected = find(repmat(r,1,h2_divisions)>cumsum(ps_j,1));
		h2(j) = double(selected.n_elem)/(h2_divisions);
	}

	return(h2);
}
// [[Rcpp::export()]]
vec sample_h2s_discrete_given_p_sparse_c (mat Y,
						int h2_divisions,
						vec h2_priors,
						vec Tot_prec,
						List invert_aI_bZAZ){

	S4 U_sparse = as<S4>(invert_aI_bZAZ["U_sparse"]);
	vec s = as<vec>(invert_aI_bZAZ["s"]);

	sp_mat U = convertSparse(U_sparse);

	int p = Y.n_cols;
	vec h2 = zeros(p);

	mat log_ps = zeros(p,h2_divisions);
	mat std_scores_b = sweep_times(Y.t() * U,1,sqrt(Tot_prec));

	mat det_mat,std_scores;
	mat det;
	for(double i =0; i < h2_divisions; i+=1){
		double h2 = (i)/(h2_divisions);
		if(h2 > 0) {
			std_scores = sweep_times(std_scores_b,2,1/sqrt(h2*s + (1-h2)));
			det = sum(sum(log((h2*s + (1-h2)))/2)) * ones(1,p);
		} else {
			std_scores = Y.t();
			det = zeros(1,p);
		}
		log_ps.col(i) = sum(dnorm_std_log(std_scores),1) - det.t() + log(h2_priors(i));
	}
	for(int j =0; j < p; j++){
		double norm_factor = max(log_ps.row(j))+log(sum(exp(log_ps.row(j)-max(log_ps.row(j)))));
		mat ps_j = exp(log_ps.row(j) - norm_factor);
		log_ps.row(j) = ps_j;
		vec r = randu(1);
		uvec selected = find(repmat(r,1,h2_divisions)>cumsum(ps_j,1));
		h2(j) = double(selected.n_elem)/(h2_divisions);
	}

	return(h2);
}

// [[Rcpp::export()]]
mat sample_F_a_ipx_c (mat F,
				mat Z,
				vec F_a_prec,
				vec F_e_prec,
				List invert_aZZt_Ainv) {
	//samples genetic effects on factors (F_a) conditional on the factor scores F:
	// F_i = F_{a_i} + E_i, E_i~N(0,s2*(1-h2)*I) for each latent trait i
	// U_i = zeros(r,1) if h2_i = 0
	// it is assumed that s2 = 1 because this scaling factor is absorbed in
	// Lambda
	// invert_aZZt_Ainv has parameters to diagonalize a*Z*Z' + b*I for fast
	// inversion:

	mat U = as<mat>(invert_aZZt_Ainv["U"]);
	vec s1 = as<vec>(invert_aZZt_Ainv["s1"]);
	vec s2 = as<vec>(invert_aZZt_Ainv["s2"]);

	int k = F.n_cols;
	int r = Z.n_cols;
	mat b = U.t() * Z.t() * sweep_times(F,2,F_e_prec);

	mat z = randn(r,k);
	// Environment stats("package:stats");
	// Function rnorm = stats["rnorm"];
	// vec z_v = as<vec>(rnorm(r*k));
	// mat z = reshape(z_v,r,k);

	mat F_a = zeros(r,k);

	for(int j =0; j<k; j++) {
		if(F_a_prec(j) > 10000) {
			F_a.col(j) = zeros(r);
		} else {
			vec d = s2*F_a_prec(j) + s1*F_e_prec(j);
			vec mlam = b.col(j) / d;
			F_a.col(j) = U * (mlam + z.col(j)/sqrt(d));
		}
	}

	return(F_a);
}


// [[Rcpp::export()]]
mat sample_F_a_ipx_sparse_c (mat F,
				S4 Z_sparse,
				vec F_a_prec,
				vec F_e_prec,
				List invert_aZZt_Ainv) {
	//samples genetic effects on factors (F_a) conditional on the factor scores F:
	// F_i = F_{a_i} + E_i, E_i~N(0,s2*(1-h2)*I) for each latent trait i
	// U_i = zeros(r,1) if h2_i = 0
	// it is assumed that s2 = 1 because this scaling factor is absorbed in
	// Lambda
	// invert_aZZt_Ainv has parameters to diagonalize a*Z*Z' + b*I for fast
	// inversion:

	sp_mat Z = convertSparse(Z_sparse);

	mat U = as<mat>(invert_aZZt_Ainv["U"]);
	vec s1 = as<vec>(invert_aZZt_Ainv["s1"]);
	vec s2 = as<vec>(invert_aZZt_Ainv["s2"]);

	int k = F.n_cols;
	int r = Z.n_cols;
	mat b = U.t() * Z.t() * sweep_times(F,2,F_e_prec);

	mat z = randn(r,k);
	// Environment stats("package:stats");
	// Function rnorm = stats["rnorm"];
	// vec z_v = as<vec>(rnorm(r*k));
	// mat z = reshape(z_v,r,k);

	mat F_a = zeros(r,k);

	for(int j =0; j<k; j++) {
		if(F_a_prec(j) > 10000) {
			F_a.col(j) = zeros(r);
		} else {
			vec d = s2*F_a_prec(j) + s1*F_e_prec(j);
			vec mlam = b.col(j) / d;
			F_a.col(j) = U * (mlam + z.col(j)/sqrt(d));
		}
	}

	return(F_a);
}