// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
// #include <iostream>
using namespace Rcpp;
using namespace arma;

mat sweep_times(mat x, int MARGIN, vec STATS){
	int m = x.n_rows;
	int n = x.n_cols;
	mat sweep_mat;
	if(MARGIN == 1) sweep_mat = repmat(STATS,1,n);
	if(MARGIN == 2) sweep_mat = repmat(STATS.t(),m,1);

	x = x % sweep_mat;

	return(x);
}    

// // [[Rcpp::export()]]
// vec sample_MME_single_diagR_c(
// 		vec y,
// 		sp_mat W,
// 		sp_mat CiWt,
// 		double pe,
// 		vec prior_mean,
// 		sp_mat chol_A){

// 	int n = y.n_elem;
// 	int n_theta = chol_A.n_rows;
// 	// Environment stats("package:stats");
// 	// Function rnorm = stats["rnorm"];

// 	// vec t = as<vec>(rnorm(n_theta));
// 	// vec e = as<vec>(rnorm(n));
// 	// vec theta_star = prior_mean + chol_A * t;
// 	// vec e_star = e / sqrt(pe);
// 	vec theta_star = prior_mean + chol_A * randn(n_theta);
// 	vec e_star = randn(n) / sqrt(pe);
// 	vec W_theta_star = W * theta_star;

// 	vec y_resid = y - W_theta_star - e_star;
// 	// vec WtRiy = W.t() * (y_resid * pe);

// 	// vec theta_tilda = spsolve(C,WtRiy);
// 	vec theta_tilda = CiWt * y_resid * pe;

// 	vec theta = theta_tilda + theta_star;
// 	return(theta);
// }

// // [[Rcpp::export()]]
// vec sample_MME_single_diagA_c(
// 								vec y,
// 								mat W,
// 								mat C,
// 								mat RinvSqW,
// 								vec prior_mean,
// 								vec prior_prec,
// 								sp_mat chol_R) {

// 	int n = y.n_elem;
// 	int n_theta = prior_prec.n_elem;
// 	Environment Matrix("package:Matrix");
// 	Function solve_sparse = Matrix["solve"];
// 	// Environment stats("package:stats");
// 	// Function rnorm = stats["rnorm"];

// 	vec theta_star = prior_mean + randn(n_theta) / sqrt(prior_prec);
// 	// vec t = as<vec>(rnorm(n_theta));
// 	// vec theta_star = prior_mean + t / sqrt(prior_prec);
// 	// vec e = as<vec>(rnorm(n));
// 	vec e_star = chol_R * randn(n);
// 	// vec e_star = chol_R * e;
// 	vec W_theta_star = W * theta_star;

// 	vec y_resid = y - W_theta_star - e_star;
// 	// mat WtRinvy = RinvSqW.t() * spsolve(chol_R.t(),y_resid);
// 	sp_mat chol_R_t = chol_R.t();
// 	S4 RinvSqY_S4 = solve_sparse(chol_R_t,y_resid);
// 	vec RinvSqY = as<vec>(RinvSqY_S4.slot("x"));

// 	vec WtRinvy = RinvSqW.t() * RinvSqY;

// 	// vec theta_tilda = spsolve(C,WtRinvy);
// 	vec theta_tilda = solve(C,WtRinvy);

// 	vec theta = theta_tilda + theta_star;
// 	return(theta);
// }

// [[Rcpp::export()]]
mat sample_factors_scores_sparse_c(mat Y_tilde,
							sp_mat Z,
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

	mat F = trans(solve(S,trans(Meta + Zlams)));

	return(F);
}
         

// [[Rcpp::export()]]
vec sample_delta_c(
					vec delta,
					vec tauh,
					mat Lambda_prec,
					double delta_1_shape,
					double delta_1_rate,
					double delta_2_shape,
					double delta_2_rate,
					mat Lambda2,
					int times = 1
					) {
	int k = tauh.n_elem;
	mat scores_mat = Lambda_prec % Lambda2;
	int p = Lambda2.n_rows;
	rowvec scores = sum(scores_mat,0);

	double shape, rate;
	vec delta_h;
	for(int i = 0; i < times; i++){
		shape = delta_1_shape + 0.5 * p * k;
		rate = delta_1_rate + 0.5 * (1/delta(0)) * dot(tauh,scores);
		delta_h = randg(1, distr_param(shape, 1/rate));
		delta(0) = delta_h(0);
		tauh = cumprod(delta);

		for(int h = 1; h < k-1; h++) {
			shape = delta_2_shape + 0.5*p*(k-h);
			rate = delta_2_rate + 0.5*(1/delta(h))*dot(tauh.subvec(h, k-1),scores.subvec(h,k-1));
			delta_h = randg(1, distr_param(shape, 1/rate));
			delta(h) = delta_h(0);
			tauh = cumprod(delta);
		}
	}
	return(delta);
}

