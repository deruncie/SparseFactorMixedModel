# include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]

// Note: functions contain commented code to use R's random number generator for testing to ensure identical results to the R functions

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
mat sample_factors_scores_c(mat Y_tilde,
							mat Z_1,
							mat Lambda,
							vec resid_Y_prec,
							mat F_a,
							vec F_h2 ) {
//Sample factor scores given factor loadings (F_a), factor heritabilities (F_h2) and
//phenotype residuals

	mat Lmsg = sweep_times(Lambda,1,resid_Y_prec);
	vec tau_e = 1.0/(1.0-F_h2);
	mat S = chol(Lambda.t() * Lmsg + diagmat(tau_e));
	mat tS = S.t();

	mat Meta = trans(solve(tS,trans(Y_tilde * Lmsg + sweep_times(Z_1 * F_a,2,tau_e))));

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
vec sample_h2s_discrete_c (mat F,
						int h2_divisions,
						vec h2_priors,
						List invert_aI_bZAZ){
	// sample factor heritibilties from a discrete set on [0,1)
	// prior places 50% of the weight at h2=0
	// samples conditional on F, marginalizes over F_a.
	// uses invert_aI_bZAZ.U and invert_aI_bZAZ.s to not have to invert aI + bZAZ
	// each iteration.

	mat U = as<mat>(invert_aI_bZAZ["U"]);
	vec s = as<vec>(invert_aI_bZAZ["s"]);

	int k = F.n_cols;
	vec F_h2 = zeros(k);

	mat log_ps = zeros(k,h2_divisions);
	mat std_scores_b = F.t() * U;

	mat det_mat,std_scores;
	mat det;
	for(double i =0; i < h2_divisions; i+=1){
		double h2 = (i)/(h2_divisions);
		if(h2 > 0) {
			std_scores = 1/sqrt(h2) * sweep_times(std_scores_b,2,1/sqrt(s+(1-h2)/h2));
			det = sum(sum(log((s+(1-h2)/h2)*h2)/2)) * ones(1,k);
		} else {
			std_scores = F.t();
			det = zeros(1,k);
		}
		log_ps.col(i) = sum(dnorm_std_log(std_scores),1) - det.t() + log(h2_priors(i));
	}
	for(int j =0; j < k; j++){
		double norm_factor = max(log_ps.row(j))+log(sum(exp(log_ps.row(j)-max(log_ps.row(j)))));
		mat ps_j = exp(log_ps.row(j) - norm_factor);
		log_ps.row(j) = ps_j;
		vec r = randu(1);
		uvec selected = find(repmat(r,1,h2_divisions)>cumsum(ps_j,1));
		F_h2(j) = double(selected.n_elem)/(h2_divisions);
	}

	return(F_h2);
}

// [[Rcpp::export()]]
vec sample_prec_discrete_conditional_c(mat Y,
									   int h2_divisions,
									   vec h2_priors,
									   List invert_aI_bZAZ,
									   vec res_prec) {
	//sample factor heritibilties conditional on a given residual precision
	//(res_precision)
	//prior given as relative weights on each of h2_divisions points. Doesn't
	//have to be normalized
	//samples conditional on F, marginalizes over F_a.
	//uses invert_aI_bZAZ.U and invert_aI_bZAZ.s to not have to invert aI + bZAZ
	//each iteration.
	//ident_prec is a in the above equation.

	mat U = as<mat>(invert_aI_bZAZ["U"]);
	vec s = as<vec>(invert_aI_bZAZ["s"]);

	int p = Y.n_cols;
	int n = Y.n_rows;
	vec Trait_h2 = zeros(p,1);

	mat log_ps = zeros(p,h2_divisions);
	mat std_scores_b = Y.t() * U;
	mat det_mat,std_scores;
	mat det;
	for(double i =0; i < h2_divisions; i+=1){
		double h2 = (i)/(h2_divisions);
		if(h2 > 0) {
			std_scores = sweep_times(sweep_times(std_scores_b,2,1/sqrt(s+(1-h2)/h2)),1,1/sqrt(h2/(res_prec*(1-h2))));
			det_mat = log((s+(1-h2)/h2) * reshape(h2/(res_prec*(1-h2)),1,p))/2;
			det = sum(det_mat,0).t();
		} else {
			std_scores = sweep_times(Y.t(),1,sqrt(res_prec));
			det = n/2*log(1/res_prec);
		}
		log_ps.col(i) = sum(dnorm_std_log(std_scores),1) - det + log(h2_priors(i));
	}
	for(int j =0; j < p; j++){
		double norm_factor = max(log_ps.row(j))+log(sum(exp(log_ps.row(j)-max(log_ps.row(j)))));
		mat ps_j = exp(log_ps.row(j) - norm_factor);
		log_ps.row(j) = ps_j;
		vec r = randu(1);
		uvec selected = find(repmat(r,1,h2_divisions)>cumsum(ps_j,1));
		Trait_h2(j) = double(selected.n_elem)/(h2_divisions);
	}
	vec Prec = (res_prec % (1-Trait_h2))/Trait_h2;

	return(Prec);
}


// [[Rcpp::export()]]
mat sample_F_a_c (mat F,
				mat Z_1,
				vec F_h2,
				List invert_aZZt_Ainv) {
	//samples genetic effects on factors (F_a) conditional on the factor scores F:
	// F_i = F_{a_i} + E_i, E_i~N(0,s2*(1-h2)*I) for each latent trait i
	// U_i = zeros(r,1) if h2_i = 0
	// it is assumed that s2 = 1 because this scaling factor is absorbed in
	// Lambda
	// invert_aZZt_Ainv has parameters to diagonalize a*Z_1*Z_1' + b*I for fast
	// inversion:

	mat U = as<mat>(invert_aZZt_Ainv["U"]);
	vec s1 = as<vec>(invert_aZZt_Ainv["s1"]);
	vec s2 = as<vec>(invert_aZZt_Ainv["s2"]);

	int k = F.n_cols;
	int r = Z_1.n_cols;
	vec tau_e = 1/(1-F_h2);
	vec tau_u = 1/F_h2;
	mat b = U.t() * Z_1.t() * sweep_times(F,2,tau_e);

	mat z = randn(r,k);
	// Environment stats("package:stats");
	// Function rnorm = stats["rnorm"];
	// vec z_v = as<vec>(rnorm(r*k));
	// mat z = reshape(z_v,r,k);

	mat F_a = zeros(r,k);

	for(int j =0; j<k; j++) {
		if(tau_e(j)==1) {
			F_a.col(j) = zeros(r);
		} else if(tau_e(j) > 10000) {
			F_a.col(j) = F.col(j);
		} else {
			vec d = s2*tau_u(j) + s1*tau_e(j);
			vec mlam = b.col(j) / d;
			F_a.col(j) = U * (mlam + z.col(j)/sqrt(d));
		}
		if(tau_e(j) == 1) {
		   //  disp(F_a.col(j))
		}
	}

	return(F_a);
}