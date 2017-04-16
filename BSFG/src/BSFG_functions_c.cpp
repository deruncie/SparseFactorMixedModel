// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include "BSFG_types.h"
// #include <iostream>
using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;

// Note: functions contain commented code to use R's random number generator for testing to ensure identical results to the R functions
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
mat sample_coefs_parallel_sparse_c(
    mat UtEta,
    mat UtW,
    vec h2,
    vec tot_Eta_prec,
    vec s,
    mat prior_mean,
    mat prior_prec,
    int grainSize) {

  // Sample regression coefficients
  // columns of matrices are independent
  // each column conditional posterior is a MVN due to conjugacy


  struct sampleColumn : public Worker {
    mat UtW, UtEta, prior_mean, prior_prec, randn_theta, randn_e;
    vec h2, tot_Eta_prec, s;
    int b,n;
    mat &coefs;

    sampleColumn(mat UtW, mat UtEta, mat prior_mean, mat prior_prec, vec h2, vec tot_Eta_prec,
                 mat randn_theta, mat randn_e,
                 vec s, int b, int n,
                 mat &coefs) :
      UtW(UtW), UtEta(UtEta), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
      h2(h2), tot_Eta_prec(tot_Eta_prec),
      s(s), b(b), n(n),
      coefs(coefs) {}

    void operator()(std::size_t begin, std::size_t end) {
      mat RinvSqUtW, WtURinvy;
      vec R_sq_diag, theta_star, e_star, UtW_theta_star, eta_resid, theta_tilda;
      for(std::size_t j = begin; j < end; j++){
        R_sq_diag = sqrt((h2(j) * s + (1-h2(j)))/tot_Eta_prec(j));
        theta_star = randn_theta.col(j)/sqrt(prior_prec.col(j)) + prior_mean.col(j);
        e_star = randn_e.col(j) % R_sq_diag;
        UtW_theta_star = UtW * theta_star;
        eta_resid = UtEta.col(j) - UtW_theta_star - e_star;
        RinvSqUtW = sweep_times(UtW,1,1.0/R_sq_diag);
        WtURinvy = RinvSqUtW.t() * (eta_resid/R_sq_diag);

        if(b < n) {
          mat C = RinvSqUtW.t() * RinvSqUtW ;
          for(int i = 0; i < b; i++) {
            C(i,i) += prior_prec(i,j);
          }
          theta_tilda = solve(C,WtURinvy);
        } else{
          mat VAi = sweep_times(UtW,2,1.0/prior_prec.col(j)); // using Binomial Inverse Theorem if b > n
          mat inner = VAi*UtW.t();
          for(int i = 0; i < n; i++) {
            inner(i,i) += (h2(j) * s(i) + (1-h2(j))) / tot_Eta_prec(j);
          }
          mat outer = VAi.t() * solve(inner,VAi);
          theta_tilda = WtURinvy / prior_prec.col(j) - (outer * WtURinvy);
        }

        coefs.col(j) = theta_tilda + theta_star;
      }
    }
  };

  int p = tot_Eta_prec.n_elem;
  int b = UtW.n_cols;
  int n = UtW.n_rows;

  mat randn_theta = randn(b,p);
  mat randn_e = randn(n,p);
  mat coefs = zeros(b,p);

  sampleColumn sampler(UtW, UtEta, prior_mean, prior_prec, h2, tot_Eta_prec, randn_theta,randn_e, s, b, n, coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return(coefs);
}

// [[Rcpp::export()]]
rowvec sample_tot_prec_sparse_c (mat UtEta,
					   vec h2,
					   vec s,
					   double tot_Eta_prec_shape,
					   double tot_Eta_prec_rate
					  ) {

	int n = UtEta.n_rows;
	int p = UtEta.n_cols;

	vec tot_Eta_prec = zeros(p);

	for(int i = 0; i < p; i++){
		vec Sigma_sqrt = sqrt(h2(i) * s + (1.0 - h2(i)));
		vec SiUtEta_i = UtEta.col(i) / Sigma_sqrt;
		vec prec = randg(1,distr_param(tot_Eta_prec_shape + n/2, 1.0/(tot_Eta_prec_rate + 0.5 * dot(SiUtEta_i,SiUtEta_i))));
		tot_Eta_prec(i) = prec(0);
	}
	return(tot_Eta_prec.t());
}

// [[Rcpp::export()]]
rowvec sample_h2s_discrete_given_p_sparse_c (mat UtEta,
            				int h2_divisions,
            				vec h2_priors,
            				vec Tot_prec,
            				vec s
              ){

	int p = UtEta.n_cols;
	int n = UtEta.n_rows;
	vec h2_index = zeros(p);

	mat log_ps = zeros(p,h2_divisions);
	mat std_scores_b = sweep_times(UtEta.t(),1,sqrt(Tot_prec));

	vec s2s;
	mat scores_2;
	for(double i =0; i < h2_divisions; i+=1){
		double h2 = (i)/(h2_divisions);
		s2s = h2*s + (1-h2);
		scores_2 = -sweep_times(std_scores_b % std_scores_b,2,0.5/s2s);
		double det = -n/2 * log(2.0*M_PI) - 0.5*sum(log(s2s));
		log_ps.col(i) = sum(scores_2,1) + det + log(h2_priors(i));
	}
	for(int j =0; j < p; j++){
		double norm_factor = max(log_ps.row(j))+log(sum(exp(log_ps.row(j)-max(log_ps.row(j)))));
		mat ps_j = exp(log_ps.row(j) - norm_factor);
		log_ps.row(j) = ps_j;
		vec r = randu(1);
		uvec selected = find(repmat(r,1,h2_divisions)>cumsum(ps_j,1));
		// h2(j) = double(selected.n_elem)/(h2_divisions);
		h2_index(j) = selected.n_elem;
	}

	return(h2_index.t() + 1);
}

// [[Rcpp::export()]]
mat sample_randomEffects_parallel_sparse_c (mat Eta,
				sp_mat Z,
				vec tot_prec,
				vec h2,
				List invert_aZZt_Kinv,
				int grainSize ) {
	//samples genetic effects on factors (F_a) conditional on the factor scores F:
	// F_i = F_{a_i} + E_i, E_i~N(0,s2*(1-h2)*I) for each latent trait i
	// U_i = zeros(r,1) if h2_i = 0
	// it is assumed that s2 = 1 because this scaling factor is absorbed in
	// Lambda
	// invert_aZZt_Kinv has parameters to diagonalize a*Z*Z' + b*K for fast
	// inversion:

	vec a_prec = tot_prec / h2;
	vec e_prec = tot_prec / (1-h2);

	sp_mat U = as<sp_mat>(invert_aZZt_Kinv["U"]);
	vec s1 = as<vec>(invert_aZZt_Kinv["s1"]);
	vec s2 = as<vec>(invert_aZZt_Kinv["s2"]);

	int p = Eta.n_cols;
	int r = Z.n_cols;
	mat b = U.t() * Z.t() * sweep_times(Eta,2,e_prec);

	mat z = randn(r,p);
	// Environment stats("package:stats");
	// Function rnorm = stats["rnorm"];
	// vec z_v = as<vec>(rnorm(r*k));
	// mat z = reshape(z_v,r,k);

	mat effects = zeros(r,p);

	struct sampleColumn : public Worker {
		vec s1, s2, a_prec, e_prec;
		sp_mat U;
		mat b, z;
		mat &effects;

		sampleColumn(vec s1, vec s2, vec a_prec, vec e_prec, sp_mat U, mat b, mat z, mat &effects)
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
	RcppParallel::parallelFor(0,p,sampler,grainSize);

	return(effects);
}


// [[Rcpp::export()]]
mat sample_factors_scores_sparse_c(mat Eta_tilde,
                                         mat prior_mean,
                                         mat Lambda,
                                         vec resid_Eta_prec,
                                         vec F_e_prec
) {
  //Sample factor scores given factor loadings (F_a), factor residual variances (F_e_prec) and
  //phenotype residuals
  mat Lmsg = sweep_times(Lambda,1,resid_Eta_prec);
  mat S = chol(Lambda.t() * Lmsg + diagmat(F_e_prec));
  mat tS = S.t();

  mat Meta = trans(solve(tS,trans(Eta_tilde * Lmsg + sweep_times(prior_mean,2,F_e_prec))));

  mat Zlams = randn(Meta.n_rows,Meta.n_cols);

  mat F = trans(solve(S,trans(Meta + Zlams)));

  return(F);
}

// [[Rcpp::export()]]
List GSVD_2_c(mat K, mat B){
	int n = B.n_cols;
	mat U,V;
	vec d;
	svd(U,d,V,K * inv(B));

	vec norm_factor = sqrt(ones(n) + d%d);

	mat C = diagmat(d / norm_factor);
	mat S = diagmat(1/norm_factor);
	mat X = sweep_times(B.t() * V,2,norm_factor);

	return(List::create(_["U"]=U,_["V"]=V,_["X"]=X,_["C"]=C,_["S"]=S));
}

// [[Rcpp::export()]]
rowvec sample_delta_c(
    vec delta,
    vec tauh,
    mat Lambda_prec,
    double delta_1_shape,
    double delta_1_rate,
    double delta_2_shape,
    double delta_2_rate,
    mat randg_draws,  // all done with rate = 1;
    mat Lambda2
) {
  int times = randg_draws.n_rows;
  int k = tauh.n_elem;
  mat scores_mat = Lambda_prec % Lambda2;
  rowvec scores = sum(scores_mat,0);

  double rate;
  vec delta_h;
  for(int i = 0; i < times; i++){
    rate = delta_1_rate + 0.5 * (1/delta(0)) * dot(tauh,scores);
    delta_h = randg_draws(i,1) / rate;
    delta(0) = delta_h(0);
    tauh = cumprod(delta);

    for(int h = 1; h < k-1; h++) {
      rate = delta_2_rate + 0.5*(1/delta(h))*dot(tauh.subvec(h, k-1),scores.subvec(h,k-1));
      delta_h = randg_draws(i,h) / rate;
      delta(h) = delta_h(0);
      tauh = cumprod(delta);
    }
  }
  return(delta.t());
}

// [[Rcpp::export()]]
double log_binom_c(
                      vec beta,
                      mat X,
                      vec y,
                      vec N,
                      vec mu,
                      vec sigma2
) {
  vec Xbeta = X * beta;
  vec e_Xbeta = exp(Xbeta);
  vec p = e_Xbeta / (1.0 + e_Xbeta);
  double ll = sum(y%log(p) + (N-y) % log(1.0-p)) - sum((beta - mu) % (beta - mu) / (2*sigma2));
  return ll;
}


