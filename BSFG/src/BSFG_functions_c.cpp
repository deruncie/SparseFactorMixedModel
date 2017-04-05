// [[Rcpp::depends(RcppParallel)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppParallel.h>
#include <iostream>
using namespace Rcpp;
using namespace RcppParallel;
using namespace arma;

// Note: functions contain commented code to use R's random number generator for testing to ensure identical results to the R functions
arma::mat sweep_times(arma::mat x, int MARGIN, arma::vec STATS){
  int m = x.n_rows;
  int n = x.n_cols;
  arma::mat sweep_mat;
  if(MARGIN == 1) sweep_mat = repmat(STATS,1,n);
  if(MARGIN == 2) sweep_mat = repmat(STATS.t(),m,1);

  x = x % sweep_mat;

  return(x);
}
// Note: functions contain commented code to use R's random number generator for testing to ensure identical results to the R functions
arma::sp_mat sweep_times_sp(arma::sp_mat x, int MARGIN, arma::vec STATS){
  int m = x.n_rows;
  int n = x.n_cols;
  arma::mat sweep_mat;
  if(MARGIN == 1) sweep_mat = repmat(STATS,1,n);
  if(MARGIN == 2) sweep_mat = repmat(STATS.t(),m,1);

  x = x % sweep_mat;

  return(x);
}

// [[Rcpp::export()]]
arma::mat sample_coefs_parallel_sparse_c(
    arma::mat UtEta,
    arma::mat UtW,
    arma::vec h2,
    arma::vec tot_Eta_prec,
    arma::vec s,
    arma::mat prior_mean,
    arma::mat prior_prec,
    int grainSize) {

  // Sample regression coefficients
  // columns of matrices are independent
  // each column conditional posterior is a MVN due to conjugacy


  struct sampleColumn : public Worker {
    arma::mat UtW, UtEta, prior_mean, prior_prec, randn_theta, randn_e;
    arma::vec h2, tot_Eta_prec, s;
    int b,n;
    arma::mat &coefs;

    sampleColumn(arma::mat UtW, arma::mat UtEta, arma::mat prior_mean, arma::mat prior_prec, arma::vec h2, arma::vec tot_Eta_prec,
                 arma::mat randn_theta, arma::mat randn_e,
                 arma::vec s, int b, int n,
                 arma::mat &coefs) :
      UtW(UtW), UtEta(UtEta), prior_mean(prior_mean), prior_prec(prior_prec), randn_theta(randn_theta), randn_e(randn_e),
      h2(h2), tot_Eta_prec(tot_Eta_prec),
      s(s), b(b), n(n),
      coefs(coefs) {}

    void operator()(std::size_t begin, std::size_t end) {
      arma::mat RinvSqUtW, WtURinvy;
      arma::vec R_sq_diag, theta_star, e_star, UtW_theta_star, eta_resid, theta_tilda;
      for(std::size_t j = begin; j < end; j++){
        R_sq_diag = sqrt((h2(j) * s + (1-h2(j)))/tot_Eta_prec(j));
        theta_star = randn_theta.col(j)/sqrt(prior_prec.col(j)) + prior_mean.col(j);
        e_star = randn_e.col(j) % R_sq_diag;
        UtW_theta_star = UtW * theta_star;
        eta_resid = UtEta.col(j) - UtW_theta_star - e_star;
        RinvSqUtW = sweep_times(UtW,1,1.0/R_sq_diag);
        WtURinvy = RinvSqUtW.t() * (eta_resid/R_sq_diag);

        if(b < n) {
          arma::mat C = RinvSqUtW.t() * RinvSqUtW ;
          for(int i = 0; i < b; i++) {
            C(i,i) += prior_prec(i,j);
          }
          theta_tilda = solve(C,WtURinvy);
        } else{
          arma::mat VAi = sweep_times(UtW,2,1.0/prior_prec.col(j)); // using Binomial Inverse Theorem if b > n
          arma::mat inner = VAi*UtW.t();
          for(int i = 0; i < n; i++) {
            inner(i,i) += (h2(j) * s(i) + (1-h2(j))) / tot_Eta_prec(j);
          }
          arma::mat outer = VAi.t() * solve(inner,VAi);
          theta_tilda = WtURinvy / prior_prec.col(j) - (outer * WtURinvy);
        }

        coefs.col(j) = theta_tilda + theta_star;
      }
    }
  };

  int p = tot_Eta_prec.n_elem;
  int b = UtW.n_cols;
  int n = UtW.n_rows;

  arma::mat randn_theta = randn(b,p);
  arma::mat randn_e = randn(n,p);
  arma::mat coefs = zeros(b,p);

  sampleColumn sampler(UtW, UtEta, prior_mean, prior_prec, h2, tot_Eta_prec, randn_theta,randn_e, s, b, n, coefs);
  RcppParallel::parallelFor(0,p,sampler,grainSize);
  return(coefs);
}

// [[Rcpp::export()]]
arma::rowvec sample_tot_prec_sparse_c (arma::mat Eta,
					   arma::vec h2,
					   double tot_Eta_prec_shape,
					   double tot_Eta_prec_rate,
					   List invert_aI_bZKZ
					  ) {

	arma::sp_mat U = as<arma::sp_mat>(invert_aI_bZKZ["U"]);
	arma::vec s = as<arma::vec>(invert_aI_bZKZ["s"]);

	arma::mat UtEta = U.t() * Eta;

	int n = Eta.n_rows;
	int p = Eta.n_cols;

	arma::vec tot_Eta_prec = zeros(p);

	for(int i = 0; i < p; i++){
		arma::vec Sigma_sqrt = sqrt(h2(i) * s + (1.0 - h2(i)));
		arma::vec SiUtEta_i = UtEta.col(i) / Sigma_sqrt;
		arma::vec prec = randg(1,distr_param(tot_Eta_prec_shape + n/2, 1.0/(tot_Eta_prec_rate + 0.5 * dot(SiUtEta_i,SiUtEta_i))));
		tot_Eta_prec(i) = prec(0);
	}
	return(tot_Eta_prec.t());
}

// [[Rcpp::export()]]
arma::rowvec sample_h2s_discrete_given_p_sparse_c (arma::mat Eta,
						int h2_divisions,
						arma::vec h2_priors,
						arma::vec Tot_prec,
						List invert_aI_bZKZ){

	arma::sp_mat U = as<arma::sp_mat>(invert_aI_bZKZ["U"]);
	arma::vec s = as<arma::vec>(invert_aI_bZKZ["s"]);

	int p = Eta.n_cols;
	int n = Eta.n_rows;
	arma::vec h2_index = zeros(p);

	arma::mat log_ps = zeros(p,h2_divisions);
	arma::mat std_scores_b = sweep_times(Eta.t() * U,1,sqrt(Tot_prec));

	arma::vec s2s;
	arma::mat scores_2;
	for(double i =0; i < h2_divisions; i+=1){
		double h2 = (i)/(h2_divisions);
		s2s = h2*s + (1-h2);
		scores_2 = -sweep_times(std_scores_b % std_scores_b,2,0.5/s2s);
		double det = -n/2 * log(2.0*M_PI) - 0.5*sum(log(s2s));
		log_ps.col(i) = sum(scores_2,1) + det + log(h2_priors(i));
	}
	for(int j =0; j < p; j++){
		double norm_factor = max(log_ps.row(j))+log(sum(exp(log_ps.row(j)-max(log_ps.row(j)))));
		arma::mat ps_j = exp(log_ps.row(j) - norm_factor);
		log_ps.row(j) = ps_j;
		arma::vec r = randu(1);
		uvec selected = find(repmat(r,1,h2_divisions)>cumsum(ps_j,1));
		// h2(j) = double(selected.n_elem)/(h2_divisions);
		h2_index(j) = selected.n_elem;
	}

	return(h2_index.t() + 1);
}

// [[Rcpp::export()]]
arma::mat sample_randomEffects_parallel_sparse_c (arma::mat Eta,
				arma::sp_mat Z,
				arma::vec tot_prec,
				arma::vec h2,
				List invert_aZZt_Kinv,
				int grainSize ) {
	//samples genetic effects on factors (F_a) conditional on the factor scores F:
	// F_i = F_{a_i} + E_i, E_i~N(0,s2*(1-h2)*I) for each latent trait i
	// U_i = zeros(r,1) if h2_i = 0
	// it is assumed that s2 = 1 because this scaling factor is absorbed in
	// Lambda
	// invert_aZZt_Kinv has parameters to diagonalize a*Z*Z' + b*I for fast
	// inversion:

	arma::vec a_prec = tot_prec / h2;
	arma::vec e_prec = tot_prec / (1-h2);

	arma::sp_mat U = as<arma::sp_mat>(invert_aZZt_Kinv["U"]);
	arma::vec s1 = as<arma::vec>(invert_aZZt_Kinv["s1"]);
	arma::vec s2 = as<arma::vec>(invert_aZZt_Kinv["s2"]);

	int p = Eta.n_cols;
	int r = Z.n_cols;
	arma::mat b = U.t() * Z.t() * sweep_times(Eta,2,e_prec);

	arma::mat z = randn(r,p);
	// Environment stats("package:stats");
	// Function rnorm = stats["rnorm"];
	// arma::vec z_v = as<arma::vec>(rnorm(r*k));
	// arma::mat z = reshape(z_v,r,k);

	arma::mat effects = zeros(r,p);

	struct sampleColumn : public Worker {
		arma::vec s1, s2, a_prec, e_prec;
		arma::sp_mat U;
		arma::mat b, z;
		arma::mat &effects;

		sampleColumn(arma::vec s1, arma::vec s2, arma::vec a_prec, arma::vec e_prec, arma::sp_mat U, arma::mat b, arma::mat z, arma::mat &effects)
			: s1(s1), s2(s2), a_prec(a_prec), e_prec(e_prec), U(U), b(b), z(z), effects(effects) {}

      	void operator()(std::size_t begin, std::size_t end) {
			arma::vec d, mlam;
			for(std::size_t j = begin; j < end; j++){
				arma::vec d = s2*a_prec(j) + s1*e_prec(j);
				arma::vec mlam = b.col(j) / d;
				effects.col(j) = U * (mlam + z.col(j)/sqrt(d));
			}
		}
	};

	sampleColumn sampler(s1, s2, a_prec, e_prec, U, b, z, effects);
	RcppParallel::parallelFor(0,p,sampler,grainSize);

	return(effects);
}


// [[Rcpp::export()]]
arma::mat sample_factors_scores_sparse_c(arma::mat Eta_tilde,
                                         arma::mat prior_mean,
                                         arma::mat Lambda,
                                         arma::vec resid_Eta_prec,
                                         arma::vec F_e_prec
) {
  //Sample factor scores given factor loadings (F_a), factor residual variances (F_e_prec) and
  //phenotype residuals
  arma::mat Lmsg = sweep_times(Lambda,1,resid_Eta_prec);
  arma::mat S = chol(Lambda.t() * Lmsg + diagmat(F_e_prec));
  arma::mat tS = S.t();

  arma::mat Meta = trans(solve(tS,trans(Eta_tilde * Lmsg + sweep_times(prior_mean,2,F_e_prec))));

  arma::mat Zlams = randn(Meta.n_rows,Meta.n_cols);

  arma::mat F = trans(solve(S,trans(Meta + Zlams)));

  return(F);
}

// [[Rcpp::export()]]
List GSVD_2_c(arma::mat K, arma::mat B){
	int n = B.n_cols;
	arma::mat U,V;
	arma::vec d;
	svd(U,d,V,K * inv(B));

	arma::vec norm_factor = sqrt(ones(n) + d%d);

	arma::mat C = diagmat(d / norm_factor);
	arma::mat S = diagmat(1/norm_factor);
	arma::mat X = sweep_times(B.t() * V,2,norm_factor);

	return(List::create(_["U"]=U,_["V"]=V,_["X"]=X,_["C"]=C,_["S"]=S));
}

// [[Rcpp::export()]]
arma::rowvec sample_delta_c(
					arma::vec delta,
					arma::vec tauh,
					arma::mat Lambda_prec,
					double delta_1_shape,
					double delta_1_rate,
					double delta_2_shape,
					double delta_2_rate,
					arma::mat Lambda2,
					int times = 1
					) {
	int k = tauh.n_elem;
	arma::mat scores_mat = Lambda_prec % Lambda2;
	int p = Lambda2.n_rows;
	arma::rowvec scores = sum(scores_mat,0);

	double shape, rate;
	arma::vec delta_h;
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
	return(delta.t());
}

// [[Rcpp::export()]]
double log_binom_c(
                      arma::vec beta,
                      arma::mat X,
                      arma::vec y,
                      arma::vec N,
                      arma::vec mu,
                      arma::vec sigma2
) {
  arma::vec Xbeta = X * beta;
  arma::vec e_Xbeta = exp(Xbeta);
  arma::vec p = e_Xbeta / (1.0 + e_Xbeta);
  double ll = sum(y%log(p) + (N-y) % log(1.0-p)) - sum((beta - mu) % (beta - mu) / (2*sigma2));
  return ll;
}


