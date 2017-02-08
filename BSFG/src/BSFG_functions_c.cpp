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


// [[Rcpp::export()]]
arma::mat sample_coefs_parallel_sparse_c(
					arma::mat Y,
					arma::mat W,
					arma::vec h2,
					arma::vec tot_Y_prec,
					arma::mat prior_mean,
					arma::mat prior_prec,
					List invert_aI_bZAZ,
					int grainSize) {

	// Sample regression coefficients
	// columns of matrices are independent
	// each column conditional posterior is a MVN due to conjugacy


	struct sampleColumn : public Worker {
		arma::mat WtU, UtY, prior_prec, Z;
		arma::vec h2, tot_Y_prec, s;
		int b;
		arma::mat &coefs;

		sampleColumn(arma::mat WtU, arma::mat UtY, arma::mat prior_prec, arma::mat Z, arma::vec h2, arma::vec tot_Y_prec, arma::vec s, int b, arma::mat &coefs) :
			WtU(WtU), UtY(UtY), prior_prec(prior_prec), Z(Z), h2(h2), tot_Y_prec(tot_Y_prec), s(s), b(b), coefs(coefs) {}

      	void operator()(std::size_t begin, std::size_t end) {
			arma::mat WtUDi, Q, cholQ;
			arma::vec means,v,m,y,z;
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

	arma::sp_mat U = as<arma::sp_mat>(invert_aI_bZAZ["U"]);
	arma::vec s = as<arma::vec>(invert_aI_bZAZ["s"]);

	arma::mat WtU = W.t() * U;
	arma::mat UtY = U.t() * Y;

	arma::mat Z = randn(b,p);
	arma::mat coefs = zeros(b,p);

	sampleColumn sampler(WtU, UtY, prior_prec, Z, h2, tot_Y_prec, s, b, coefs);
	RcppParallel::parallelFor(0,p,sampler,grainSize);
	return(coefs);
}

// [[Rcpp::export()]]
arma::rowvec sample_tot_prec_sparse_c (arma::mat Y,
					   arma::vec h2,
					   double tot_Y_prec_shape,
					   double tot_Y_prec_rate,
					   List invert_aI_bZAZ
					  ) {

	arma::sp_mat U = as<arma::sp_mat>(invert_aI_bZAZ["U"]);
	arma::vec s = as<arma::vec>(invert_aI_bZAZ["s"]);

	arma::mat UtY = U.t() * Y;

	int n = Y.n_rows;
	int p = Y.n_cols;

	arma::vec tot_Y_prec = zeros(p);

	for(int i = 0; i < p; i++){
		arma::vec Sigma_sqrt = sqrt(h2(i) * s + (1.0 - h2(i)));
		arma::vec SiUtY_i = UtY.col(i) / Sigma_sqrt;
		arma::vec prec = randg(1,distr_param(tot_Y_prec_shape + n/2, 1.0/(tot_Y_prec_rate + 0.5 * dot(SiUtY_i,SiUtY_i))));
		tot_Y_prec(i) = prec(0);
	}
	return(tot_Y_prec.t());
}

// [[Rcpp::export()]]
arma::rowvec sample_h2s_discrete_given_p_sparse_c (arma::mat Y,
						int h2_divisions,
						arma::vec h2_priors,
						arma::vec Tot_prec,
						List invert_aI_bZAZ){

	arma::sp_mat U = as<arma::sp_mat>(invert_aI_bZAZ["U"]);
	arma::vec s = as<arma::vec>(invert_aI_bZAZ["s"]);

	int p = Y.n_cols;
	int n = Y.n_rows;
	arma::vec h2_index = zeros(p);

	arma::mat log_ps = zeros(p,h2_divisions);
	arma::mat std_scores_b = sweep_times(Y.t() * U,1,sqrt(Tot_prec));

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
arma::mat sample_randomEffects_parallel_sparse_c (arma::mat Y,
				arma::sp_mat Z,
				arma::vec tot_prec,
				arma::vec h2,
				List invert_aZZt_Ainv,
				int grainSize ) {
	//samples genetic effects on factors (F_a) conditional on the factor scores F:
	// F_i = F_{a_i} + E_i, E_i~N(0,s2*(1-h2)*I) for each latent trait i
	// U_i = zeros(r,1) if h2_i = 0
	// it is assumed that s2 = 1 because this scaling factor is absorbed in
	// Lambda
	// invert_aZZt_Ainv has parameters to diagonalize a*Z*Z' + b*I for fast
	// inversion:

	arma::vec a_prec = tot_prec / h2;
	arma::vec e_prec = tot_prec / (1-h2);

	arma::sp_mat U = as<arma::sp_mat>(invert_aZZt_Ainv["U"]);
	arma::vec s1 = as<arma::vec>(invert_aZZt_Ainv["s1"]);
	arma::vec s2 = as<arma::vec>(invert_aZZt_Ainv["s2"]);

	int p = Y.n_cols;
	int r = Z.n_cols;
	arma::mat b = U.t() * Z.t() * sweep_times(Y,2,e_prec);

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
arma::mat sample_means_parallel_c(arma::mat Y_tilde,
				   arma::vec resid_Y_prec,
				   arma::vec E_a_prec,
				   List invert_aA_bZtZ,
				   int grainSize ) {
	// when used to sample [B;E_a]:
	//  W - F*Lambda' = X*B + Z_1*E_a + E, arma::vec(E)~N(0,kron(Psi_E,In)).
	//  Note: conditioning on F, Lambda and W.
	// The vector [b_j;E_{a_j}] is sampled simultaneously. Each trait is sampled separately because their
	// conditional posteriors factor into independent MVNs.
	// note:invert_aA_bZtZ has parameters to diagonalize mixed model equations for fast inversion:
	// inv(a*blkdiag(fixed_effects_prec*eye(b),Ainv) + b*[X Z_1]'[X Z_1]) = U*diag(1./(a.*s1+b.*s2))*U'
	// Z_U = [X Z_1]*U, which doesn't change each iteration.

	struct sampleColumn : public Worker {
		arma::vec E_a_prec, resid_Y_prec, s1, s2;
		arma::mat means, Zlams, U;

		arma::mat &location_sample;

		sampleColumn(arma::vec E_a_prec, arma::vec resid_Y_prec, arma::vec s1, arma::vec s2, arma::mat means, arma::mat Zlams, arma::mat U, arma::mat &location_sample)
			: E_a_prec(E_a_prec), resid_Y_prec(resid_Y_prec), s1(s1), s2(s2), means(means), Zlams(Zlams),U(U), location_sample(location_sample) {}

      	void operator()(std::size_t begin, std::size_t end) {
			arma::vec d, mlam;
			for(std::size_t j = begin; j < end; j++){
				d = s1*E_a_prec(j) + s2*resid_Y_prec(j);
				mlam = means.col(j) /d;
				location_sample.col(j) = U * (mlam + Zlams.col(j)/sqrt(d));
			}
		}
	};

	arma::mat U = as<arma::mat>(invert_aA_bZtZ["U"]);
	arma::vec s1 = as<arma::vec>(invert_aA_bZtZ["s1"]);
	arma::vec s2 = as<arma::vec>(invert_aA_bZtZ["s2"]);
	arma::mat Z_U = as<arma::mat>(invert_aA_bZtZ["Z_U"]);

	// int n = Y_tilde.n_rows;
	int p = Y_tilde.n_cols;
	int br = Z_U.n_cols;

	arma::mat means = sweep_times(Z_U.t() * Y_tilde,2,resid_Y_prec);
	arma::mat location_sample = zeros(br,p);

	arma::mat Zlams = randn(br,p);

	sampleColumn sampler(E_a_prec,resid_Y_prec,s1,s2,means,Zlams,U,location_sample);
	RcppParallel::parallelFor(0,p,sampler,grainSize);

	return(location_sample);
}

// [[Rcpp::export()]]
arma::mat sample_factors_scores_sparse_c(arma::mat Y_tilde,
                                         arma::mat prior_mean,
                                         arma::mat Lambda,
                                         arma::vec resid_Y_prec,
                                         arma::vec F_e_prec
) {
  //Sample factor scores given factor loadings (F_a), factor residual variances (F_e_prec) and
  //phenotype residuals
  arma::mat Lmsg = sweep_times(Lambda,1,resid_Y_prec);
  arma::mat S = chol(Lambda.t() * Lmsg + diagmat(F_e_prec));
  arma::mat tS = S.t();

  arma::mat Meta = trans(solve(tS,trans(Y_tilde * Lmsg + sweep_times(prior_mean,2,F_e_prec))));

  arma::mat Zlams = randn(Meta.n_rows,Meta.n_cols);

  arma::mat F = trans(solve(S,trans(Meta + Zlams)));

  return(F);
}

// [[Rcpp::export()]]
List GSVD_2_c(arma::mat A, arma::mat B){
	int n = B.n_cols;
	arma::mat U,V;
	arma::vec d;
	svd(U,d,V,A * inv(B));

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

