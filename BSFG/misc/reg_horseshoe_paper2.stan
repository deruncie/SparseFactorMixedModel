data {
  int<lower=0> n;
  int<lower=0> d;
  vector[n] y;
  matrix[n,d] x;
  real<lower=0> scale_icept;
  real<lower=0> scale_global;
  real<lower=1> nu_global;
  real<lower=1> nu_local;
  real<lower=0> slab_scale;
  real<lower=0> slab_df;
  real<lower=1> shape_sigma2;
  real<lower=0> scale_sigma2;
  // real<lower=0> c2;
  // real<lower=0> tau;
  real<lower=0> sigma2;
}
parameters {
  // real<lower=0> sigma2;
  // real<lower=0> sigma3;
  real beta0;
  vector[d] z;
  real<lower=0> tau;
  vector<lower=0>[d] lambda;
  // vector<lower=0>[d] nu;
  // vector<lower=0>[d] lambda3;
  real<lower=0> caux;
}
transformed parameters {
  // vector<lower=0>[d] lambda_tilde;
  vector<lower=0>[d] lambda2;
  // vector<lower=0>[d] lambda32;
  vector<lower=0>[d] beta_prec;
  real<lower=0> c2;
  real<lower=0> tau2;
  vector[d] beta;
  vector[n] f;
  // sigma = exp(logsigma);
  c2 = slab_scale^2 * caux;
  lambda2 = square(lambda);
  tau2 = tau^2;
  // lambda32 = square(lambda3);
  beta_prec = (c2 + (tau2 * lambda2)) ./ (c2 * tau2 * lambda2);
  // beta_prec = 1.0 ./ (tau^2 * lambda2 ./ nu);
  beta = sqrt(sigma2) * (z ./ sqrt(beta_prec));
  // beta = sqrt(tau)*sqrt(sigma2) * (z .* sqrt(lambda2 ./ nu));
  // lambda_tilde = sqrt( c2 * tau^2*square(lambda) ./ (c2 + tau^2*square(lambda)) );
  // beta = z .* lambda_tilde*sqrt(sigma2);
  f = beta0 + x*beta;
  // f = x*beta;
}
model {  // half-t priors for lambdas and tau, and inverse-gamma for c^2
  // sigma2 ~ inv_gamma(shape_sigma2,scale_sigma2);
  // sigma3 ~ inv_gamma(shape_sigma2,scale_sigma2);
  z ~ normal(0, 1);
  // lambda ~ student_t(nu_local , 0,1);
  // tau ~ student_t(nu_global , 0, scale_global);
  lambda ~ cauchy(0,1);
  // lambda2 ~ inv_gamma(1.0/2.0,1.0/2.0);
  // nu ~ inv_gamma(1.0/2.0,1.0);
  // lambda3 ~ cauchy(0,1);
  tau ~ cauchy(0,scale_global);
  caux ~ inv_gamma (0.5* slab_df , 0.5* slab_df);
  beta0 ~ normal(0, scale_icept);
  y ~ normal(f, sqrt(sigma2));
}
