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
}
parameters {
  real logsigma;
  real beta0;
  vector[d] z;
  real<lower=0> tau;
  vector<lower=0>[d] lambda;
  real<lower=0> caux;
}
transformed parameters {
  real<lower=0> sigma;
  vector<lower=0>[d] lambda_tilde;
  real<lower=0> c;
  vector[d] beta;
  vector[n] f;
  sigma = exp(logsigma);
  c = slab_scale * sqrt(caux);
  lambda_tilde = sqrt( c^2 * square(lambda) ./ (c^2 + tau^2*square(lambda)) );
  beta = z .* lambda_tilde*tau;
  f = beta0 + x*beta;
}
model {  // half-t priors for lambdas and tau, and inverse-gamma for c^2
  z ~ normal(0, 1);
  lambda ~ student_t(nu_local , 0,1);
  tau ~ student_t(nu_global , 0, scale_global*sigma);
  caux ~ inv_gamma (0.5* slab_df , 0.5* slab_df);
  beta0 ~ normal(0, scale_icept);
  y ~ normal(f, sigma);
}
