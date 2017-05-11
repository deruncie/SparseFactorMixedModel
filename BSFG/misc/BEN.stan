functions {
	vector sqrt_vec(vector X){
		vector[rows(X)] res;
		for(i in 1:rows(X)){
			res[i] = sqrt(X[i]);
		}
		return(res);
	}
}
data {
  int<lower=1> n;
  int<lower=1> p;
  matrix[n,p] X;
  vector[n] y;
}
parameters {
  real mu;
  vector[p] beta_innov;
  vector<lower=1>[p] tau;
  real<lower=0> sigma;
  real<lower=0> lambda1;
  real<lower=0> lambda2;
}
transformed parameters{
  vector[p] beta_prec = lambda2 / sigma^2 * (tau ./(tau-1.0));
  vector[p] beta = beta_innov ./ sqrt_vec(beta_prec);
}
model {
  sigma ~ cauchy(0,5);

  lambda1 ~ gamma(3,1);
  lambda2 ~ gamma(3,1);

  tau ~ gamma(0.5,lambda1^2/(8.0*lambda2*sigma^2));

  beta_innov ~ normal(0,1);

  y ~ normal(mu + X * beta,sigma);

}
// generated quantities{
//   vector[p] beta;
//   beta = beta_innov ./ sqrt_vec(beta_prec);
// }
