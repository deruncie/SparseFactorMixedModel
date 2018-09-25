# testing regularized horseshoe
# testing if I can sample it with a Gibbs sampler
# Gibbs sampler from Makalic and Schmidt (2015)
# regularization from Piironen and Vehtari (2017)

# also want to check induced prior on var(X\beta)


n = 1e6
tau_0 = 1

# direct horseshoe
lambda = abs(rcauchy(n))
# tau = abs(rcauchy(n,0,tau_0))
tau = 0.1
vars_1 = lambda^2*tau^2
beta_1 = rnorm(n,0,sqrt(vars_1))

# Makalic and Schmidt algorithm
nu = 1/rgamma(n,shape=1/2,rate=1)
lambda2 = 1/rgamma(n,shape=1/2,rate=1/nu)
xi = 1/rgamma(n,shape=1/2,rate=1)
# tau2 = 1/rgamma(n,shape=1/2,rate=1/xi)
tau2 = tau^2
vars_2 = lambda2*tau2
beta_2 = rnorm(n,0,sqrt(vars_2))

# qqplot(beta_1,beta_2);abline(0,1)


k = 100;plot(quantile(vars_1,seq(0,1,length=k)),quantile(vars_2,seq(0,1,length=k)),log='xy');abline(0,1)


# direct regularized horseshoe
lambda = abs(rcauchy(n))
c = .1
tau = 0.1
lambda2_tilde = c^2*lambda^2/(c^2 + tau^2*lambda^2)
vars_1r = tau^2*lambda2_tilde
beta_1r = rnorm(n,0,sqrt(vars_1r))
k = 100;plot(quantile(vars_1,seq(0,1,length=k)),quantile(vars_1r,seq(0,1,length=k)),log='xy');abline(0,1)

# indirect regularized horseshoe
nu = 1/rgamma(n,shape=1/2,rate=1)
lambda2 = 1/rgamma(n,shape=1/2,rate=1/nu)
# xi = 1/rgamma(n,shape=1/2,rate=1)
# tau2 = 1/rgamma(n,shape=1/2,rate=1/xi)
tau2 = tau^2
c2 = c^2
lambda2_tilde = c^2*lambda2/(c^2 + tau2*lambda2)
vars_2r = tau2*lambda2_tilde
beta_2r = rnorm(n,0,sqrt(vars_2r))
k = 100;plot(quantile(vars_1r,seq(0,1,length=k)),quantile(vars_2r,seq(0,1,length=k)),log='xy');abline(0,1)


# regularized horseshoe total variance
# total variance is sum(beta_j^2), j = 1:p
n = 1000
n_c = 10
n_r = 2000
p = 1000
cs = matrix(rep(seq(0,1,length=n_c+1)[-1],each = n_r),p,n_c*n_r,byrow=T)
p_0 = 9
tau_0 = p_0 / (p - p_0) * 1 / sqrt(n)
tau = matrix(abs(rcauchy(n_c*n_r,0,tau_0)),p,n_c*n_r,byrow=T)
# tau = matrix(abs(rnorm(n_c*n_r,0,tau_0)),p,n_c*n_r,byrow=T)
# tau = tau_0
# tau = matrix(sapply(9^seq(3,1,length=n_c),function(x) x/(p-x)),p,n_c,byrow=T)
lambda = matrix(abs(rcauchy(n_c*n_r*p)),ncol = n_c*n_r)
# b = 1/(1+n*cs^2)
# k = 1/(1+n*tau^2*lambda^2);k = (1-b)*k+b
# plot(colMeans(1-b)*p_0);abline(h=p_0)
# plot(colSums(1-k));abline(h=p_0)
lambda2_tilde = cs^2*lambda^2/(cs^2 + tau^2*lambda^2)
# k = 1/(1+n*tau^2*lambda2_tilde)
# plot(colSums(1-k));abline(h=p_0)
# hist(colSums(1-k),breaks=c(seq(0,20,length=100),1e6),xlim = c(0,20))
# lambda2_tilde = lambda^2
vars = tau^2*lambda2_tilde
beta = matrix(rnorm(n_c*p,0,sqrt(vars)),nrow = p)
# plot(cs[1,],colSums(1-k));abline(h=p_0)
# plot(cs[1,],colSums(beta^2));lines(cs[1,],cs[1,]*p_0)
# vars[vars<.5*cs^2]=0
# vars[vars > cs^2] = NA
# plot(cs[1,],colSums(vars));lines(cs[1,],cs[1,]^2*p_0)
plot(jitter(cs[1,]),log(colSums(vars)));lines(cs[1,],log(cs[1,]^2*p_0))
# plot(cs[1,seq(1,n_c*n_r,by=n_r)],colMeans(matrix(colSums(vars),n_r)>.1*cs[1,seq(1,n_c*n_r,by=n_r)]^2*p_0))
boxplot((matrix(colSums(vars),n_r)))
