# testing the method of Makalik for sampling sigma, beta in a block
library(BSFG)
library(Matrix)
library(MCMCpack)
library(rstan)
library(rstanarm)
n = 50
b = 100
Var = 1
X = rstdnorm_mat(n,b)
# Beta=rep(0,b)
# Beta = c(1,rep(0,b-1))
Beta = c(rnorm(10),rep(0,b-10))
y = X %*% Beta + rnorm(n,0,sqrt(Var))

p0 = sum(Beta != 0)
tau_0 = p0/(b-p0)/sqrt(n)
tau_0 = 1
options(mc.cores = parallel::detectCores())

colnames(X) = paste0('X',1:ncol(X))
formula = paste(c('y~X1',colnames(X)[-1]),collapse='+')
data = data.frame(y=y,X)
res3 = stan_glm.fit(X,y,prior = hs(global_df = 1,global_scale = tau_0,slab_df = 10e3,slab_scale = 1),prior_intercept = normal(location = 0,scale = 1e6, autoscale = FALSE),
                    prior_aux = normal(1,1e-3,autoscale = FALSE),#cauchy(0,1,autoscale = FALSE))
                    pars = c('global','local','beta','aux','caux'))
# prior_summary(res3)

res = stan('reg_horseshoe_paper.stan',data = list(
  n = n,
  d = b,
  y=c(y),
  x=X,
  scale_icept = 1e6,
  scale_global = tau_0,
  nu_global=1,
  nu_local=1,
  slab_scale=1,
  slab_df=10
))
summary(res)

res2 = stan('reg_horseshoe_paper2.stan',data = list(
  n = n,
  d = b,
  y=c(y),
  x=X,
  scale_icept = 1e6,
  scale_global = tau_0,
  nu_global=1,
  nu_local=1,
  slab_scale=1,
  slab_df=10e6,
  shape_sigma2=1,
  scale_sigma2=1,
  sigma2=1
  # c2 = 1.3,
  # tau = 0.038
))

res4 = stan('reg_horseshoe_paper3.stan',data = list(
  n = n,
  d = b,
  y=c(y),1,
  x=X,
  scale_icept = 1e6,
  scale_global = tau_0,
  nu_global=1,
  nu_local=1,
  slab_scale=1,
  slab_df=10e6,
  shape_sigma2=1,
  scale_sigma2=1,
  sigma2=1,
  lambda2 = lambda2,
  c2 =1,
  beta2 = Beta
),iter=10000)

# , control = list(adapt_delta = 0.9999,max_treedepth=20))

beta = Beta
var = Var
lambda2 = rep(1,b)
nus = rep(1,b)
tau2 = 0.038^2
xi = 1
c2 = 1
c2_shape = 10/2
c2_rate = 10/2


nI = 1000
thin = 10
r1 = matrix(0,nI,b+4)
r2 = matrix(0,nI,b)
chol_R = as(diag(1,n),'dgCMatrix')
beta_prec = 1/(tau2*lambda2)
XtX = crossprod(X)
yty = sum(y^2)
Xty = crossprod(X,y)
for(i in -1000:(thin*nI)){
  # beta_prec = 1/(tau2*lambda2)
  # new_samples = regression_sampler_parallel(
  #   y,
  #   matrix(1,n,1),
  #   list(),
  #   X,
  #   NULL,
  #   1,
  #   list(chol_R),
  #   1/var,
  #   1e10, 1e10,
  #   matrix(1/1e6,nr = 1,ncol = 1),
  #   rep(0,0),
  #   rep(0,b),
  #   beta_prec,
  #   1
  # )
  # var = 1/new_samples$Y_prec
  # beta = new_samples$beta
  # A = XtX
  # diag(A) = diag(A) + beta_prec/var
  # chol_A = chol(A)
  # var=1
  # e = rnorm(b,0,1)
  # beta = backsolve(chol_A,forwardsolve(t(chol_A),Xty/var) + e)

  beta2 = beta^2
  # nus = 1/rgamma(b,shape = 1,rate = 1+1/lambda2)
  # # lambda2 = 1/rgamma(b,shape=1,rate=1/nus + beta2/(2*var))
  # lambda2 = 1/rgamma(b,shape=1,rate=1/nus + beta2/(2*var*tau2))
  tau2 = 1/rgamma(1,shape = (b+1)/2,rate = 1/xi + 1/(2*var)*sum(beta2/lambda2))
  xi = 1/rgamma(1,shape=1,rate=1/(tau_0) + 1/tau2)
  # c2 = 1/rgamma(1,shape = c2_shape + b/2,rate = c2_rate + sum(beta2/var)/2)
  c2=1^2#.3^2
  beta_prec = (c2+tau2*lambda2)/(c2*tau2*lambda2)
  if(i>0 && i %% thin == 0) {
    r1[i/thin,1:b] = beta
    r1[i/thin,b+1] = tau2
    r1[i/thin,b+2] = c2
    r1[i/thin,b+3] = new_samples$alpha1
    r1[i/thin,b+4] = var
    r2[i/thin,] = beta_prec
  }
}

stan_betas = extract(res2,'beta')[[1]]
plot(colMeans(stan_betas),Beta);abline(0,1)
plot(colMeans(r1[,1:b]),Beta);abline(0,1)
plot(colMeans(stan_betas),colMeans(r1[,1:b]));abline(0,1)
plot(apply(stan_betas,2,sd),apply(r1[,1:b],2,sd));abline(0,1)
boxplot(stan_betas)
boxplot(r1[,1:b])

plot(colMeans(1/(r2)),colMeans(1/(extract(res2,'beta_prec')[[1]])),log='xy');abline(0,1)
colSDs = function(x) apply(x,2,sd)
plot(colSDs(1/(r2)),colSDs(1/(extract(res2,'beta_prec')[[1]])),log='xy');abline(0,1)

plot(colMeans(extract(res2,'beta')[[1]]),colMeans(extract(res3,'beta')[[1]]));abline(0,1)
plot(colSDs(extract(res2,'beta')[[1]]),colSDs(extract(res3,'beta')[[1]]));abline(0,1)
plot(colMeans(r1[,1:b]),colMeans(extract(res3,'beta')[[1]]));abline(0,1)
plot(colSDs(r1[,1:b]),colSDs(extract(res3,'beta')[[1]]));abline(0,1)

qqplot(extract(res2,'c2')[[1]],extract(res3,'caux')[[1]]);abline(0,1)
qqplot(log(extract(res2,'tau2')[[1]]),log(r1[,b+1]));abline(0,1)
qqplot((extract(res2,'tau2')[[1]]),(r1[,b+1]));abline(0,1)



plot(colMeans(r1[,1:b]),colMeans(extract(res4,'beta')[[1]]));abline(0,1)
qqplot((extract(res4,'tau2')[[1]]),(r1[,b+1]));abline(0,1)

summary(res2,'c2')$summary
summary(r1[,b+2])
summary(res2,'tau')$summary
summary(sqrt(r1[,b+1]))
summary(res2,'sigma2')$summary
summary(r1[,b+4])
summary(res2,'beta0')$summary
summary(r1[,b+3]);sd(r1[,b+3])

a = extract(res4,'tau_raw')[[1]]
na = length(a)
qqplot(log(a),log(abs(rcauchy(na,0,1))));abline(0,1)
qqplot(log(sqrt(1/rgamma(na,1/2,rgamma(na,1/2,1/2)))),log(abs(rcauchy(length(a),0,1))));abline(0,1)
qqplot(log(sqrt(1/rgamma(na,1/2,rgamma(na,1/2,1)))),log(a));abline(0,1)
qqplot(log(sqrt((1/rgamma(na,1/2,1))/(1/rgamma(na,1/2,1)))),log(a));abline(0,1)

nI = 1000
r1 = matrix(0,nI,b+1)
Sigma_R = crossprod(rstdnorm_mat(n,n)) + diag(1,n)
chol_R = chol(Sigma_R)
# chol_R = as(diag(1,n),'CsparseMatrix')
for(i in -100:(10*nI)){
  nus = 1/rgamma(b,shape = 1,rate = 1+1/lambda2)
  lambda2 = 1/rgamma(b,shape=1,rate=1/nus + beta^2/(2*var))
  if(b<n) {
    beta = sample_MME_single_diagK(y,X,rep(0,b),1/(lambda2*var),chol_R,1/var,rnorm(b),rep(0,0))
  } else{
    beta = sample_MME_single_diagK(y,X,rep(0,b),1/(lambda2*var),chol_R,1/var,rnorm(b),rnorm(n))
  }
  eta = y - X %*% beta
  var = 1/rgamma(1,shape=(n+b)/2,rate = sum(eta^2)/2 + sum(beta^2/lambda2)/2)
  if(i>0 && i %% 10 == 0) r1[i/10,] = c(beta,var)
}



beta = Beta
var = Var
lambda2 = rep(1,b)
nus = rep(1,b)
r2 = matrix(0,nI,b+1)
# chol_R = as(diag(1,n),'CsparseMatrix')
chol_R_ = as(chol_R,'CsparseMatrix')
yty = sum(y^2)
xtx = t(X) %*% X
sX = svd(X)
U = sX$u[,1:50]
V = (sX$d*t(sX$v))[1:50,]
for(i in -100:nI){
  for(j in 1:1){
    nus = 1/rgamma(b,shape = 1,rate = 1+1/lambda2)
    lambda2 = 1/rgamma(b,shape=1,rate=1/nus + beta^2/(2*var))
  }
  if(b<n) {
    res = sample_MME_block(y,X,rep(0,b),1/(lambda2),chol_R_,rnorm(b),rep(0,0),rgamma(1,shape=n/2,rate=1))
  } else {
    set.seed(1)
    rb = rnorm(b)
    if(b>n) {
      re = rnorm(n)
      rg = rgamma(1,shape=n/2,rate=1)
    } else{
      rg = rgamma(1,shape=n/2,rate=1)
      re = rnorm(n)
    }
    res0 = sample_MME_block(y,X,rep(0,b),1/(lambda2),chol_R_,rb,rep(0,0),rg)
    res = sample_MME_block(y,X,rep(0,b),1/(lambda2),chol_R_,rb,re,rg)
    # R = as.matrix(crossprod(chol_R))
    Rinv = chol2inv(chol_R)
    # U=X
    # V = diag(1,b)
    RinvU = Rinv %*% U
    UtRinvU = t(U) %*% RinvU
    W = matrix(0,n,0)
    ra = rnorm(0)
    prec_b = 0
    res2 = regression_sampler_v2(y,W,X,matrix(0,0,0),rep(0,b),1/(lambda2),chol_R,
                                 Sigma_R,1/var,ra,rb,re,rg,prec_b)
    res3 = regression_sampler_v3(y,W,U,V,matrix(0,0,0),rep(0,b),1/(lambda2),chol_R,
                                 Rinv,RinvU,UtRinvU,1/var,ra,rb,re,rg,prec_b)
    RinvSqX = as.matrix(solve(t(chol_R),X))
    res2a = regression_sampler_v2a(y,W,RinvSqX,matrix(0,0,0),rep(0,b),1/(lambda2),chol_R,1/var,ra,rb,re,rg,prec_b)
    C = crossprod(solve(t(chol_R),X))
    res1 = regression_sampler_v1(y,W,RinvSqX,C,matrix(0,0,0),rep(0,b),1/(lambda2),chol_R,1/var,ra,rb,rg,prec_b)

    set.seed(1)
    res4 = regression_sampler_parallel(y,W,NULL,X,NULL,1,list(chol_R),1/var,0,prec_b,rep(0,ncol(W)),rep(0,1),matrix(0,b,1),matrix(1/lambda2,b,1),1)
    set.seed(1)
    res5 = regression_sampler_parallel(y,W,NULL,U,V,1,list(chol_R),1/var,0,prec_b,rep(0,ncol(W)),rep(0,1),matrix(0,b,1),matrix(1/lambda2,b,1),1)
    Wa = matrix(1,n,1)
    res4 = regression_sampler_parallel(y,Wa,NULL,X,NULL,1,list(chol_R),1/var,0,prec_b,rep(0,ncol(Wa)),rep(0,1),matrix(0,b,1),matrix(1/lambda2,b,1),1)
    microbenchmark::microbenchmark(
      # sample_MME_block(y,X,rep(0,b),1/(lambda2),chol_R_,rb,rep(0,0),rg),
      # sample_MME_block(y,X,rep(0,b),1/(lambda2),chol_R_,rb,re,rg),
      regression_sampler_v2(y,W,X,matrix(0,0,0),rep(0,b),1/(lambda2),chol_R,
      Sigma_R,1/var,ra,rb,re,rg,prec_b),
      regression_sampler_v3(y,W,U,V,matrix(0,0,0),rep(0,b),1/(lambda2),chol_R,
                            Sigma_R,Rinv,RinvU,UtRinvU,1/var,ra,rb,re,rg,prec_b),
      regression_sampler_v2a(y,W,RinvSqX,matrix(0,0,0),rep(0,b),1/(lambda2),chol_R,1/var,ra,rb,re,rg,prec_b),
      regression_sampler_v1(y,W,RinvSqX,C,matrix(0,0,0),rep(0,b),1/(lambda2),chol_R,1/var,ra,rb,rg,prec_b)
    )
  }
  var = 1/res[1]
  beta = res[-1]
  # Ainv = solve(t(X) %*% X + diag(1/lambda2))
  # XAinvXt = X %*% Ainv %*% t(X)
  # rate = (yty - t(y) %*% (XAinvXt) %*% y)/2
  # var = 1/rgamma(1,shape=(n)/2,rate = rate)
  # beta = sample_MME_single_diagK(y,X,rep(0,b),1/(lambda2*var),chol_R,1/var,rnorm(b),rep(0,0))
  if(i>0) r2[i,] = c(beta,var)
}
apply(r1,2,effectiveSize)
apply(r2,2,effectiveSize)
plot(colMeans(r1),colMeans(r2));abline(0,1)


# plot(as.mcmc(r2))
