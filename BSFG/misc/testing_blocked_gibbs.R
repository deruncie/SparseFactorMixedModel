# testing the method of Makalik for sampling sigma, beta in a block
library(BSFG)
library(Matrix)
library(MCMCpack)
n = 300
b = 60
Var = 1
X = rstdnorm_mat(n,b)
Beta = c(1,rep(0,b-1))
y = X %*% Beta + rnorm(n,0,sqrt(Var))

beta = Beta
var = Var
lambda2 = rep(1,b)
nus = rep(1,b)

nI = 1000
r1 = matrix(0,nI,b+1)
Sigma_R = crossprod(rstdnorm_mat(n,n)) + diag(1,n)
chol_R = chol(Sigma_R)
# chol_R = as(diag(1,n),'CsparseMatrix')
for(i in -100:nI){
  nus = 1/rgamma(b,shape = 1,rate = 1+1/lambda2)
  lambda2 = 1/rgamma(b,shape=1,rate=1/nus + beta^2/(2*var))
  if(b<n) {
    beta = sample_MME_single_diagK(y,X,rep(0,b),1/(lambda2*var),chol_R,1/var,rnorm(b),rep(0,0))
  } else{
    beta = sample_MME_single_diagK(y,X,rep(0,b),1/(lambda2*var),chol_R,1/var,rnorm(b),rnorm(n))
  }
  eta = y - X %*% beta
  var = 1/rgamma(1,shape=(n+b)/2,rate = sum(eta^2)/2 + sum(beta^2/lambda2)/2)
  if(i>0) r1[i,] = c(beta,var)
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
    } else{
      re = rnorm(0)
    }
    rg = rgamma(1,shape=n/2,rate=1)
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
                                 Sigma_R,Rinv,RinvU,UtRinvU,1/var,ra,rb,re,rg,prec_b)
    RinvSqX = as.matrix(solve(t(chol_R),X))
    res2a = regression_sampler_v2a(y,W,RinvSqX,matrix(0,0,0),rep(0,b),1/(lambda2),chol_R,1/var,ra,rb,re,rg,prec_b)
    C = crossprod(solve(t(chol_R),X))
    res1 = regression_sampler_v1(y,W,RinvSqX,C,matrix(0,0,0),rep(0,b),1/(lambda2),chol_R,1/var,ra,rb,rg,prec_b)

    set.seed(1)
    res4 = regression_sampler_parallel(y,W,NULL,X,NULL,1,list(chol_R),1/var,0,prec_b,rep(0,ncol(W)),rep(0,1),matrix(0,b,1),matrix(1/lambda2,b,1),1)
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
