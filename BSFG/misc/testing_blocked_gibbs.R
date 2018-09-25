# testing the method of Makalik for sampling sigma, beta in a block
library(BSFG)
library(Matrix)
library(MCMCpack)
n = 50
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
chol_R = as(diag(1,n),'CsparseMatrix')
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
chol_R = as(diag(1,n),'CsparseMatrix')
yty = sum(y^2)
xtx = t(X) %*% X
for(i in -100:nI){
  for(j in 1:1){
    nus = 1/rgamma(b,shape = 1,rate = 1+1/lambda2)
    lambda2 = 1/rgamma(b,shape=1,rate=1/nus + beta^2/(2*var))
  }
  if(b<n) {
    res = sample_MME_block(y,X,rep(0,b),1/(lambda2),chol_R,rnorm(b),rep(0,0),rgamma(1,shape=n/2,rate=1))
  } else {
    res = sample_MME_block(y,X,rep(0,b),1/(lambda2),chol_R,rnorm(b),rnorm(n),rgamma(1,shape=n/2,rate=1))
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
