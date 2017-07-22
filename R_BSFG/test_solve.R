library(BSFG)
library(Matrix)

n = 100
b = 20000
p = 5
X = matrix(rnorm(n*b),n)
X = matrix(0,n,b)
diag(X)=1
bs = matrix(rnorm(b*p),b)
Y = .1*matrix(rnorm(n*p),n) + X %*% bs

X1 = X[,1:n]
b1 = solve(crossprod(X1)) %*% crossprod(X1,Y)


Ut = as(diag(1,n),'dgCMatrix')
h2 = rep(0,p)
s = rep(1,n)

prior_mean = matrix(0,b,p)
prior_prec = matrix(1e-10,b,p)
randn_theta = matrix(rnorm(b*p),b)
randn_e = matrix(rnorm(n*p),n)
tot_Y_prec = matrix(1,p)

W = as(X,'dgCMatrix')

b2=sample_coefs_parallel_sparse_c_Eigen(Ut,Y,W,h2,tot_Y_prec,s,prior_mean,prior_prec,randn_theta, randn_e,1)
