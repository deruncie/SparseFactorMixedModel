library(Matrix)
library(BSFG)
library(rbenchmark)
n = 1000
p = 2000
s2 = 40

# R = tcrossprod(rstdnorm_mat(n,3)) + diag(1,n)
R = diag(1,n)
chol_R = drop0(as(chol(R),'dgCMatrix'),tol=1e-10)
inv_chol_Rt = drop0(as(solve(t(chol_R)),'dgCMatrix'),tol=1e-10)

str(chol_R)
str(inv_chol_Rt)

X1 = rstdnorm_mat(n,p)
# X = matrix(sample(c(0,1),n*p,replace = T),n,p)

# X = sweep(X,2,colMeans(X),'-')
b = c(rep(4,5),2^((6-6:15)/2),rep(0,p-15))

y = X1 %*% b + rnorm(n,0,sqrt(s2))#t(chol_R) %**% rnorm(n) * sqrt(s2)
X = X1#.7*X1 + .3*rstdnorm_mat(n,p)


p = ncol(X)
n = nrow(X)
X = sweep(X,2,colMeans(X),'-')
b = sample(c(rep(4,5),0*2^((6-6:15)/2),rep(0,p-15)))
y = X %*% b + rnorm(n,0,sqrt(s2))

# y = y - mean(y)
#
#
# randn_theta = rnorm(p)
# randn_e = rnorm(n)
# prior_mean = rep(1,p)
# prior_prec = rep(1,p)
# tot_Eta_prec = 1/s2
#
# benchmark(
#   b1 <- sample_MME_single_diagK2(y,X,prior_mean,prior_prec,chol_R,tot_Eta_prec,randn_theta,matrix(0,0,0)),
#   b1b <- sample_MME_single_diagK2b(y,X,prior_mean,prior_prec,inv_chol_Rt,tot_Eta_prec,randn_theta,matrix(0,0,0)),
# # b2 <- sample_MME_single_diagK3(y,X,prior_mean,prior_prec,chol_R,tot_Eta_prec,randn_theta,randn_e),
# # b3 <- sample_MME_single_diagK2(y,X,prior_mean,prior_prec,chol_R,tot_Eta_prec,randn_theta,randn_e),
# # b3b <- sample_MME_single_diagK2b(y,X,prior_mean,prior_prec,inv_chol_Rt,tot_Eta_prec,randn_theta,randn_e),
# replications = 100
# )
#
# b1s = sapply(1:1000,function(x) {
#   randn_theta = rnorm(p)
#   sample_MME_single_diagK3(y,X,prior_mean,prior_prec,chol_R,tot_Eta_prec,randn_theta,matrix(0,0,0))
# })
# b2s = sapply(1:1000,function(x) {
#   randn_theta = rnorm(p)
#   randn_e = rnorm(n)
#   sample_MME_single_diagK3(y,X,prior_mean,prior_prec,chol_R,tot_Eta_prec,randn_theta,randn_e)
# })
# b3s = sapply(1:1000,function(x) {
#   randn_theta = rnorm(p)
#   randn_e = rnorm(n)
#   sample_MME_single_diagK2(y,X,prior_mean,prior_prec,chol_R,tot_Eta_prec,randn_theta,randn_e)
# })

# y = read.csv('y.csv',h=F)[,1]
# X = as.matrix(read.csv('X.csv',h=F))
# p = ncol(X)
res = horseshoe2(y,X,1/2,1/2,2000,00,T)
plot(colMeans(res[,1:p]));abline(v=which(b>0))

res2 = horseshoe1(y,X,1/2,1/2,2000,00,1,T)
plot(colMeans(res2[,1:p]));abline(v=which(b>0))

res3 = ard(y,X,1/4,1/4,3/2,1,1,1,2000,100,T)
plot(colMeans(res3[,1:p]));abline(v=which(b>0))

hist(res[,ncol(res)]);abline(v=s2)
hist(res2[,ncol(res2)]);abline(v=s2)
hist(res3[,ncol(res3)]);abline(v=s2)

posterior_plot(res[,1:p])
posterior_plot(res2[,1:p])
posterior_plot(res3[,1:p])
#
# res2 = blasso(X = X,y = y,T = 1000,thin = 1,case='hs',RJ = F)
# res2 = bhs(X = X,y=y)
# plot(colMeans(res2$beta[-c(1:200),]))
