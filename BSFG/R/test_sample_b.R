n = 100
p = 200
s2 = 1/2

R = crossprod(rstdnorm_mat(n,n)) + diag(1,n)
chol_R = as(chol(R),'dgCMatrix')
X = matrix(sample(c(0,1),n*p,replace = T),n,p)

b = rnorm(p)

y = X %*% b + t(chol_R) %**% rnorm(n) * sqrt(s2)


randn_theta = rnorm(p)
randn_e = rnorm(n)
prior_mean = rep(0,p)
prior_prec = rep(1,p)
tot_Eta_prec = 1/s2

b1 <- sample_MME_single_diagK2(y,X,prior_mean,prior_prec,chol_R,tot_Eta_prec,randn_theta,matrix(0,0,0))
b2 <- sample_MME_single_diagK2(y,X,prior_mean,prior_prec,chol_R,tot_Eta_prec,randn_theta,randn_e)

b1s = sapply(1:1000,function(x) {
  randn_theta = rnorm(p)
  sample_MME_single_diagK2(y,X,prior_mean,prior_prec,chol_R,tot_Eta_prec,randn_theta,matrix(0,0,0))
})
b2s = sapply(1:1000,function(x) {
  randn_theta = rnorm(p)
  randn_e = rnorm(n)
  sample_MME_single_diagK2(y,X,prior_mean,prior_prec,chol_R,tot_Eta_prec,randn_theta,randn_e)
})
