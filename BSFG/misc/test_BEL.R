library(Matrix)
library(statmod)
library(BSFG)

n = 200
p = 500
b = 20

# X = cbind(1,matrix(sample(c(0,1),n*p,replace=T),n,p))
# X = cbind(1,X)
n = nrow(X)
p = ncol(X)-1
X = sweep(X,2,colMeans(X),'-')
b_act = c(0,rep(0,p))
b_act[sample(1:p,b)] = rnorm(b)
# b_act[1:b] = rnorm(b)

y = X %*% b_act + .1*rnorm(n)


l1 = 1
l2 = 10

b = b_act+1e-10
prec_y = 10

Ut = as(diag(1,n),'dgCMatrix')
s = rep(1,n)
h2 = 0
prior_mean = rep(0,p+1)

bs = c()
prec_ys = c()
prec_bs = c()
l1s = c()
l2s = c()

t = 1

log_post = function(log_l2,tau,l1,b,prec_y,t,shape,rate){
  l2 = exp(log_l2)
  (shape-1 + length(b)/2)*log_l2 - 1/2*sum(tau/(tau-1)*b^2)*l2*t*prec_y - length(b)/2*log_l2 - prec_y*t*l1^2*sum(tau)/(8*l2) - l2*rate
}

for(i in 1:1000){
  if(i %% 100 == 0) print(i)
  # inverse Gaussian Bayesian Elastic Net from Li and Lin 2010
  # x  = pmax(1e-10,SuppDists::rinvGauss(n = p, nu = sqrt(l1)/(2*l2*abs(b)), lambda =  prec_y*t*l1/(4*l2)))  # wrong I think
  x  = pmax(1e-10,SuppDists::rinvGauss(n = p, nu = l1/(2*l2*abs(b)), lambda =  prec_y*t*l1^2/(4*l2)))
  tau = 1/x+1
  prec_b = c(1e-1,tau/(tau-1)*l2*t*prec_y)

  randn_theta = rnorm(p+1)
  randn_e = rnorm(n)
  b = sample_coefs_parallel_sparse_c_Eigen(Ut,y,X,h2, prec_y,s, prior_mean,prec_b,randn_theta,randn_e,1)

  l1 = sqrt(rgamma(1,shape = 3+p/2, rate = 1+prec_y*t*sum(tau)/(8*l2)))

  old_state = log(l2)
  trial_state = old_state + rnorm(1,0,.5)
  logps = log_post(c(old_state,trial_state),tau,l1,b[-1],prec_y,t,shape=3,rate=1)
  if(log(runif(1)) < diff(logps)) l2 = exp(trial_state)

  scores = tot_prec_scores_c(y - X %*% b,h2,s)
  alpha = n/2 + p
  beta = 1/2*(scores + l2*sum(tau/(tau-1)*b[-1]^2*t) + l1^2/(4*l2)*sum(tau))
  var_y = rinvgamma(1,shape = alpha,scale = beta)
  if(log(runif(1)) < p*log(gamma(0.5)) - p*log(pracma::gammainc(x=l1^2/(8*var_y*l2),a=0.5)['uppinc'])) {
    prec_y = 1/var_y
  }
  # prec_y = rgamma(1,shape = 2 + n/2,rate = 1.5 + 0.5*scores)

  bs = rbind(bs,t(b))
  prec_bs = rbind(prec_bs,prec_b)
  prec_ys = c(prec_ys,prec_y)
  l1s = c(l1s,l1)
  l2s = c(l2s,l2)

}

trace_plot(bs[,-1])
trace_plot(cbind(l1s,l2s))
posterior_plot(bs[,-1])
plot(colMeans(bs[,-1]),b_act[-1]);abline(0,1)


prec_y = 10

bs = c()
prec_ys = c()
prec_bs = c()
l1s = c()
l2s = c()

t = 1

log_post = function(log_l2,tau,l1,b,shape,rate){
  l2 = exp(log_l2)
  (shape-1 + length(b)/2)*log_l2 - 1/2*sum(tau/(tau-1)*b^2)*l2*t - length(b)/2*log_l2 - t*l1^2*sum(tau)/(8*l2) - l2*rate
}

for(i in 1:1000){
  if(i %% 100 == 0) print(i)
  # inverse Gaussian Bayesian Elastic Net from Li and Lin 2010
  # x  = pmax(1e-10,SuppDists::rinvGauss(n = p, nu = sqrt(l1)/(2*l2*abs(b)), lambda =  t*l1/(4*l2)))
  x  = pmax(1e-10,SuppDists::rinvGauss(n = p, nu = l1/(2*l2*abs(b)), lambda =  t*l1^2/(4*l2)))
  tau = 1/x+1
  prec_b = c(1e-1,tau/(tau-1)*l2*t)

  randn_theta = rnorm(p+1)
  randn_e = rnorm(n)
  b = sample_coefs_parallel_sparse_c_Eigen(Ut,y,X,h2, prec_y,s, prior_mean,prec_b,randn_theta,randn_e,1)

  l1 = sqrt(rgamma(1,shape = 3+p/2, rate = 1+t*sum(tau)/(8*l2)))

  old_state = log(l2)
  trial_state = old_state + rnorm(1)
  logps = log_post(c(old_state,trial_state),tau,l1,b[-1],shape=3,rate=1)
  if(log(runif(1)) < diff(logps)) l2 = exp(trial_state)

  scores = tot_prec_scores_c(y - X %*% b,h2,s)
  prec_y = rgamma(1,shape = 2 + n/2,rate = 1.5 + 0.5*scores)

  bs = rbind(bs,t(b))
  prec_bs = rbind(prec_bs,prec_b)
  prec_ys = c(prec_ys,prec_y)
  l1s = c(l1s,l1)
  l2s = c(l2s,l2)

}

trace_plot(bs[,-1])
trace_plot(cbind(l1s,l2s))
posterior_plot(bs[,-1])
plot(colMeans(bs[,-1]),b_act[-1]);abline(0,1)



bs = matrix(0,1000,length(b_act))
prec_ys = c()
prec_bs = c()
l1s = c()
l2s = c()

t = 1000
df=1
for(i in 1:1000){
  if(i %% 100 == 0) print(i)
  # inverse Gaussian Bayesian Elastic Net from Li and Lin 2010
  prec_b = c(1e-1,rgamma(p,shape = (df+1)/2,rate = (df+b^2*t)/2))
  prec_b = prec_b*t

  randn_theta = rnorm(p+1)
  randn_e = rnorm(n)
  b = sample_coefs_parallel_sparse_c_Eigen(Ut,y,X,h2, prec_y,s, prior_mean,prec_b,randn_theta,randn_e,1)

  scores = tot_prec_scores_c(y - X %*% b,h2,s)
  prec_y = rgamma(1,shape = 2 + n/2,rate = 1.5 + 0.5*scores)

  bs[i,] = c(b)
  # bs = rbind(bs,t(b))
  # prec_bs = rbind(prec_bs,prec_b)
  prec_ys = c(prec_ys,prec_y)
}

trace_plot(bs[,-1])
p1=posterior_plot(bs[,-1]);p1 + geom_point(data = data.frame(x=1:p,y=b_act[-1]),aes(x=x,y=y))#geom_vline(xintercept = which(b_act[-1] != 0))
plot(colMeans(bs[,-1]),b_act[-1]);abline(0,1)


## Bayesian Lasso from BLR
lambda = 50
bs = matrix(0,1000,length(b_act))
prec_ys = c()
prec_bs = c()
lambdas = c()

t = 1
# lambda = 3.2
df=1
for(i in 1:1000){
  if(i %% 100 == 0) print(i)
  # inverse Gaussian Bayesian Elastic Net from Li and Lin 2010

  b_std = b*sqrt(t)
  nu <- lambda/abs(b_std[-1]) / sqrt(prec_b[-1])
  tmp  = SuppDists::rinvGauss(n = p, nu = nu, lambda = lambda^2)
  tau_j = 1/tmp
  prec_b[-1] <- prec_y/tau_j
  lambda = sqrt(rgamma(1,shape = 3+p, rate = 1+sum(tau_j)/2))

  randn_theta = rnorm(p+1)
  randn_e = rnorm(n)
  b = sample_coefs_parallel_sparse_c_Eigen(Ut,y,X,h2, prec_y,s, prior_mean,prec_b*t,randn_theta,randn_e,1)

  scores = tot_prec_scores_c(y - X %*% b,h2,s)
  prec_y = rgamma(1,shape = 2 + n/2 + p/2,rate = 1.5 + 0.5*(scores + sum(b[-1]^2/tau_j)))

  bs[i,] = c(b)
  lambdas = c(lambdas,lambda)
  # bs = rbind(bs,t(b))
  # prec_bs = rbind(prec_bs,prec_b)
  prec_ys = c(prec_ys,prec_y)
}

trace_plot(bs[,-1])
p1=posterior_plot(bs[,-1]);p1 + geom_point(data = data.frame(x=1:p,y=b_act[-1]),aes(x=x,y=y)) + ylim(c(-3,2))#geom_vline(xintercept = which(b_act[-1] != 0))
plot(colMeans(bs[,-1]),b_act[-1]);abline(0,1)



library(glmnet)
res = cv.glmnet(X,y,alpha=0)
plot(b_act,res$glmnet.fit$beta[,order(res$cvm)[1]])
res = glmnet(X,y,alpha=l1)
plot(res$beta,colMeans(bs));abline(0,1)

library(EBEN)
res2 = EBelasticNet.Gaussian(X[,-1],y,lambda = .24,alpha = .5)
plot(b_act[-1][res2$weight[,1]],res2$weight[,3])

library(BLR)
res3 = BLR(y,XL = X,prior=list(varE=list(df=3,S=0.25),
                               varU=list(df=3,S=0.63),
                               lambda=list(shape=0.52,rate=1e-4,
                                           type='random',value=30)),nIter=5500,burnIn=500,thin=1)


grav2 = grav
grav2$pheno = y
out = scanone(grav2,pheno.col = 1,method='hk')
plot(out$lod)

