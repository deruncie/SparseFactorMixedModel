library(Matrix)
library(statmod)
library(BSFG)
library(MCMCpack)
library(ggplot2)

summarize_posterior = function(X) {
  X = mcmc(X)
  Xdata = data.frame(ID = 1:ncol(X),Mean = colMeans(X),Median = apply(X,2,median))
  Xi = HPDinterval(X,prob = 0.95)
  Xdata$low_95 = Xi[,1]
  Xdata$high_95 = Xi[,2]
  Xi = HPDinterval(X,prob = 0.8)
  Xdata$low_80 = Xi[,1]
  Xdata$high_80 = Xi[,2]
  return(Xdata)
}
posterior_plot = function(X,xlab='',ylab='',colorSig = T,ylim = NULL,colorGroup = NULL) {
  Xdata = summarize_posterior(X)
  if(!is.null(colorGroup)) {
    Xdata$color = colorGroup
    if(colorSig) {
      Xdata$color[sign(Xdata$low_95) == sign(Xdata$high_95)] = 'Significant'
      Xdata$color = factor(Xdata$color,levels = c('Significant',unique(colorGroup)))
    }
  } else{
    Xdata$color = 'NS'
    if(colorSig) {
      Xdata$color[sign(Xdata$low_95) == sign(Xdata$high_95)] = 'Sig'
    }
  }
  if(is.null(ylim)) ylim = range(Xdata[,4:5])

  p = ggplot(Xdata,aes(x=ID)) + geom_hline(yintercept = 0) +
    xlab(xlab) + ylab(ylab) + ylim(ylim)+
    geom_segment(aes(xend = ID,y = low_95,yend=high_95,color=color),size=.5) +
    geom_segment(aes(xend = ID,y = low_80,yend=high_80,color=color),size = .9) +
    geom_point(aes(y=Median,color = color)) +
    theme(legend.position = 'none')
  p
}

n = 200
p = 1000
b = 20

X = cbind(1,matrix(sample(c(0,1),n*p,replace=T),n,p))

K = matrix(rnorm(p*3),p); K = tcrossprod(K) + diag(1,p); cK = chol(K); X = matrix(rnorm(n*p),n,p) %*% (cK)
X = cbind(1,X)
n = nrow(X)
p = ncol(X)-1
X = sweep(X,2,colMeans(X),'-')
b_act = c(0,rep(0,p))
b_act[sample(1:p,b)] = 1*rnorm(b)
# b_act[1:b] = rnorm(b)

y = X %*% b_act + sqrt(.1)*rnorm(n)
Xm = as(X,'dgCMatrix')


l1 = 1
l2 = 10

b = b_act+1e-10
prec_y = 10

Ut = as(diag(1,n),'dgCMatrix')
s = rep(1,n)
h2 = 0
prior_mean = rep(0,p+1)


bs = matrix(0,1000,p+1)
prec_bs = matrix(0,1000,p+1)
prec_ys = c()
l1s = c()
l2s = c()
mu_hat = rep(0,length(y))

t = 1

log_post = function(log_l2,tau,l1,b,prec_y,t,shape,rate){
  l2 = exp(log_l2)
  (shape-1 + length(b)/2)*log_l2 - 1/2*sum(tau/(tau-1)*b^2)*l2*t*prec_y - length(b)/2*log_l2 - prec_y*t*l1^2*sum(tau)/(8*l2) - l2*rate
}

for(i in (-100):1000){
  if(i %% 100 == 0) print(i)
  # inverse Gaussian Bayesian Elastic Net from Li and Lin 2010
  # x  = pmax(1e-10,SuppDists::rinvGauss(n = p, nu = sqrt(l1)/(2*l2*abs(b)), lambda =  prec_y*t*l1/(4*l2)))  # wrong I think
  x  = pmax(1e-10,SuppDists::rinvGauss(n = p, nu = l1/(2*l2*abs(b)), lambda =  prec_y*t*l1^2/(4*l2)))
  tau = 1/x+1
  prec_b = c(1e-1,tau/(tau-1)*l2*t*prec_y)

  randn_theta = rnorm(p+1)
  randn_e = rnorm(n)
  b = sample_coefs_parallel_sparse_c_Eigen(Ut,y,Xm,h2, prec_y,s, prior_mean,prec_b,randn_theta,randn_e,1)

  l1 = sqrt(rgamma(1,shape = 3+p/2, rate = 1+prec_y*t*sum(tau)/(8*l2)))

  for(j in 1:20){
    old_state = log(l2)
    trial_state = old_state + rnorm(1,0,.5)
    logps = log_post(c(old_state,trial_state),tau,l1,b[-1],prec_y,t,shape=3,rate=1)
    if(log(runif(1)) < diff(logps)) l2 = exp(trial_state)
  }

  scores = tot_prec_scores_c(y - X %*% b,h2,s)
  alpha = n/2 + p
  beta = 1/2*(scores + l2*sum(tau/(tau-1)*b[-1]^2*t) + l1^2/(4*l2)*sum(tau))
  var_y = rinvgamma(1,shape = alpha,scale = beta)
  if(log(runif(1)) < p*log(gamma(0.5)) - p*log(pracma::gammainc(x=l1^2/(8*var_y*l2),a=0.5)['uppinc'])) {
    prec_y = 1/var_y
  }
  # prec_y = rgamma(1,shape = 2 + n/2,rate = 1.5 + 0.5*scores)

  if(i > 0){
    bs[i,] = c(b)
    prec_bs[i,] = c(prec_b)
    prec_ys = c(prec_ys,prec_y)
    l1s = c(l1s,l1)
    l2s = c(l2s,l2)
    mu_hat = mu_hat + X %*% b / 1000
  }

}

trace_plot(bs[,-1])
plot(1/prec_ys)
trace_plot(cbind(l1s,l2s))
posterior_plot(bs[,-1])
plot(colMeans(bs[,-1]),b_act[-1]);abline(0,1)
plot(X %*% b_act,mu_hat);abline(0,1)


### TPB prior:
library(GIGrvg)
bs = matrix(0,1000,p+1)
prec_bs = matrix(0,1000,p+1)
prec_ys = c()
l1s = c()
l2s = c()
mu_hat = rep(0,length(y))


t = 1
A=1/2
B=1/2
omega = 1/10
N=1e4
lambda = rgamma(N,shape=B,rate=omega)
tau = rgamma(N,shape=A,rate=lambda)
b = rnorm(N,0,sqrt(tau))
hist(b,breaks=c(-Inf,seq(-1,1,length=100),Inf),xlim=c(-1,1))



# omega = rgamma(1,shape=1/2,rate=1)
psi = rgamma(1,shape=B,rate=omega)
lambda = rgamma(p,shape=B,rate=psi)
tau = rgamma(p,shape=A,rate=lambda)
s2 = 1
b = c(0,rnorm(p,0,sqrt(s2*tau)))
c0=0;d0=0


prec_ys=c();mu_hat=0
for(i in (-100):1000){
  if(i %% 100 == 0) print(i)
  ### TPB prior:

  prec_y = 1/s2
  prec_b = c(1e-1,1/(tau))  # *s2
  randn_theta = rnorm(p+1)
  randn_e = rnorm(n)
  b = sample_coefs_parallel_sparse_c_Eigen(Ut,y,Xm,h2, prec_y,s, prior_mean,prec_b,randn_theta,randn_e,1)

  scores = tot_prec_scores_c(y - X %*% b,h2,s)  + d0 #+ sum(b[-1]^2/tau)
  s2 = 1/rgamma(1,shape = (n+p+c0)/2,rate = scores/2)
  
  b2 = b[-1]^2
  tau = sapply(1:p,function(j) rgig(n=1,lambda = A-1/2, chi = b2[j], psi = 2*lambda[j])) #/s2
  tau[tau==Inf]=1e-10

  lambda = rgamma(p,shape = A + B, rate = tau + psi)

  psi = rgamma(1,shape=p*B + 1/2, rate = sum(lambda) + omega)

  if(i > 0){
    bs[i,] = c(b)
    prec_bs[i,] = c(prec_b)
    prec_ys = c(prec_ys,prec_y)
    mu_hat = mu_hat + X %*% b / 1000
  }
}

trace_plot(bs[,-1])
plot(1/prec_ys)
posterior_plot(bs[,-1])
p1=posterior_plot(bs[,-1]);p1 + geom_point(data = data.frame(x=1:p,y=b_act[-1]),aes(x=x,y=y)) + ylim(c(-3,2))#geom_vline(xintercept = which(b_act[-1] != 0))
plot(colMeans(bs[,-1]),b_act[-1]);abline(0,1)
plot(X %*% b_act,mu_hat);abline(0,1)


prec_y = 10

bs = matrix(0,1000,p+1)
prec_bs = matrix(0,1000,p+1)
prec_ys = c()
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
  b = sample_coefs_parallel_sparse_c_Eigen(Ut,y,Xm,h2, prec_y,s, prior_mean,prec_b,randn_theta,randn_e,1)

  l1 = sqrt(rgamma(1,shape = 3+p/2, rate = 1+t*sum(tau)/(8*l2)))

  old_state = log(l2)
  trial_state = old_state + rnorm(1)
  logps = log_post(c(old_state,trial_state),tau,l1,b[-1],shape=3,rate=1)
  if(log(runif(1)) < diff(logps)) l2 = exp(trial_state)

  scores = tot_prec_scores_c(y - X %*% b,h2,s)
  prec_y = rgamma(1,shape = 2 + n/2,rate = 1.5 + 0.5*scores)

  bs[i,] = c(b)
  prec_bs[i,] = c(prec_b)
  prec_ys = c(prec_ys,prec_y)
  l1s = c(l1s,l1)
  l2s = c(l2s,l2)

}

trace_plot(bs[,-1])
trace_plot(cbind(l1s,l2s))
plot(1/prec_ys)
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
  b = sample_coefs_parallel_sparse_c_Eigen(Ut,y,Xm,h2, prec_y,s, prior_mean,prec_b,randn_theta,randn_e,1)

  scores = tot_prec_scores_c(y - X %*% b,h2,s)
  prec_y = rgamma(1,shape = 2 + n/2,rate = 1.5 + 0.5*scores)

  bs[i,] = c(b)
  # bs = rbind(bs,t(b))
  # prec_bs = rbind(prec_bs,prec_b)
  prec_ys = c(prec_ys,prec_y)
}

trace_plot(bs[,-1])
p1=posterior_plot(bs[,-1]);p1 + geom_point(data = data.frame(x=1:p,y=b_act[-1]),aes(x=x,y=y)) + ylim(-1,1)#geom_vline(xintercept = which(b_act[-1] != 0))
plot(colMeans(bs[,-1]),b_act[-1]);abline(0,1)


## Bayesian Lasso from BLR
lambda = .13
prec_b = c(1e-1,rep(1,p))
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
  b = sample_coefs_parallel_sparse_c_Eigen(Ut,y,Xm,h2, prec_y,s, prior_mean,prec_b*t,randn_theta,randn_e,1)

  scores = tot_prec_scores_c(y - X %*% b,h2,s)
  prec_y = rgamma(1,shape = 1 + n/2 + p/2,rate = 0 + 0.5*(scores + sum(b[-1]^2/tau_j)))

  bs[i,] = c(b)
  lambdas = c(lambdas,lambda)
  # bs = rbind(bs,t(b))
  # prec_bs = rbind(prec_bs,prec_b)
  prec_ys = c(prec_ys,prec_y)
}
plot(lambdas/prec_ys)

trace_plot(bs[,-1])
plot(lambdas)
p1=posterior_plot(bs[,-1]);p1 + geom_point(data = data.frame(x=1:p,y=b_act[-1]),aes(x=x,y=y)) + ylim(c(-3,2))#geom_vline(xintercept = which(b_act[-1] != 0))
plot(colMeans(bs[,-1]),b_act[-1]);abline(0,1)



library(glmnet)
res = cv.glmnet(X,y,alpha=1,keep=T)
res$lambda[order(res$cvm)[1]]
alphas = seq(0.1,1,length=20)
res_a = sapply(alphas,function(alpha) min(cv.glmnet(X,y,alpha=alpha,foldid = res$foldid)$cvm ))
plot(alphas,res_a)
plot(b_act,res$glmnet.fit$beta[,order(res$cvm)[1]]);abline(0,1)
res = glmnet(X,y,alpha=l1)
plot(res$beta,colMeans(bs));abline(0,1)

library(EBEN)
res2 = EBelasticNet.Gaussian(X[,-1],y,lambda = .24,alpha = .5)
plot(b_act[-1][res2$weight[,1]],res2$weight[,3]);abline(0,1)

library(BLR)
res3 = BLR(y,XL = X,prior=list(varE=list(df=3,S=0.25),
                               varU=list(df=3,S=0.63),
                               lambda=list(shape=0.52,rate=1e-4,
                                           type='random',value=30)),nIter=5500,burnIn=500,thin=1)


grav2 = grav
grav2$pheno = y
out = scanone(grav2,pheno.col = 1,method='hk')
plot(out$lod)

library(rstan)
rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

stan_model = stan('BEN.stan',chains=0)
stan_data = list(
  n = n,
  p = p,
  X = X[,-1],
  y=y[,1]
)

Nchains <- 1
Nwarm <- 200
Niter <- Nwarm+100

stan_fit = sampling(  object = get_stanmodel(stan_model),
                      data =stan_data,
                      iter=Niter,
                      warmup = Nwarm,
                      # init='0',
                      init = list(list(beta = b[-1],lambda1=l1,lambda2=l2,
                                       sigma = sqrt(1/prec_y),tau = tau)),
                      chains = Nchains, verbose=F,
                      refresh = 10,
                      pars=c('beta','lambda1','lambda2','sigma','tau'), #
                      open_progress = FALSE,
                      show_messages = FALSE
)
plot(stan_fit)
print(stan_fit)
la = extract(stan_fit)
plot(colMeans(la$beta),colMeans(bs[,-1]))
trace_plot(cbind(la$lambda1,la$lambda2))
