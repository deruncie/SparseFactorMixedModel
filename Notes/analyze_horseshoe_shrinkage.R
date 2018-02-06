
get_k = function(n,sigma2,tau2,lambda2,s2=1){
  1/(1+n/sigma2 * tau2 * s2 * lambda2)
}

# note:
# n = # individuals
# sigma2 = 1/tot_Eta_prec
# tau2 = omega2/tauh
# lambda2 = phi2,
# s2 = 1


res = get_posterior_FUN(BSFG_state,{
  po = colSums(abs(Lambda) > 0.01)
  po/(nrow(Lambda)-po)/sqrt(n)})

tau = get_posterior_FUN(BSFG_state,sqrt(Lambda_omega2[1]/tauh[1,]))
plot(colMeans(tau),colMeans(res));abline(0,1)
meff = get_posterior_FUN(BSFG_state,{
  a = sqrt(Lambda_omega2[1]/tauh[1,])*sqrt(nrow(data))
  a/(1+a)*nrow(Lambda)})
n_eff= get_posterior_FUN(BSFG_state,colSums(abs(Lambda) > 0.1))
plot(colMeans(n_eff),colMeans(meff));abline(0,1)

p = ncol(Y)
n = nrow(data)
k = BSFG_state$current_state$k
meff2 = get_posterior_FUN(BSFG_state,
                       {
                       ks = get_k(nrow(data),
                                  1/t(tot_Eta_prec[rep(1,k),]),
                                  (Lambda_omega2[1]/tauh[rep(1,p),]),
                                  Lambda_phi2)
                       colSums(1-ks)
                       })
Lambda = get_posterior_mean(BSFG_state,Lambda)
phi2 = get_posterior_mean(BSFG_state,Lambda_phi2)
tau2 = get_posterior_mean(BSFG_state,(Lambda_omega2[1]/tauh[rep(1,p),]))
meff3 = colSums(1-get_k(nrow(data),1,tau2,phi2))

a=sqrt(tau2[1,])*sqrt(n)/sqrt(4);plot(a/(1+a)*p,ylim=c(0,p))
a=sqrt(tau_0^2/(cumprod(c(1,rep(3,20)))))*sqrt(n)/sqrt(4);plot(a/(1+a)*p,ylim=c(0,p0))

# all above tests seem to be reasonable, with counting genes as included at ~abs(lambda) = 0.01


# This is concerning: we need very large delta values to get small meff
# if we start non-sparse
a=sqrt(.3/(cumprod(c(1,rep(3,20)))))*sqrt(n);plot(a/(1+a)*p,ylim=c(0,p))
a=sqrt(.3/(cumprod(c(1,rep(30,20)))))*sqrt(n);plot(a/(1+a)*p,ylim=c(0,p))
a=sqrt(1e-1/(cumprod(c(1,rep(3,20)))))*sqrt(n);plot(a/(1+a)*p,ylim=c(0,p))


# calibrating delta1_mean prior
p = ncol(Y)
n = nrow(Y)
p0 = .2*p
tau_0 = p0/(p-p0)*sqrt(1)/sqrt(n)
delta1_mean = 1/tau_0^2
nR = 1e4
k = 20
meff = t(sapply(1:nR,function(i) {
  tauh = rgamma(1,shape=5,rate=5/delta1_mean) * cumprod(c(delta1_mean,rgamma(k-1,shape=1,rate=1/delta1_mean)))
  tau = sqrt(1/tauh)
  l = matrix(rcauchy(p*k),nc=k)
  ki = get_k(n,sigma2 = 4,tau2 = matrix(tau^2,p,k,byrow=T),lambda2 = l^2)
  meff = colSums(1-ki)
  meff
}))
boxplot(meff)
tau = sqrt(1/rgamma(nR,shape=1,rate=1/delta1_mean))# * 1/rgamma(nR,shape=1,rate=1/3))
# tau = abs(rnorm(nR,0,tau_0))
tau = abs(rcauchy(nR,0,.1)) * 1/sqrt(apply(matrix(rgamma(nR*10,shape=1,rate=1/3),nr=nR),1,prod))
# tau = abs(rcauchy(nR,0,scale=tau_0^2))
# tau = rep(tau_0,nR)
l = matrix(rcauchy(p*nR),nc=nR)
k = get_k(n,sigma2 = 1,tau2 = matrix(tau^2,p,nR,byrow=T),lambda2 = l^2)
meff = colSums(1-k)
hist(meff,breaks=c(seq(0,100,length=100),Inf),xlim=c(0,100));abline(v=p0)

# plot delta against p0 - it's non-linear, with changes in delta causing biggest
# changes in p0 at 0.5*p
# This suggests the prior is for some number of non-sparse factors,
# and then som number of very sparse factors
p0 = seq(p-1,1)
p0 = seq(70,30)
tau_0 = p0/(ncol(Y)-p0)*sqrt(4)/sqrt(nrow(Y))
delta1_mean = 1/tau_0^2
plot(p0[-1],exp(diff(log(delta1_mean))))

# the effect of delta on tau_0 should be independent of p (ie number of genes)
