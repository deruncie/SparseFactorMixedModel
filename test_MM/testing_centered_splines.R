library(rrBLUP)
library(splines)
library(MASS)

d = data.frame(time = seq(1,24,length=30))
d$y = sin(d$time/24*6) + d$time*1
d$y = d$y-mean(d$y)+1
e = rnorm(nrow(d),0,.5);e=e-mean(e)
d$y2 = d$y + e
with(d,plot(time,y2))



Z1 = model.matrix(~0+bs(d$time,df=10,intercept = F))
Image(Z1,F)

D1 = diag(1,ncol(Z1))
D1[lower.tri(D1)] = 1

Z2 = Z1 %*% D1

# Z3 = Z1 %*% contr.sum(ncol(Z1))
# D3 = diag(1,ncol(Z3))
# D3[lower.tri(D3)] = 1
# diag(D3) = 0

# Z4 = Z3 %*% D3
# Z3 = model.matrix(~0+bs_diff(d$time,df=10,center=T,intercept = T)) %*% contr.sdif(9)
# Z3 = model.matrix(~0+bs_diff(d$time,df=10,center=T,intercept = T)) %*% cbind(1,contr.sdif(9)) #%*% cbind(1,contr.sdif(9))
Z3 = model.matrix(~0+bs_diff(d$time,df=10,center=T,intercept = T,differences=2))

Z4 = model.matrix(~0+bs_diff(d$time,df=10,center=T,intercept = T,periodic = F))

res1 = mixed.solve(d$y2,Z=Z1)
res2 = mixed.solve(d$y2,Z=Z2)
res3 = mixed.solve(d$y2,Z=Z3)
res4 = mixed.solve(d$y2,Z=Z4)

with(d,plot(time,y2))
# lines(d$time,Z1 %*% res1$u + res1$beta[1],col=1)
# lines(d$time,Z2 %*% res2$u + res2$beta[1],col=2)
lines(d$time,Z3 %*% res3$u + res3$beta[1],col=3)
lines(d$time,Z4 %*% res4$u + res4$beta[1],col=4)


# mean(Z1 %*% res1$u)
# mean(Z2 %*% res2$u)
# mean(Z3 %*% res3$u)
mean(Z4 %*% res4$u)
res1$beta
res2$beta
res3$beta
res4$beta

obj_LL2 = function(sigma2,log_det,e2){
  -1/2*(length(e2)*log(2*pi*sigma2)+ log_det) - sum(e2/sigma2)/2
}

logLL = function(y,Z,h2,mu = mean(y)){
  y = y-mu
  Sigma = h2*tcrossprod(Z) + (1-h2)*diag(1,nrow(Z))
  chol_Sigma = chol(Sigma)
  e2 = solve(t(chol_Sigma),y)^2
  log_det = 2*sum(log(diag(chol_Sigma)))
  res = optimize(obj_LL2,interval = c(1e-10,10),maximum=T,log_det = log_det,e2=e2)
  return(res)
}
h2s = seq(0,.99,length=20)
mu = mean(d$y2)
# mu = res1$beta[1]
Z = Z5
res = sapply(h2s,function(h2) logLL(d$y2,Z,h2,mu)$obj)
# plot(h2s,res);abline(v=which(res==max(res))/length(res))

h2 = h2s[res==max(res)]
s2 = logLL(d$y2,Z5,h2,mu)$max
s2e = (1-h2)*s2
s2a = h2*s2

# with(d,plot(time,y2))

lines(d$time,Z %*% solve(crossprod(Z)/s2e + diag(1/s2a,ncol(Z)),crossprod(Z,(d$y2-mu))/s2e) + mu,col=1)

library(lme4qtl)
d$ID = 1:nrow(d)
rownames(Z5) = d$ID
lme1 = relmatLmer(y2~(1|ID),d,relmat = list(ID = tcrossprod(Z5)))
lme1
lines(d$time,Z5 %*% ginv(Z5) %*% ranef(lme1)[[1]][,1] + fixef(lme1)[1],col=5)
