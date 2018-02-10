library(rrBLUP)

d = data.frame(time = seq(1,24,length=30))
d$y = sin(d$time/24*6)
d$y = d$y-mean(d$y)+1
e = rnorm(nrow(d),0,.5);e=e-mean(e)
d$y2 = d$y + e
with(d,plot(time,y2))



Z1 = model.matrix(~0+pbs(d$time,df=10,intercept = F,Boundary.knots = c(0,24)))
Image(Z1,F)

D1 = diag(1,ncol(Z1))
D1[lower.tri(D1)] = 1

Z2 = Z1 %*% D1

Z3 = Z1 %*% contr.sum(ncol(Z1))
D3 = diag(1,ncol(Z3))
D3[lower.tri(D3)] = 1
diag(D3) = 0

Z4 = Z3 %*% D3

res1 = mixed.solve(d$y2,Z=Z1)
res2 = mixed.solve(d$y2,Z=Z2)
res3 = mixed.solve(d$y2,Z=Z3)
res4 = mixed.solve(d$y2,Z=Z4)

with(d,plot(time,y2))
lines(d$time,Z1 %*% res1$u + res1$beta[1],col=1)
lines(d$time,Z2 %*% res2$u + res2$beta[1],col=2)
lines(d$time,Z3 %*% res3$u + res3$beta[1],col=3)
# lines(d$time,Z4 %*% res4$u + res4$beta[1],col=4)

b_spline = function(x, df = NULL, knots = NULL, degree = 3, intercept = FALSE,
                    Boundary.knots = range(x),
                    differences = TRUE,
                    center = FALSE
) {
  # following code from https://github.com/SurajGupta/r-source/blob/master/src/library/splines/R/splines.R
  bs_X = bs(x,df,knots,degree,intercept,Boundary.knots)
  X = bs_X
  if(center){
    X = X %*% contr.sum(ncol(X))
  }
  if(differences){
    D = diag(1,ncol(X))
    D[lower.tri(D)] = 1
    diag(D) = 1
    X = X %*% D
  }
  # X
  bs_X_attributes = attributes(bs_X)
  bs_X_attributes = bs_X_attributes[names(bs_X_attributes) %in% c('dim','dimnames') == F]
  attributes(X) = c(attributes(X),bs_X_attributes)
  attr(X,'differences') = differences
  attr(X,'center') = center
  class(X) = c('b_spline',class(X))
  X
}
makepredictcall.b_spline <- function(var, call)
{
  if(as.character(call)[1L] != "b_spline") return(call)
  at <- attributes(var)[c("degree", "knots", "Boundary.knots", "intercept","differences","center")]
  xxx <- call[1L:2]
  xxx[names(at)] <- at
  xxx
}

model=y2~b_spline(time,df=10)
# model=y2~bs(time,df=10)
terms = delete.response(terms(model.frame(model,d)))
m1 = model.matrix(model,d)
m2 = model.matrix(terms,data.frame(time=d$time[1:4]))
m2-m1[1:4,]


Z5 = model.matrix(~0+b_spline(d$time,df=ncol(Z1),intercept = T))
# Z5 = sweep(Z5,1,apply(Z5,1,sd),'/')
res5 = mixed.solve(d$y2,Z=Z5,bounds = c(1e-9,1e1))
lines(d$time,Z5 %*% res5$u + res5$beta[1],col=5)

mean(Z1 %*% res1$u)
mean(Z2 %*% res2$u)
mean(Z3 %*% res3$u)
mean(Z4 %*% res4$u)
mean(Z5 %*% res5$u)
res1$beta
res2$beta
res3$beta
res4$beta
res5$beta

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
