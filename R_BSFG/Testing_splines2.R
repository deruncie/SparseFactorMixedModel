library(mgcv)
library(BSFG)
library(splines)
t = seq(1,100)

X = model.matrix(~1+poly(t,4))
b = rnorm(ncol(X))

y1 = X %*% b
y = X %*% b + rnorm(length(t),0,.1)

plot(t,y)

Xb = model.matrix(~1+bs(t,df=30,intercept = T))
bh = coef(lm(y~0+Xb))
# lines(t,Xb %*% bh,col=2)

# bh

D = matrix(0,ncol = ncol(Xb)-1,nrow = ncol(Xb)-2)
# D[row(D) == col(D)] = 1
# D[row(D)-1 == col(D)] = -1
diag(D) = -1
D[row(D)+1==col(D)]=1

Pinv = crossprod(D)
svdP = svd(Pinv)
i = svdP$d> 1e-10
XX = svdP$u[,i] %*% diag(1/sqrt(svdP$d[i]))




# Image(XX)

# Image(Xb,aspect=1.5)
# Image(cbind(Xb,Xb %*% XX),aspect=1.5)
# Image(Xb %*% XX,aspect=1.5)
# Image(model.matrix(~0+bs(t,df=9,degree=3)),aspect=1.5)

# lines(t,bh2[1]+Xb %*% XX %*% bh2[-1],col=3)

# library(glmnet)
# # XX = XX*100
# res1 = glmnet(cbind(Xb[,-1] %*% XX,Xb[,-1]),y,alpha = 1,intercept=T,standardize=F,lambda.min.ratio = 0.01)
# #
#
# i = ncol(res1$beta)
# bh2 = c(res1$a0[i],XX %*% res1$beta[1:ncol(XX),i] + res1$beta[-c(1:ncol(XX)),i])
# # par(mfrow=c(1,2))
# # plot(res$beta[,i])
# # plot(XX %*% res$beta[1:ncol(XX),i],res$beta[-c(1:ncol(XX)),i],xlim = range(res$beta[,i]),ylim = range(res$beta[,i]));abline(0,1)
# lines(t,Xb %*% bh2,col=3)
#
# res3 = glmnet(Xb[,-1] %*% XX,y,alpha = 0,lambda.min.ratio = 0.00001,intercept=T,standardize=F)
# i = ncol(res3$beta)
# bh4 = c(res3$a0[i],XX %*% res3$beta[,i])
# lines(t,Xb%*% bh4,col=5)
#
# res2 = glmnet(Xb[,-1],y,alpha = 1,lambda.min.ratio = 0.001,intercept=T)
# i = ncol(res2$beta)
# bh3 = c(res2$a0[i],res2$beta[,i])
# lines(t,Xb %*% bh3,col=4)


res4 = gam(y~s(t,k = ncol(Xb),bs='cs'))
lines(t,predict(res4,newdata = list(t=t)))

# I learned that the interpretation of b-spline coefficients is more like
# the moving average than the moving slope!
# also, incorporating a moving average may be reasonable.
# Using the difference operator in the covariance works great!


library(rrBLUP)
res5 = mixed.solve(y,Z=Xb[,-1],X = Xb[,1])
bh5 = c(res5$beta,res5$u)
# lines(t,Xb %*% bh5,col=2)

# i=c(1:82)
res6 = mixed.solve(y,Z=Xb[,-1] %*% XX,X = Xb[,1])
bh6 = c(res6$beta,XX %*% res6$u)
lines(t,Xb %*% bh6,col=3)

p=ncol(Xb)-1
res7 = mixed.solve(y,Z=cbind(Xb[,-1], Xb[,-1] %*% XX),X = Xb[,1])
bh7 = c(res7$beta,res7$u[1:p] + c(XX %*% res7$u[-c(1:p)]))
# lines(t,Xb %*% bh7,col=3)

#
# X3 = matrix(0,ncol = ncol(Xb)-2,nrow = ncol(Xb)-1)
# diag(X3) = 1
# X3[row(X3)-1 == col(X3)] = -1

X3 = matrix(0,ncol(Xb)-1,ncol(Xb)-1)
diag(X3) = 1
X3[lower.tri(X3)] = 1


p=ncol(Xb)-1
res8 = mixed.solve(y,Z=cbind(0*Xb[,-1],Xb[,-1] %*% X3),X = Xb[,1])
bh8 = c(res8$beta,res8$u[1:p] + c(X3 %*% res8$u[-c(1:p)]))
# lines(t,Xb %*% bh8,col=4)

p=ncol(Xb)-1
res8 = mixed.solve(y,Z=Xb[,-1] %*% X3,X = Xb[,1])
bh8 = c(res8$beta,c(X3 %*% res8$u))
lines(t,Xb %*% bh8,col=4)

res3 = lm(y~0+Xb[,-1] %*% X3)
bh9 = X3 %*% coef(res3)
lines(t,Xb[,-1] %*% bh9,col=5)
#
# library(BGLR)
# res9= BGLR(y,ETA = list(list(X=Xb[,1],model='FIXED'),list(X=Xb[,-1],model='BRR'),list(X=Xb[,-1] %*% XX,model='BRR')),verbose=F)
# bh9 = c(res9$mu+res9$ETA[[1]]$b,res9$ETA[[2]]$b+XX %*% res9$ETA[[3]]$b)
# lines(t,Xb %*% bh9,col=5)

