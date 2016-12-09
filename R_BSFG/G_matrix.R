CovToCor = function(X){
  # Normalize a covariance matrix into a correlation matrix
  corX=diag(1/sqrt(diag(X))) %*% X %*% diag(1/sqrt(diag(X)));
  return(corX)
}

run_variables  = BSFG_state$run_variables
Posterior      = BSFG_state$Posterior
E_a_prec       = BSFG_state$current_state$E_a_prec
resid_Y_prec   = BSFG_state$current_state$resid_Y_prec
traitnames     = BSFG_state$traitnames

sp_num = ncol(BSFG_state$Posterior$Lambda)   
p = run_variables$p
k = nrow(Posterior$Lambda)/p;

G_Lambdas = array(0,dim = dim(Posterior$Lambda))
Lambda_est = matrix(0,p,k)
G_est = E_est = matrix(0,p,p)
traces_G = matrix(,p*(p+1)/2,sp_num)
diag_G = matrix(,p,sp_num)
traces_G_cor = matrix(,p*(p+1)/2,sp_num)
traces_E = matrix(,p*(p+1)/2,sp_num)
diag_E = matrix(,p,sp_num)



for(j in 1:sp_num) {
  Lj = matrix(Posterior$Lambda[,j],p,k)
  h2j = Posterior$F_h2[,j]
  G_Lj = Lj %*%  diag(sqrt(h2j))
  G_Lambdas[,j] = c(G_Lj)
  Gj = G_Lj %*%  t(G_Lj) + diag(1/Posterior$E_a_prec[,j])
  G_est = G_est + Gj/sp_num
  traces_traits[,j] = Gj%*%beta_v
  library(gdata)
  traces_G[,j] = lowerTriangle(Gj,diag = TRUE)
  diag_G[,j] = diag(Gj)
  traces_G_cor[,j] = lowerTriangle(CovToCor(Gj),diag = TRUE)
  
  E_Lj = Lj  %*% diag(1-h2j) %*%  t(Lj) + diag(1/Posterior$resid_Y_prec[,j])
  E_est = E_est + E_Lj/sp_num;
  Lambda_est = Lambda_est + matrix(Posterior$Lambda[,j],p,k)/sp_num;
  traces_E[,j] = lowerTriangle(E_Lj,diag = TRUE)
  diag_E[,j] = diag(E_Lj)
}
G_Lambda = matrix(rowMeans(G_Lambdas),p,k)



load("C:/Users/Xin~/Desktop/setup_LineCode1.RData")
Y = setup$Y
VY = apply(Y,2,var,na.rm=T)

#selection gradient
beta = 0.3609
beta_v = c(0,beta,rep(0,16))
traces_traits = matrix(,p,sp_num)
for(j in 1:sp_num) {
  Lj = matrix(Posterior$Lambda[,j],p,k)
  h2j = Posterior$F_h2[,j]
  G_Lj = Lj %*%  diag(sqrt(h2j))
  Gj = G_Lj %*%  t(G_Lj) + diag(1/Posterior$E_a_prec[,j])
  Gj = sqrt(diag(VY))%*%Gj%*%t(sqrt(diag(VY)))
  traces_traits[,j] = Gj%*%beta_v
}
boxplot(t(traces_traits))
save(traces_traits,file = 'traces_traits.RData')

#heritabiity of posterior mean
G_est1 = sqrt(diag(VY))%*%G_est%*%t(sqrt(diag(VY)))
save(G_est1,file = "G_est_posterior.RData")

E_est1 = sqrt(diag(VY))%*%E_est%*%t(sqrt(diag(VY)))
heritability = diag(G_est1)/diag(G_est1+E_est1)
names(heritability) = traitnames
plot(heritability)

#heritability confidence interval 95%  
diag_E1 = sweep(diag_E,1,VY,'*')
diag_G1 = sweep(diag_G,1,VY,'*')
P = diag_G1+diag_E1
trace_h2 = matrix(mapply(function(x,y)x/y, diag_G1,P),nr=p)

CI_h2 = t(apply(trace_h2,1,function(x)quantile(x,c(0.025,0.975))))
q_h2 = cbind(CI_h2,heritability)
colnames(q_h2)[3] = "posterior"
save(q_h2, file = "q_h2.RData")

#G_matrix confidence interval 95%
CI_G = t(apply(traces_G,1,function(x)quantile(x,c(0.025,0.975))))
q_G = cbind(CI_G,lowerTriangle(G_est,diag = TRUE))
colnames(q_G)[3] = "posterior"

#2.5%
G_low = matrix(,nr=p,nc=p)
G_low[lower.tri(G_low,diag=TRUE)]=q_G[,1]
G_low[upper.tri(G_low)]=upperTriangle(t(G_low))
G_low1 = sqrt(diag(VY))%*%G_low%*%t(sqrt(diag(VY)))
save(G_low1,file = "G_est_low.RData")
#97.5%
G_high = matrix(,nr=p,nc=p)
G_high[lower.tri(G_high,diag=TRUE)]=q_G[,2]
G_high[upper.tri(G_high)]=upperTriangle(t(G_high))
G_high1 = sqrt(diag(VY))%*%G_high%*%t(sqrt(diag(VY)))
save(G_high1,file = "G_est_high.RData")

library(lattice)
levelplot(CovToCor(G_est1))

load("C:/Users/Xin~/Desktop/Gmatrix_MCMCglmm.r")
library(gdata)
plot(upperTriangle(G_est1),upperTriangle(G),xlab = "BSFG",ylab = "MCMCglm")
lines(x = c(-2,100), y = c(-2,100))
