library(MCMCglmm)
nG = 20
nR = 8
data = data.frame(G = gl(nG,nR))
n = nrow(data)
Z = model.matrix(~0+G,data)
A = cbind(sample(c(0,1),nG,replace=T,prob = c(.7,.3)),sample(c(0,1),nG,replace=T,prob = c(.7,.3)),sample(c(0,1),nG,replace=T,prob = c(.9,.1)));
A = A * matrix(rnorm(length(A)),nrow(A))
A = A %*% t(A) + diag(1,nG)
# A = diag(1,nG)
Ainv = diag(1,nG)
svd_ZtZ = svd(t(Z) %*% Z)
Q = svd_ZtZ$u
d = svd_ZtZ$d

s2G = 2
s2R = 2
pa = 1/s2G
pe = 1/s2R


g = rnorm(nG,0,sqrt(s2G))

data$y = g[data$G] + rnorm(nrow(data),0,sqrt(s2R))

library(lme4)
lme1 = lmer(y~(1|G),data)
summary(lme1)

library(MCMCglmm)

priors = list(R = list(V = 1, nu = 4), G = list(G1 = list(V = 1,nu=4)))

m1 = MCMCglmm::MCMCglmm(y~1,random = ~G,data=data,prior = priors)


resid_prec_shape = priors$R$nu/2
resid_prec_rate = resid_prec_shape*priors$R$V
a_prec_shape = priors$G$G1$nu/2
a_prec_rate = a_prec_shape*priors$G$G1$V

pa = rgamma(1,shape = a_prec_shape,rate = a_prec_rate)
pe = rgamma(1,shape = resid_prec_shape,rate = resid_prec_rate)

a = rnorm(nG,0,1/sqrt(pa))

burn = 3
Iter = 1000
thin = 1

Zty = t(Z) %*% data$y

r = c()
as = c()

for(i in 1:(Iter+burn)) {
	if(i %% 100 == 0) print(i)
	dsqrt = 1/sqrt(d*pe + pa)
	Q_sqrt = Q * dsqrt
	a_star = Q_sqrt %*% rnorm(nG)
	a = a_star + Q_sqrt %*% t(Q_sqrt) %*% Zty * pe
	as = cbind(as,as.matrix(a))

	pa = rgamma(1,shape = a_prec_shape + nG/2, rate = a_prec_rate + 1/2*sum(a^2))
	y_resid = data$y - Z %*% a
	pe = rgamma(1,shape = resid_prec_shape + n/2, rate = resid_prec_rate + 1/2*sum(y_resid^2))

	if( (i-burn) %% thin == 0 ){
		r = rbind(r,c(1/pa,1/pe))
	}

}
plot(rowMeans(as),a);abline(0,1)
plot(r[,1])
qqplot(m1$VCV[,1],r[,1]);abline(0,1)
qqplot(m1$VCV[,2],r[,2]);abline(0,1)


## Testing Hadfield's strategy
# A = matrix(rnorm(nG*10),nG);A = A %*% t(A) + diag(1,nG)
# A = cbind(sample(c(0,1),nG,replace=T,prob = c(.9,.1)),sample(c(0,1),nG,replace=T,prob = c(.9,.1)),sample(c(0,1),nG,replace=T,prob = c(.9,.1)));A = A %*% t(A) + diag(1,nG)
# A = cbind(c(rep(0,2),rep(1,5),rep(0,nG-7)),c(rep(1,3),rep(0,nG-3)),c(rep(0,4),rep(1,5),rep(0,nG-9)));A = A %*% t(A) + diag(1,nG)
As = as(A,'sparseMatrix')
Ai = solve(As)
Aic = chol(Ai)
Ac = chol(As)

Zs = Matrix(Z,sparse=T)
W = Zs
Cbase = t(W) %*% Diagonal(n,1) %*% W + Ai
cholCbase = Cholesky(Cbase,LDL=F)
image(cholCbase)
image(update(cholCbase,Cbase))

	P = expand(cholCbase)$P


as = c()
r = c()

for(i in 1:(Iter+burn)) {
	if(i %% 100 == 0) print(i)

	thetas1 = c()
	thetas2 = c()
	C = pe*t(W) %*% W + Ai
	Ci = solve(C)
	cholC = chol(C)
	cholCi = chol(Ci)
	a = as.matrix(solve((cholC),matrix(rnorm(nG*1000),nG)))
	a = as.matrix(t(cholCi) %*% matrix(rnorm(nG*1000),nG))
	plot(cov(t(a)),Ci);abline(0,1)
for(i in 1:1000){

microbenchmark(solve(C,WtRiy), Ci %*% WtRiy)

microbenchmark({
	theta_star = Ac %*% rnorm(nG)
	# theta_star = solve(Aic,rnorm(nG))
	W_theta_star = W %*% theta_star
	# e_star = W_theta_star + rnorm(n,0,1/sqrt(pe))
	e_star = rnorm(n,0,1/sqrt(pe))



	y_resid = data$y - W_theta_star - e_star
	WtRiy = pe*t(W) %*% y_resid

	# plot(solve(cholC,WtRiy,'A'),solve(C) %*% WtRiy);abline(0,1)
	# plot(solve(cholCbase,WtRiy,'L'),solve(Cbase) %*% WtRiy);abline(0,1)

	# adsf = as.matrix(solve(cholC,matrix(rnorm(nG*1000),nG),'L'))
	# plot(solve(C),cov(t(adsf)));abline(0,1)
	# asdf = as.matrix(chol(solve(C)) %*% matrix(rnorm(nG*1000),nG))
	# plot(solve(C),cov(t(asdf)));abline(0,1)


	theta_tilda1 = solve(C,WtRiy)
	# theta_tilda1 = Ci %*% WtRiy
	theta = theta_tilda1 + theta_star
},
{#Ci = solve(C)
	# theta = pe * Ci %*% t(W) %*% data$y + solve(cholC,rnorm(nG))
	theta = pe * Ci %*% t(W) %*% data$y + t(cholCi) %*% rnorm(nG)
})
	thetas1 = cbind(thetas1,as.matrix(theta_tilda1 + theta_star))


	# thetas2 = cbind(thetas2,as.matrix(pe * Ci %*% t(W) %*% data$y + solve(cholC,rnorm(nG))))
	thetas2 = cbind(thetas2,as.matrix(pe * Ci %*% t(W) %*% data$y + t(cholCi) %*% rnorm(nG)))
}
plot(rowMeans(thetas1),rowMeans(thetas2));abline(0,1)
plot(cov(t(thetas1)),cov(t(thetas2)));abline(0,1)

	PCPt = solve(cholCbase,t(solve(cholCbase,C,'P')),'P')
	cCp = chol(PCPt)
	cC = update(cholCbase,solve(cholCbase,cCp,'Pt'))
	theta_tilda2 = solve(cC,WtRiy,'A')

	# plot(theta_tilda1,theta_tilda2);abline(0,1)

	microbenchmark(

	theta_tilda1 = solve(C,WtRiy),
	{
	PCPt = solve(cholCbase,t(solve(cholCbase,C,'P')),'P')
	cCp = chol(PCPt)
	cC = update(cholCbase,solve(cholCbase,cCp,'Pt'))
	theta_tilda2 = solve(cC,WtRiy,'A')

	}

		)

	# Ci = t(P) %*% t(solve(cCp)) %*% solve(cCp) %*% P
	# Ci = crossprod(solve(cCp) %*% P)
	# Ci = tcrossprod(solve(cholCbase,solve(cCp),'Pt'))
	# plot(solve(C,WtRiy),Ci %*% WtRiy);abline(0,1)
	# cC = update(cholCbase,solve(cholCbase,cCp,'Pt'))
	# plot(solve(C,WtRiy),solve(cC,WtRiy,'A'));abline(0,1)

	# cholC = update(cholCbase,C)

	theta_tilda = theta_tilda1

	a = theta_tilda + theta_star
	as = cbind(as,as.matrix(a))

	pa = rgamma(1,shape = a_prec_shape + nG/2, rate = a_prec_rate + 1/2*sum(a^2))
	y_resid = data$y - Zs %*% a
	pe = rgamma(1,shape = resid_prec_shape + n/2, rate = resid_prec_rate + 1/2*sum(y_resid^2))

	if( (i-burn) %% thin == 0 ){
		r = rbind(r,c(1/pa,1/pe))
	}

}
plot(rowMeans(as),a);abline(0,1)
plot(r[,1])
qqplot(m1$VCV[,1],r[,1]);abline(0,1)
qqplot(m1$VCV[,2],r[,2]);abline(0,1)



as = c()
r = c()
p = 100
for(i in 1:(Iter+burn)) {
	if(i %% 100 == 0) print(i)
		Y = matrix(rep(data$y,p),nc=p)
		Sigma_Rs = list()
		Sigma_Rs[[1]] = list(Sigma = Diagonal(n,1))
		sigma_indexes = rep(1,p)
		h2 = pe/(pe+pa)
		tot_Y_prec = rep(pe*(1-h2),p)
		h2s = matrix(h2,nr=1,nc=p)
		Ais = list()
		Ais[[1]] = Ai
		prior_mean = matrix(0,nG,p)
		as = sample_MME_matrixPrior(Y,W,Sigma_Rs,sigma_indexes,tot_Y_prec,prior_mean,Ais, h2s,Ac,1)


pas = rgamma(1e5,shape = a_prec_shape,rate = 1/a_prec_rate)
pes = rgamma(1e5,shape = resid_prec_shape,rate = 1/resid_prec_rate)
# h2 = pes/(pes+pas)
# hist(h2,prob=T)
# z = seq(0,1,length=100)
# lines(z,dbeta(seq(0,1,length=h2_divisions+1),2,2)[-(h2_divisions+1)])


h2_divisions = 10
h2_priors = dbeta(seq(0,1,length=h2_divisions+1),2,2)[-(h2_divisions+1)]
# h2_priors = h2_priors/sum(h2_priors)
result = svd(Z %*% solve(Ainv) %*% t(Z))
invert_aI_bZAZ = list(
    U = result$u,
    s = result$d
)

# h2s = sample(seq(0,1,length=h2_divisions+1)[-(h2_divisions+1)],size = 1e5,replace=T,prob = h2_priors)
# ps = 
# pas2 = pes*(1-h2s)/(h2s)
# qqplot(pas,pas2);abline(0,1)
# qqplot(h2,pes/(pes+pas2));abline(0,1)


pres_shape = 2.95 + 00
pres_rate = 1.84
pres_rate = var(data$y) * (pres_shape - 1)

# pres_shape = 200
# pres_rate = 201

# pres = rgamma(1e5,shape = pres_shape,rate = pres_rate)
# h2s = sample(seq(0,1,length=h2_divisions+1)[-(h2_divisions+1)],size = 1e5,replace=T,prob = h2_priors)
# pas2 = pres/h2s
# pes2 = pres/(1-h2s)
# hist(1/pres)
# hist(1/pas2)
# hist(h2s)
# plot(1/pres,h2s)
# qqplot(pres,1/(1/pas + 1/pes));abline(0,1)
# qqplot(1/pas,1/pas2);abline(0,1)
# plot(1-h2s,pres/pes2)
# plot(pres,h2s)
# hist(pes2/(pes2+pas2))
# hist(pes/(pes+pas))
# hist(1/pas,breaks=c(seq(0,5,length=100),Inf),xlim = c(0,1))
# hist(1/pas2,breaks=c(seq(0,5,length=100),Inf),xlim = c(0,1))

pres = rgamma(1,shape = pres_shape,rate = pres_rate)
ZZt = t(Z) %*% Z
result = GSVD_2_c(cholcov(ZZt),cholcov(Ainv))
invert_aZZt_Ainv = list(
	U = t(solve(result$X)),
		s1 = diag(result$C)^2,
		s2 = diag(result$S)^2
	)
r2 = c()
for(i in 1:(Iter+burn)) {
	if(i %% 100 == 0) print(i)

	h2 = sample_h2s_discrete_given_p_c(matrix(data$y,nr=n,nc=1),h2_divisions,h2_priors,pres,invert_aI_bZAZ)
	U = invert_aI_bZAZ$U
	s = invert_aI_bZAZ$s
	UtY_tilde = t(U) %*% data$y
	Sigma_sqrt = sqrt(h2*s + (1-h2))
	SiUtY_tilde = UtY_tilde / Sigma_sqrt
	pres = rgamma(1,shape = pres_shape + n/2, rate = pres_rate+1/2*sum(SiUtY_tilde^2))
	pa = pres / h2
	pe = pres / (1-h2)

	a = sample_F_a_ipx_c(matrix(data$y),Z,pa,pe,invert_aZZt_Ainv)

	if( (i-burn) %% thin == 0 ){
		r2 = rbind(r2,c(1/pa,1/pe,1/pres,h2))
	}

}
plot(r2[,1])
library(MCMCpack)

effectiveSize(mcmc(r))
effectiveSize(mcmc(r2))
qqplot(m1$VCV[,1],r2[,1]);abline(0,1)
qqplot(m1$VCV[,2],r2[,2]);abline(0,1)

boxplot(as.data.frame(m1$VCV))
boxplot(r)
boxplot(r2)