library(MCMCglmm)
nG = 20
nR = 8
data = data.frame(G = gl(nG,nR))
n = nrow(data)
Z = model.matrix(~0+G,data)
Ainv = diag(1,nG)
svd_ZtZ = svd(t(Z) %*% Z)
Q = svd_ZtZ$u
d = svd_ZtZ$d

s2G = 5
s2R = .1

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

for(i in 1:(Iter+burn)) {
	if(i %% 100 == 0) print(i)
	dsqrt = 1/sqrt(d*pe + pa)
	Q_sqrt = Q * dsqrt
	a_star = Q_sqrt %*% rnorm(nG)
	a = a_star + Q_sqrt %*% t(Q_sqrt) %*% Zty * pe

	pa = rgamma(1,shape = a_prec_shape + nG/2, rate = a_prec_rate + 1/2*sum(a^2))
	y_resid = data$y - Z %*% a
	pe = rgamma(1,shape = resid_prec_shape + n/2, rate = resid_prec_rate + 1/2*sum(y_resid^2))

	if( (i-burn) %% thin == 0 ){
		r = rbind(r,c(1/pa,1/pe))
	}

}
plot(r[,1])
qqplot(m1$VCV[,1],r[,1]);abline(0,1)
qqplot(m1$VCV[,2],r[,2]);abline(0,1)

pas = rgamma(1e5,shape = a_prec_shape,rate = 1/a_prec_rate)
pes = rgamma(1e5,shape = resid_prec_shape,rate = 1/resid_prec_rate)
# h2 = pes/(pes+pas)
# hist(h2,prob=T)
# z = seq(0,1,length=100)
# lines(z,dbeta(seq(0,1,length=h2_divisions+1),2,2)[-(h2_divisions+1)])


h2_divisions = 100
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