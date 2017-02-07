## 1 - my code doesn't seem to give the same variance component estimates. Maybe the priors are different
## 2 - I should try to re-write the functions in terms of Ainv rather than A. It is much more sparse at least in this example.

library(MCMCglmm)
Data <- as.data.frame(read.table(file = "./gryphon.dat", header = TRUE))
names(Data)[1] <- "animal"
Data$animal <- as.factor(Data$animal)
Data$MOTHER <- as.factor(Data$MOTHER)
Data$BYEAR <- as.factor(Data$BYEAR)
Data$SEX <- as.factor(Data$SEX)
Data$BWT <- as.numeric(Data$BWT)
Data$TARSUS <- as.numeric(Data$TARSUS)
head(Data)
Ped <- as.data.frame(read.table(file = "./gryphon.ped", header = TRUE))
for (x in 1:3) Ped[, x] <- as.factor(Ped[, x])
head(Ped)

# Data$BWT = matrix(rnorm(nrow(Data),0,sd(Data$BWT,na.rm=T)),nc=1)

prior1.1 <- list(G = list(G1 = list(V = 1, nu = 0.002)), R = list(V = 1,nu = 0.002))
prior1.1 = list(R = list(V = 1/3, nu = 3), G = list(G1 = list(V = 1/3,nu=3)))
model1.1 <- MCMCglmm(BWT ~ 1, random = ~animal, pedigree = Ped, data = Data, prior = prior1.1)

model1.2 <- MCMCglmm(BWT ~ SEX, random = ~animal, pedigree = Ped,
    data = Data, prior = prior1.1, nitt = 6500, thin = 5, burnin = 1500,
    verbose = FALSE)

prior1.3 <- list(G = list(G1 = list(V = 1, n = 0.002), G2 = list(V = 1,n = 0.002)), R = list(V = 1, n = 0.002))
model1.3 <- MCMCglmm(BWT ~ SEX, random = ~animal + BYEAR, pedigree = Ped,
    data = Data, nitt = 65000, thin = 50, burnin = 15000, prior = prior1.3,
    verbose = FALSE)
posterior.mode(model1.3$VCV)


# now BSFG
# sourceCpp(paste(model_path,'BSFG_discreteRandom_functions_c.cpp',sep='/')) #,cacheDir = model_path

# # testing methods to solve from A for only subset of individuals
# index = match(data$animal,rownames(model1.2$ginverse$animal))

# Ainv_full = forceSymmetric(model1.2$ginverse$animal)
# cAinv_full = chol(Ainv_full)
# A_full = solve(Ainv_full)
# A = A_full[index,index]
# cA = chol(A)
# cAinv = chol(solve(A))
# cAinv2t = solve(cA)

# i = sample(1:prod(dim(A)),10000)
# e1 = t(cA) %*% matrix(rnorm(nrow(A)*1000),nc=1000)
# e2 = solve(cAinv,matrix(rnorm(nrow(A)*1000),nc=1000))
# e3 = solve(t(cAinv2t),matrix(rnorm(nrow(A)*1000),nc=1000))
# e4 = solve(cAinv_full,matrix(rnorm(nrow(Ainv_full)*1000),nc=1000))[index,]
# plot(as.matrix(A)[i],cov(t(as.matrix(e1)))[i]);abline(0,1)
# plot(as.matrix(A)[i],cov(t(as.matrix(e2)))[i]);abline(0,1)
# plot(as.matrix(A)[i],cov(t(as.matrix(e3)))[i]);abline(0,1)
# plot(as.matrix(A)[i],cov(t(as.matrix(e4)))[i]);abline(0,1)


# library(microbenchmark)
# j1 = Matrix(matrix(rnorm(nrow(A)),nc=1))
# j2 = Matrix(matrix(rnorm(nrow(Ainv_full)),nc=1))
# microbenchmark(
# 	t(cA) %*% j1,
# 	solve(cAinv,j1),
# 	solve(t(cAinv2t),j1),
# 	solve(cAinv_full,j2)[index]
# 	)
# # best two options are: 
# 	# solve(cAinv,j1) - calc cAinv from A for target individuals, with A from Ainv_full
# 	# solve(cAinv_full,j2)[index] - go straight from cAinv_full, only use the target individuals


discrete_divisions = 100
priors = list(
    tot_Y_prec_shape      =   2.4,
    tot_Y_prec_rate       =   1/3.5
)

data = droplevels(subset(Data,!is.na(BWT)))
data$animal = factor(as.character(data$animal),levels = data$animal)
# data$animal = factor(as.character(data$animal),levels = Ped$ID)

index = match(data$animal,rownames(model1.2$ginverse$animal))


randomEffects = list(animal = model1.2$ginverse$animal)
randomEffects = lapply(randomEffects,forceSymmetric)

# build Z matrices from random model
RE_names = names(randomEffects)
n_RE = length(randomEffects)
r_RE = sapply(RE_names,function(x) dim(randomEffects[[x]])[1])
Z_matrices = lapply(RE_names,function(re) {
	Diagonal(r_RE[re],1)
	Z = Matrix(model.matrix(formula(sprintf('~0 + %s',re)),data),sparse = TRUE)
	Z[,paste0(re,levels(data[[re]]))]
})
names(Z_matrices) = RE_names

Z_all = do.call(cbind,Z_matrices)

for(re in RE_names[RE_names %in% names(randomEffects) == F]){
	randomEffects[[re]] = Diagonal(ncol(Z_matrices[[re]]),1)
}

# fix_A = function(x) x
fix_A = function(x) {
	a = as.matrix(x)
	a[abs(a)<1e-10] = 0
	Matrix(a)
}
A_mats = lapply(RE_names,function(re) {
	A = solve(randomEffects[[re]])
	index = match(sub(re,'',colnames(Z_matrices[[re]])),rownames(randomEffects[[re]]))
	fix_A(A[index,index])
})
names(A_mats) = RE_names
Ai_mats = lapply(A_mats,function(x) fix_A(solve(x)))

h2_divisions = expand.grid(lapply(RE_names,function(re) 0:discrete_divisions)) / discrete_divisions
colnames(h2_divisions) = RE_names
h2_divisions = t(h2_divisions[rowSums(h2_divisions) < 1,,drop=FALSE])

library(MCMCpack)
priors$discrete_priors = sapply(1:ncol(h2_divisions),function(x) {
    h2s = h2_divisions[,x]
    pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10)
})/100

# chol_Ainvs = lapply(1:n_RE,function(i) chol(Ai_mats[[i]]))

ZAZts = list()
for(i in 1:n_RE){
	# ZAZts[[i]] = crossprod(solve(t(chol_Ainvs[[i]]),t(Z_matrices[[i]])))
	# ZAZts[[i]] = Z_matrices[[i]] %*% solve(Ai_mats[[i]]) %*% t(Z_matrices[[i]])
	# ZAZts[[i]] = Diagonal(length(index),1) %*% solve(Ai_mats[[i]])[index,index] %*% Diagonal(length(index),1)
	ZAZts[[i]] = Z_matrices[[i]] %*% A_mats[[i]] %*% t(Z_matrices[[i]])
}
# ZAZts = list(animal = solve(forceSymmetric(model1.2$ginverse$animal))[index,index])

ZtZ = crossprod(Z_all)
randomEffect_Cs = lapply(1:ncol(h2_divisions),function(i) {    	
	h2s = h2_divisions[,i]
	Ai = do.call(bdiag,lapply(1:length(h2s),function(j) {
			if(h2s[j] == 0) {  # if h2==0, then we want a Diagonal matrix with Inf diagonal. This will allow Cinv = 0
				Diagonal(nrow(Ai_mats[[j]]),Inf)
			} else{
				Ai_mats[[j]]/h2s[j]  
			}
		}))
	C = ZtZ/(1-sum(h2s))
	C = C + Ai
	C
})

make_Sigma = function(ZAZts,h2s){
	R = 0
	for(i in 1:length(h2s)){
		R = R + h2s[i]*ZAZts[[i]]
	}
	R + (1-sum(h2s)) * Diagonal(nrow(R))
}

Sigma_invs = lapply(1:ncol(h2_divisions),function(i) {
	print(i)
	Sigma = make_Sigma(ZAZts,h2_divisions[,i])
	# microbenchmark({
	chol_Sigma = chol(Sigma)
	chol_Sigma_inv = t(solve(chol_Sigma))
	det = det(chol_Sigma_inv)^2
	list(det = det,chol = chol_Sigma_inv)
# 	Si = crossprod(cSi)
# },{
# 	Sigma_inv = solve(Sigma)
# 	chol_Sigma_inv = chol(Sigma_inv)
# },times=10)
	# det = det(chol_Sigma_inv)^2
	# print(c(length(Sigma@x),length(Sigma_inv@x)))
	# list(Sigma_inv = Sigma_inv, det = det,chol = chol_Sigma_inv)
})

# i=2
# Sigma_base = forceSymmetric(make_Sigma(ZAZts,h2_divisions[,i]))
# chol_Base = Cholesky(Sigma_base,perm=F,LDL=F)
# cb = update(chol_Base,Sigma_base)
# all.equal(chol_Base,cb,tol = 1e-14)

# i=20
# Sigma20 = forceSymmetric(make_Sigma(ZAZts,h2_divisions[,i]))
# d = Diagonal(nrow(Sigma20),1)
# microbenchmark({  # best if only need chol_Sigma_inv
# 	chol_Sigma20 = chol(Sigma)  
# 	chol_Sigma_inv20 = t(solve(chol_Sigma20))
# 	Sigma_inv = crossprod(chol_Sigma_inv20)
# 	},{
# Sigma_inv20 = solve(Sigma20) # generally worst
# chol_Sigma_inv20 = chol(Sigma_inv20)
# },{
# cb20 = update(chol_Base,Sigma20) # this works to update chol_Sigma_inv20, but is not as fast as the first one
# chol_Sigma_inv20 = solve(cs20,d,'Lt')
# Sigma_inv20 = solve(cs20)  # adding this and it is slightly faster than 2 and 1.
# },times=10)
# max(abs(si20 - Sigma_inv20))
# max(abs(crossprod(csi_20)-Sigma_inv20))
# max(abs(crossprod(chol_Sigma_inv20)-Sigma_inv20))


# all.equal(chol_Sigma_inv20,csi20,tol = 1e-14)

# max(abs(solve(chol_20) - solve(cb20)))




# Sigma_inv_base = solve(Sigma_base)
# i=20
# Sigma20 = make_Sigma(ZAZts,h2_divisions[,i])
# Sigma_inv20 = solve(Sigma20)
# chol_Sigma_base = Cholesky(Sigma_base,LDL=F)
# cs20 = update(chol_Sigma_base,Sigma20)
# max(abs(solve(chol_Sigma_base) - Sigma_inv_base))
# max(abs(solve(cs20) - Sigma_inv20))

Y = data[,'BWT',drop=F]
fixed = formula(~SEX)

# build X from fixed model
X = model.matrix(fixed,data)
b = ncol(X)
p=1

# animal = solve(model1.2$ginverse$animal)
# index = match(data$animal,rownames(model1.2$ginverse$animal))
# animal = animal[index,index]
# animal[upper.tri(animal)] = t(animal)[upper.tri(animal)]
# animal = as.matrix(animal)
# animal[abs(animal) < 1e-10] = 0
# randomEffects = list(animal = animal)

# # build Z matrices from random model
# RE_names = names(randomEffects)
# n_RE = length(randomEffects)
# r_RE = sapply(RE_names,function(x) dim(randomEffects[[x]])[1])
# Z_matrices = lapply(RE_names,function(re) {
# 	Diagonal(r_RE[re],1)
# 	# Matrix(model.matrix(formula(sprintf('~0 + %s',re)),data),sparse = TRUE)
# })
# names(Z_matrices) = RE_names

# randomEffects = lapply(randomEffects,function(x) Matrix(x))
# names(randomEffects) = RE_names

# ZAZts = list()
# for(i in 1:n_RE){
# 	ZAZts[[i]] = tcrossprod(Z_matrices[[i]] %*% randomEffects[[i]],Z_matrices[[i]])
# }

# make_Sigma = function(ZAZts,h2s){
# 	R = 0
# 	for(i in 1:length(h2s)){
# 		R = R + h2s[i]*ZAZts[[i]]
# 	}
# 	R + (1-sum(h2s)) * Diagonal(nrow(R))
# }
# Sigma1 = make_Sigma(ZAZts,rep(1,length(ZAZts)))

# chol_Sigma1 = Cholesky(Sigma1,LDL=F)



# data = data[apply(expand(chol_Sigma1)$P==1,1,which),]

# animal = solve(model1.2$ginverse$animal)
# index = match(data$animal,rownames(model1.2$ginverse$animal))
# animal = forceSymmetric(animal[index,index])
# animal[animal<0] = 0
# a1 = Matrix(forceSymmetric(animal),sparse=T)
# a2 = Matrix(forceSymmetric(solve(model1.2$ginverse$animal[index,index])),sparse=T)

# ca1 = chol(a1)
# ca1[ca1<0]=0
# ca1 = Cholesky(a1,LDL=F)
# ca2 = Cholesky(a2,LDL=F)
# # animal[upper.tri(animal)] = t(animal)[upper.tri(animal)]
# # animal = as.matrix(animal)
# # animal[abs(animal) < 1e-10] = 0
# randomEffects = list(animal = animal)

# # build Z matrices from random model
# RE_names = names(randomEffects)
# n_RE = length(randomEffects)
# r_RE = sapply(RE_names,function(x) dim(randomEffects[[x]])[1])
# Z_matrices = lapply(RE_names,function(re) {
# 	Matrix(model.matrix(formula(sprintf('~0 + %s',re)),data),sparse = TRUE)
# })
# names(Z_matrices) = RE_names

# randomEffects = lapply(randomEffects,function(x) Matrix(x))
# names(randomEffects) = RE_names

# ZAZts = list()
# for(i in 1:n_RE){
# 	ZAZts[[i]] = tcrossprod(Z_matrices[[i]] %*% randomEffects[[i]],Z_matrices[[i]])
# }

# Z_all = do.call(cbind,Z_matrices)

# Y = data[,'BWT',drop=F]
# fixed = formula(~SEX)

# # build X from fixed model
# X = model.matrix(fixed,data)
# b = ncol(X)
# p=1



tot_Y_prec = with(priors,rgamma(1,shape = tot_Y_prec_shape,rate = tot_Y_prec_rate))

h2_index = sample(1:ncol(h2_divisions),1,replace=T)
h2 = h2_divisions[,h2_index,drop=FALSE]


a = do.call(rbind,lapply(RE_names,function(effect){
	rnorm(r_RE[effect], 0, sqrt(h2[effect,] / tot_Y_prec))
}))

B = rnorm(b)

# chol_As = lapply(1:n_RE,function(i) chol(randomEffects[[i]]))
# Ai_mats = lapply(1:n_RE,function(i) chol2inv(chol_As[[i]]))

# ZtZ = crossprod(Z_all)
# randomEffect_Cs = lapply(1:ncol(h2_divisions),function(i) {    	
# 	h2s = h2_divisions[,i]
# 	Ai = do.call(bdiag,lapply(1:length(h2s),function(i) {
# 			if(h2s[i] == 0) {  # if h2==0, then we want a Diagonal matrix with Inf diagonal. This will allow Cinv = 0
# 				Diagonal(nrow(Ai_mats[[i]]),Inf)
# 			} else{
# 				Ai_mats[[i]]/h2s[i]  
# 			}
# 		}))
# 	C = ZtZ/(1-sum(h2s))
# 	C = C + Ai
# 	C
# })

# make_Sigma = function(ZAZts,h2s){
# 	R = 0
# 	for(i in 1:length(h2s)){
# 		R = R + h2s[i]*ZAZts[[i]]
# 	}
# 	R + (1-sum(h2s)) * Diagonal(nrow(R))
# }

# Sigmas = lapply(1:ncol(h2_divisions),function(i) {
# 	Sigma = make_Sigma(ZAZts,h2_divisions[,i])
# 	det = det(Sigma)
# 	chol = chol(Sigma)
# 	list(Sigma = Sigma, det = det,chol = chol)
# })

nIter = 1000
burn = 50
thin = 1

posterior = matrix(NA,nrow = (nIter-burn)/thin,ncol = length(B)+length(h2)+1)
sp = 0

for(i in 1:nIter){
	if(i %% 100 == 0) print(i)

	Design = X
	rows = b
	prior_mean = matrix(0,rows,p)
	prior_prec = matrix(rep(1e-10,b),nc=1)
	# recover()
	B = sample_MME_fixedEffects_inv(Y,Design,Sigma_invs, h2_index, tot_Y_prec, prior_mean, prior_prec,1)

	Y_tilde = as.matrix(Y - X %*% B)
	tot_Y_prec = sample_tot_prec_inv(Y_tilde, priors$tot_Y_prec_shape, priors$tot_Y_prec_rate, Sigma_invs, h2_index,1)
	h2_index = sample_h2s_discrete_inv(Y_tilde,tot_Y_prec, Sigma_invs,priors$discrete_priors,1)
	h2 = h2_divisions[,h2_index,drop=FALSE]
	# a_prec = tot_Y_prec / colSums(h2)

	# prior_mean = matrix(0,sum(r_RE),p)
	# a = sample_MME_ZAZts(Y_tilde, Z_all, tot_Y_prec, prior_mean, randomEffect_Cs, Ai_mats, h2, h2_index,chol_As,1)	

	if((i > burn) & (i-burn)%%thin == 0){
		sp = sp + 1
		posterior[sp,] = c(B, h2/tot_Y_prec, (1-sum(h2))/tot_Y_prec)
	}

}

boxplot(posterior)
boxplot(cbind(model1.2$Sol,model1.2$VCV))
i=2;qqplot(posterior[,i],model1.2$Sol[,i]);abline(0,1)
i=1;qqplot(posterior[,2+i],model1.2$VCV[,i]);abline(0,1)
i=2;qqplot(posterior[,2+i],model1.2$VCV[,i]);abline(0,1)

