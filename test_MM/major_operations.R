# These are the major types of problems involved in the collapsed Gibbs sampler:


library(MCMCglmm)
library(Matrix)
library(pedantics)
library(rbenchmark)

# fix_A = function(x) {
# 	a = as.matrix(x)
# 	a[abs(a)<1e-10] = 0
# 	Matrix(a)
# }
fix_A = function(x) forceSymmetric(drop0(x,tol = 1e-10))

# types of A matrices:

# type 1: diagonal A, bigZ
nG = 50
nR = 5
data = data.frame(Group = gl(nG,nR))
A1 = Diagonal(nG,1)
A1_inv = solve(A1)
Z1 = Matrix(model.matrix(~0+Group,data))
ZtZ_1 = crossprod(Z1)
ZAZ_1 = forceSymmetric(Z1 %*% A1 %*% t(Z1))

# type 2: big sparse A
nM = 50
nI = 5
pedigree = data.frame(ind=nM*nI + nM + 1:(nM*nI),dam=1:(nM*nI) + nM, sire = gl(nM,nI))
pedigree<-fixPedigree(pedigree)
children = !is.na(pedigree[,3])

#generate A matrix as 2* kinship matrix from whole pedigree
A2 = Matrix(2*kinship(pedigree[,1], pedigree[,2], pedigree[,3])[children,children])
A2_inv = forceSymmetric(solve(A2))

Z2 = Diagonal(nM*nI,1)
ZtZ_2 = crossprod(Z2)
ZAZ_2 = Z2 %*% A2 %*% t(Z2)

# type 3: sparse Ainv
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

data = droplevels(subset(Data,!is.na(BWT)))
data$animal = factor(as.character(data$animal),levels = data$animal)
# data$animal = factor(as.character(data$animal),levels = Ped$ID)

index = match(data$animal,Ped[,1])


A3_inv_full = forceSymmetric(inverseA(Ped)[[1]])
A3 = solve(A3_inv_full)[index,index]
A3 = fix_A(A3)
A3_inv = solve(A3)
A3_inv = fix_A(A3_inv)
Z3 = Diagonal(nrow(A3),1)
ZtZ_3 = crossprod(Z3)

ZAZ_3 = Z3 %*% A3 %*% t(Z3)

# constructed matrices:
# C = ZtZ + Ai
# Sigma = ZAZ + I
# chol_Sigma_inv

ZAZ = ZAZ_1
Ai = A1_inv
ZtZ = ZtZ_1

ZAZ = ZAZ_2
Ai = A2_inv
ZtZ = ZtZ_2

ZAZ = ZAZ_3
Ai = A3_inv
ZtZ = ZtZ_3

n = nrow(ZAZ)

C = ZtZ + Ai
chol_C = chol(C)
Cholesky_C = Cholesky(C)
Sigma = ZAZ + Diagonal(nrow(ZAZ))
chol_Sigma = chol(Sigma)
chol_Sigma_inv = chol(solve(Sigma))
chol_Sigma_inv2 = t(solve(chol_Sigma))
Cholesky_Sigma = Cholesky(Sigma,LDL=F,perm=F)  # note: This is LLt = A (not LtL!)
Cholesky_Sigma2 = Cholesky(Sigma,LDL=F)  # note: This is LLt = A (not LtL!)
chol_Sigma_inv3 = solve(Cholesky_Sigma,Diagonal(nrow(Sigma)),'L')

# testing
Sigma_inv = forceSymmetric(solve(Sigma))
all.equal(Sigma_inv,crossprod(chol_Sigma_inv),tol = 1e-14)
all.equal(Sigma_inv,crossprod(chol_Sigma_inv2),tol = 1e-14)
all.equal(Sigma_inv,crossprod(chol_Sigma_inv3),tol = 1e-14)

# calculations:

# from sample_MME_fixedEffects_inv
b = 30
W = matrix(rnorm(n*b),nrow=n)
prior_prec = runif(b,0,10)
tot_Y_prec = 10

# all.equal(chol_Sigma_inv*sqrt(tot_Y_prec),chol(Sigma_inv*tot_Y_prec))
Wp = solve(Cholesky_Sigma2,W,'P')
benchmark(
RinvSqW1 <<- sqrt(tot_Y_prec)*solve(t(chol_Sigma),W),
RinvSqW2 <<- sqrt(tot_Y_prec)*chol_Sigma_inv %*% W,
RinvSqW2b <<- sqrt(tot_Y_prec)*chol_Sigma_inv2 %*% W,
RinvSqW2c <<- sqrt(tot_Y_prec)*chol_Sigma_inv3 %*% W,
RinvSqW3 <<- sqrt(tot_Y_prec)*solve(Cholesky_Sigma,W,'L'),
# RinvSqW3b <<- sqrt(tot_Y_prec)*solve(Cholesky_Sigma2,solve(Cholesky_Sigma2,solve(Cholesky_Sigma2,W,'P'),'L'),'Pt'),
RinvSqW3b <<- sqrt(tot_Y_prec)*solve(Cholesky_Sigma2,solve(Cholesky_Sigma2,Wp,'L'),'Pt'),
RinvSqW3c <- sqrt(tot_Y_prec)*solve(Cholesky_Sigma2,Wp,'L'),
replications = 100)

# A1: 2's, 3 best. 3 1.5x, 1/3b 3x
# A2: 2's best. 3 1.5x, 1 2.5x, 3b 3.5x
# A3: 3, 4 ~3-4x faster than rest

# conclusion: 3b much slower if A is simple, but much faster if A is complex. 
# Do I need to un-permute? Then 3c is not much slower even if A is simple. 
#    not for Cw below.

all.equal(crossprod(RinvSqW1),crossprod(RinvSqW2),tol = 1e-14)
all.equal(crossprod(RinvSqW1),crossprod(RinvSqW2b),tol = 1e-14)
all.equal(crossprod(RinvSqW1),crossprod(RinvSqW2c),tol = 1e-14)
all.equal(crossprod(RinvSqW1),crossprod(RinvSqW3),tol = 1e-14)
all.equal(crossprod(RinvSqW1),crossprod(RinvSqW3c),tol = 1e-14)


Cw = crossprod(RinvSqW1)
Cw@x = Cw@x + prior_prec

# from sample_MME_single_diagA_inv
e_star = matrix(rnorm(n*100),n)
y = rnorm(n)

cc1 = expand(Cholesky_Sigma)$L
cc2 = t(expand(Cholesky_Sigma2)$P) %*% expand(Cholesky_Sigma2)$L 
# microbenchmark(
benchmark(
e1 <<- 1/sqrt(tot_Y_prec)*t(chol_Sigma) %*% e_star,
e2 <<- 1/sqrt(tot_Y_prec)*solve(chol_Sigma_inv,e_star),
e2b <<- 1/sqrt(tot_Y_prec)*solve(chol_Sigma_inv2,e_star),
e2c <<- 1/sqrt(tot_Y_prec)*solve(chol_Sigma_inv3,e_star),
e3 <<- 1/sqrt(tot_Y_prec)*cc1 %*% e_star,
e3b <<- 1/sqrt(tot_Y_prec)*cc2 %*% e_star,
replications = 100)
# times = 10)

# A1: 3b ~1.5x others
# A2: 2's best. 1, 3s ~1.5x
# A3: 3b ~3x faster than 3, 1. 2's another 1.5x slower than 3/1

# conclusion: use 3b. It's rarely worse, and sometimes much better

i = sample(1:prod(dim(e1)),10000)
plot(c(cov(as.matrix(t(e1))))[i],c(as.matrix(Sigma/tot_Y_prec))[i]);abline(0,1)
plot(c(cov(as.matrix(t(e2))))[i],c(as.matrix(Sigma/tot_Y_prec))[i]);abline(0,1)
plot(c(cov(as.matrix(t(e2b))))[i],c(as.matrix(Sigma/tot_Y_prec))[i]);abline(0,1)
plot(c(cov(as.matrix(t(e2c))))[i],c(as.matrix(Sigma/tot_Y_prec))[i]);abline(0,1)
plot(c(cov(as.matrix(t(e3))))[i],c(as.matrix(Sigma/tot_Y_prec))[i]);abline(0,1)
plot(c(cov(as.matrix(t(e3b))))[i],c(as.matrix(Sigma/tot_Y_prec))[i]);abline(0,1)


P = as(Cholesky_Sigma2,'pMatrix')
Py = P %*% y
benchmark(
e1 <<- solve(t(chol_Sigma),y) * sqrt(tot_Y_prec),
e2 <<- sqrt(tot_Y_prec)*chol_Sigma_inv %*% y,
e2b <<- sqrt(tot_Y_prec)*chol_Sigma_inv2 %*% y,
e2c <<- sqrt(tot_Y_prec)*chol_Sigma_inv3 %*% y,
e3 <<- solve(Cholesky_Sigma,y,'L') * sqrt(tot_Y_prec),
e3b <<- solve(Cholesky_Sigma2,solve(Cholesky_Sigma2,solve(Cholesky_Sigma2,y,'P'),'L'),'Pt') * sqrt(tot_Y_prec),
# e3c <<- solve(Cholesky_Sigma2,solve(Cholesky_Sigma2,y,'P'),'L') * sqrt(tot_Y_prec),
e3c <<- solve(Cholesky_Sigma2,Py,'L') * sqrt(tot_Y_prec),
replications = 100)


# A1: 2's best. 3/1 ~1.6x, 3b 3.5x
# A2: 2's best. 3 ~1.6x, 1 2x, 3b 3.5x
# A3: 2's best. 3's ~1.8x, 1 ~5x

# use 3c. it is compatible with RinvSqW3c above, which is much faster. Plus if you pre-calculate P and Py, it's only marginally slower than 2's

all.equal(crossprod(RinvSqW1,e1),crossprod(RinvSqW2,e2))
all.equal(crossprod(RinvSqW2,e2),crossprod(RinvSqW2,e2))
all.equal(crossprod(RinvSqW2b,e2b),crossprod(RinvSqW2,e2))
all.equal(crossprod(RinvSqW2c,e2c),crossprod(RinvSqW2,e2))
all.equal(crossprod(RinvSqW3,e3),crossprod(RinvSqW2,e2))
all.equal(crossprod(RinvSqW3b,e3b),crossprod(RinvSqW2,e2))
all.equal(crossprod(RinvSqW3c,e3c),crossprod(RinvSqW2,e2))

WtRinvy = crossprod(RinvSqW1,e1)
theta1 = solve(Cw,WtRinvy)  # no options here

# from sample_h2s_discrete_inv
chol_Sigma_inv3c = solve(Cholesky_Sigma2,Diagonal(nrow(Sigma)),'L')
Y = matrix(rnorm(n*100),n,100)
Yp = solve(Cholesky_Sigma2,Y,'P')
benchmark(
s2_1 <<- colSums(solve(t(chol_Sigma),Y)^2),
s2_2 <<- colSums((chol_Sigma_inv %*% Y)^2),
s2_2b <<- colSums((chol_Sigma_inv2 %*% Y)^2),
s2_2c <<- colSums((chol_Sigma_inv3 %*% Y)^2),
s2_3 <<- colSums(solve(Cholesky_Sigma,Y,'L')^2),
s2_3b <<- colSums(solve(Cholesky_Sigma2,Yp,'L')^2),
s2_3c <<- colSums((chol_Sigma_inv3c %*% Yp)^2),
s2_4 <<- diag(t(Y) %*% Sigma_inv %*% Y),
replications = 100)


# A1: 2's 1.2x, 3 1.5x, 1/3b 2.5x, 4 much worse if Y has multiple columns
# A2: 2's 1.2x, 1/3 1.5x, 3b 2.5x, 4 much worse if Y has multiple columns
# A3: 2's, 3b best. 4 1.5x, 3 2x, 1 4x,4 much worse if Y has multiple columns

# conclusion: use 2 (any)
# conclusion:  actually, use 3b with pre-calculated Yp (test if needed first)

all.equal(s2_1,s2_2,s2_2b,s2_2c,s2_3,s2_3b,s2_4,s2_3c)

# from sample_MME_single_diagR
w = rnorm(nrow(C))
benchmark(
theta1 <<- solve(C,w)/tot_Y_prec,
theta2 <<- solve(Cholesky_C,w)/tot_Y_prec,
theta3 <<- solve((chol_C),solve(t(chol_C),w)) / tot_Y_prec,
replications = 100)

# A1: all about the same
# A2: all about the same
# A3: 2 best, 1 5x, 3 7x

# conclusion: use 2

all.equal(theta1,theta2,tol = 1e-14)
all.equal(theta1,theta3,tol = 1e-14)




