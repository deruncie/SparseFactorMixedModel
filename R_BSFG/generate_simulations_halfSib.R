new_halfSib_simulation = function(name, nSire,nRep,p, b, k, k_G, i_Va = 0.2, i_Ve = 0.2){
  require(MCMCglmm)
  require(pedantics)
  # build pedigree
  pedigree = data.frame(ind=nSire*nRep + nSire + 1:(nSire*nRep),dam=1:(nSire*nRep) + nSire, sire = gl(nSire,nRep))
  pedigree<-fixPedigree(pedigree)
  children = !is.na(pedigree[,3])

  #generate A matrix as 2* kinship matrix from whole pedigree
  Ainv = forceSymmetric(inverseA(pedigree)$Ainv)
  A = forceSymmetric(solve(Ainv))
  rownames(A) = rownames(Ainv)
  A = A[children,children]

  A_chol = chol(A)

  # Lambda matrix
  factor_h2s = rep(0,k)
  factor_h2s[2+1:k_G] = runif(k_G)
  Lambda = matrix(0,p,k)
  numeff = sample((p/30):(p/4),k,replace=T)
  for(h in 1:k){
    Lambda[sample(1:p,numeff[h]),h] = rnorm(numeff[h])
  }
  Lambda = Lambda[,order(-diag(t(Lambda) %*% Lambda))]
  cols=1:k
  g_cols = factor_h2s>0
  Lambda = Lambda[do.call("order", unname(split(-abs(Lambda[,cols[factor_h2s>0]]), col(Lambda[,cols[factor_h2s>0]])))),]

  # resid variances
  tot_Y_prec = rep(0.5,p)
  resid_h2 = runif(p,0,0.6)

  # G, R matrices
  G = tcrossprod(sweep(Lambda,2,sqrt(factor_h2s),'*')) + diag(resid_h2/tot_Y_prec)
  R = tcrossprod(sweep(Lambda,2,sqrt(1-factor_h2s),'*')) + diag((1-resid_h2)/tot_Y_prec)

  # fixed effect design
  n = nrow(A)
  X = cbind(1,matrix(rnorm(n*(b-1)),nrow=n))
  colnames_X = 'intercept'
  if(b > 1) colnames_X = c(colnames_X,paste0('Fixed',1:max(1,(b-1))))
  colnames(X) = colnames_X
  X_F = X[,-1]

  # RE design
  data = data.frame(X,Sire=factor(pedigree$sire[children]),animal = factor(pedigree$id[children]))
  Z = diag(1,n)
  # Z = model.matrix(~0+Sire,data)
  # r = ncol(Z)

  B_F = matrix(rnorm((b-1)*k),nrow = (b-1),ncol=k)
  B = rbind(rnorm(p),matrix(0,(b-1),p))

  F_a = A_chol %*% matrix(rnorm(n*k,0,sqrt(factor_h2s)),n,k,byrow=T)

  F = X_F %*% B_F + Z %*% F_a + matrix(rnorm(n*k,0,sqrt(1-factor_h2s)),n,k,byrow=T)

  E_a = A_chol %*% matrix(rnorm(n*p,0,sqrt(resid_h2/tot_Y_prec)),n,p,byrow=T)

  Y = X %*% B + F %*% t(Lambda) + E_a + matrix(rnorm(n*p,0,sqrt((1-resid_h2)/tot_Y_prec)),n,p,byrow=T)
  Y = as.matrix(Y)


  setup = list(
    Y = Y,
    data = data,
    A = A,
    B = B,
    B_F = B_F,
    error_factor_Lambda = Lambda,
    h2 = resid_h2,
    factor_h2s = factor_h2s,
    G = G,
    R = R,
    X = X,
    X_F = X_F,
    name = name
  )
  save(setup,file='setup.RData')
}




prepare_simulation_factorModel = function(A_chol,A,Lambda,factor_h2s,E_a_prec,resid_Y_prec,name){
	n = dim(A_chol)[1]
	p = nrow(Lambda)

	r = n
	# n = 2*n

	k = ncol(Lambda)
	F_a = t(A_chol) %*% sweep(matrix(rnorm(n*k),n),2,sqrt(factor_h2s),'*')
	F_r = sweep(matrix(rnorm(n*k),n),2,sqrt(1-factor_h2s),'*')
	F = F_a + F_r

	E_a_act = t(A_chol) %*% sweep(matrix(rnorm(n*p),n),2,sqrt(1/E_a_prec),'*')
	E_r_act = sweep(matrix(rnorm(n*p),n),2,sqrt(1/resid_Y_prec),'*')

	Y = F %*% t(Lambda) + E_a_act + E_r_act

	X = matrix(1,nrow=1,ncol=n)
	B = rep(0,p)
	Z_1 = diag(r)[,rep(1:r,n/r)]
	Va = colSums(Lambda^2) * factor_h2s + 1/E_a_prec
	Ve = colSums(Lambda^2) * (1-factor_h2s) + 1/resid_Y_prec
	h2 = Va/(Va + Ve)
	G = Lambda %*% diag(factor_h2s) %*% t(Lambda) + diag(1/E_a_prec)
	R = Lambda %*% diag(1-factor_h2s) %*% t(Lambda) + diag(1/resid_Y_prec)
	h2 = diag(G)/(diag(G)+diag(R))

	# Z_2 = rbind(rep(c(1,0,0,0),n/4),rep(c(0,1,0,0),n/4),rep(c(0,0,1,0),n/4),rep(c(0,0,0,1),n/4))

	# plot(h2,type="l")

	print(name)
	setup = list(Y=Y,F_act=F,E_a_act=E_a_act,F_a_act = F_a,Z_1=Z_1,A=A,X=X,n=n,r=r,p=p,G=G,R=R,
		gen_factor_Lambda = Lambda,error_factor_Lambda = Lambda,h2=h2,B=B, factor_h2s=factor_h2s,name=name)

	save(setup,file = 'setup.RData')
	# do.call(writeMat,c("setup.mat",setup))
	save(pedigree,file="pedigree.Robj")
}

prepare_simulation = function(G_chol,R_chol,A_chol,A,gen_factor_Lambda,error_factor_Lambda,factor_h2s,name){
	n = dim(A_chol)[1]
	p = dim(G_chol)[1]

	r = n
	# n = 2*n

	Y=U_act=E_act=c()
	U_act = t(A_chol) %*% matrix(rnorm(p*r),r,p) %*% G_chol
	E_act = matrix(rnorm(p*n),n,p) %*% R_chol
	Y=rbind(U_act) + E_act
	# Y=rbind(U_act,U_act) + E_act


	X = matrix(1,nrow=1,ncol=n)
	B = rep(0,p)
	Z_1 = diag(r)[,rep(1:r,n/r)]
	h2 = diag(G)/(diag(G)+diag(R))

	# Z_2 = rbind(rep(c(1,0,0,0),n/4),rep(c(0,1,0,0),n/4),rep(c(0,0,1,0),n/4),rep(c(0,0,0,1),n/4))

	plot(h2,type="l")

	print(name)
	setup = list(
		Y=Y,U_act=U_act,E_act=E_act,Z_1=Z_1,A=A,X=X,
		n=n,r=r,p=p,
		gen_factor_Lambda = gen_factor_Lambda,error_factor_Lambda = error_factor_Lambda,
		h2=h2,G=G,R=R,B=B, factor_h2s=factor_h2s,name=name
	)
	save(setup,file = 'setup.RData')

	do.call(writeMat,c("setup.mat",setup))
	save(pedigree,file="pedigree.Robj")
}


CovToCor=function(x) diag(1/sqrt(diag(x))) %*% x %*% diag(1/sqrt(diag(x)))


# set seed for repeatability
set.seed(1)

counter = 0

#sample_size simulations
num_reps = 10
p = 100
k=10
factor_h2s = seq(1,0,by=-0.1)[1:k]
factor_h2s[seq(1,k,by=2)]=0

sample_size = list(Small = c(50,5), Medium = c(100,10), Large = c(500,10))


for(size in names(sample_size)){
	nM = sample_size[[size]][1]
	nI = sample_size[[size]][2]
	pedigree = data.frame(ind=nM*nI + nM + 1:(nM*nI),dam=1:(nM*nI) + nM, sire = gl(nM,nI))
	pedigree<-fixPedigree(pedigree)
	children = !is.na(pedigree[,3])

	#generate A matrix as 2* kinship matrix from whole pedigree
	A = 2*kinship(pedigree[,1], pedigree[,2], pedigree[,3])[children,children]
	A_chol = chol(A)
	for(rep in 1:num_reps){
		counter = counter+1
		Lambda = matrix(0,p,k)
		numeff = sample((p/30):(p/4),k,replace=T)
		for(h in 1:k){
			Lambda[sample(1:p,numeff[h]),h] = rnorm(numeff[h])
		}
		Lambda = Lambda[,order(-diag(t(Lambda) %*% Lambda))]
		cols=1:k
		g_cols = factor_h2s>0
		Lambda =Lambda[do.call("order", unname(split(-abs(Lambda[,cols[factor_h2s>0]]), col(Lambda[,cols[factor_h2s>0]])))),]


		E_a_prec = rep(1/0.2,p)
		resid_Y_prec = rep(1/0.2,p)

		name = paste('SampleSize',size,rep,sep='_')
		dir = paste('Sim',counter,sep="_")
		try(dir.create(dir))
		setwd(dir)
		prepare_simulation_factorModel(A_chol,A,Lambda,factor_h2s,E_a_prec,resid_Y_prec,name)
		setwd('..')
	}
}

#k simulations
num_reps = 1
p = 100
factor_h2 = 0.5

ks_size = list(Small = c(5,10), Medium = c(15,25), Large = c(30,50))

nM = 100
nI = 10
pedigree = data.frame(ind=nM*nI + nM + 1:(nM*nI),dam=1:(nM*nI) + nM, sire = gl(nM,nI))
pedigree<-fixPedigree(pedigree)
children = !is.na(pedigree[,3])
#generate A matrix as 2* kinship matrix from whole pedigree
A = 2*kinship(pedigree[,1], pedigree[,2], pedigree[,3])[children,children]
A_chol = chol(A)

for(size in names(ks_size)){
	for(rep in 1:num_reps){
		k = ks_size[[size]][2]
		counter = counter+1
		Lambda = matrix(0,p,k)
		numeff = sample((p/30):(p/4),k,replace=T)
		for(h in 1:k){
			Lambda[sample(1:p,numeff[h]),h] = rnorm(numeff[h])
		}
		Lambda = Lambda[,order(-diag(t(Lambda) %*% Lambda))]
		cols=1:k
		factor_h2s = rep(0,k)
		factor_h2s[1:ks_size[[size]][1]] = factor_h2
		# factor_h2s[sample(1:k,ks_size[[size]][1])] = factor_h2
		# g_cols = factor_h2s>0
		# Lambda =Lambda[do.call("order", unname(split(-abs(Lambda[,cols[factor_h2s>0]]), col(Lambda[,cols[factor_h2s>0]])))),]

		E_a_prec = rep(1/0.2,p)
		resid_Y_prec = rep(1/0.2,p)

		name = paste('k',size,rep,sep='_')
		dir = paste('Sim',counter,sep="_")
		try(dir.create(dir))
		setwd(dir)
		prepare_simulation_factorModel(A_chol,A,Lambda,factor_h2s,E_a_prec,resid_Y_prec,name)
		setwd('..')
	}
}


#unconstrained R
library(MCMCpack)
num_reps = 10
p = 100
kG = 5

nM = 100
nI = 10
pedigree = data.frame(ind=nM*nI + nM + 1:(nM*nI),dam=1:(nM*nI) + nM, sire = gl(nM,nI))
pedigree<-fixPedigree(pedigree)
children = !is.na(pedigree[,3])
#generate A matrix as 2* kinship matrix from whole pedigree
A = 2*kinship(pedigree[,1], pedigree[,2], pedigree[,3])[children,children]
A_chol = chol(A)

	for(rep in 1:num_reps){
		k = kG
		counter = counter+1
		Lambda = matrix(0,p,k)
		numeff = sample((p/30):(p/4),k,replace=T)
		for(h in 1:k){
			Lambda[sample(1:p,numeff[h]),h] = rnorm(numeff[h])
		}
		Lambda = Lambda[,order(-diag(t(Lambda) %*% Lambda))]
		cols=1:k
		factor_h2s = rep(1,k)
		g_cols = factor_h2s>0
		Lambda =Lambda[do.call("order", unname(split(-abs(Lambda[,cols[factor_h2s>0]]), col(Lambda[,cols[factor_h2s>0]])))),]


		G = Lambda %*% diag(factor_h2s) %*% t(Lambda) + 0.2*diag(p)
		R = rwish(p+1,diag(diag(G))/p)
		G_chol = chol(G)
		R_chol = chol(R)
		h2 = diag(G)/(diag(G)+diag(R))
		hist(h2,breaks=seq(0,1,length=30))
		image(CovToCor(G)[,p:1])
		image(CovToCor(R)[,p:1])

		gen_factor_Lambda = (Lambda %*% diag(sqrt(factor_h2s)))[,factor_h2s>0]
		error_factor_Lambda = Lambda

		name = paste('Rwish',rep,sep='_')
		dir = paste('Sim',counter,sep="_")
		try(dir.create(dir))
		setwd(dir)
		prepare_simulation(G_chol,R_chol,A_chol,A,gen_factor_Lambda,error_factor_Lambda,factor_h2s,name)
		setwd('..')
	}



#non-sparse R
num_reps = 10
p = 100
kG = 5
kR = 5

nM = 100
nI = 10
pedigree = data.frame(ind=nM*nI + nM + 1:(nM*nI),dam=1:(nM*nI) + nM, sire = gl(nM,nI))
pedigree<-fixPedigree(pedigree)
children = !is.na(pedigree[,3])
#generate A matrix as 2* kinship matrix from whole pedigree
A = 2*kinship(pedigree[,1], pedigree[,2], pedigree[,3])[children,children]
A_chol = chol(A)

	for(rep in 1:num_reps){
		k = kG
		counter = counter+1
		Lambda = matrix(0,p,k)
		numeff = sample((p/30):(p/4),k,replace=T)
		for(h in 1:k){
			Lambda[sample(1:p,numeff[h]),h] = rnorm(numeff[h])
		}
		Lambda = Lambda[,order(-diag(t(Lambda) %*% Lambda))]
		cols=1:k
		factor_h2s = rep(0.5,k)
		g_cols = factor_h2s>0
		Lambda =Lambda[do.call("order", unname(split(-abs(Lambda[,cols[factor_h2s>0]]), col(Lambda[,cols[factor_h2s>0]])))),]

		Lambda = cbind(Lambda,matrix(rnorm(kR*p),p,kR))
		factor_h2s = c(factor_h2s,rep(0,kR))


		G = Lambda %*% diag(factor_h2s) %*% t(Lambda) + 0.2*diag(p)
		R = Lambda %*% diag(1-factor_h2s) %*% t(Lambda) + 0.2*diag(p)
		G_chol = chol(G)
		R_chol = chol(R)
		h2 = diag(G)/(diag(G)+diag(R))
		hist(h2,breaks=seq(0,1,length=30))
		image(CovToCor(G)[,p:1])
		image(CovToCor(R)[,p:1])

		gen_factor_Lambda = (Lambda %*% diag(sqrt(factor_h2s)))[,factor_h2s>0]
		error_factor_Lambda = Lambda

		name = paste('R_notSparse',rep,sep='_')
		dir = paste('Sim',counter,sep="_")
		try(dir.create(dir))
		setwd(dir)
		prepare_simulation(G_chol,R_chol,A_chol,A,gen_factor_Lambda,error_factor_Lambda,factor_h2s,name)
		setwd('..')
	}



#p simulations
num_reps = 10
factor_h2 = 0.5
k = 10
k_G = 5

ks_size = list(Small = 20, Large = 1000)

nM = 100
nI = 10
pedigree = data.frame(ind=nM*nI + nM + 1:(nM*nI),dam=1:(nM*nI) + nM, sire = gl(nM,nI))
pedigree<-fixPedigree(pedigree)
children = !is.na(pedigree[,3])
#generate A matrix as 2* kinship matrix from whole pedigree
A = 2*kinship(pedigree[,1], pedigree[,2], pedigree[,3])[children,children]
A_chol = chol(A)

for(size in names(ks_size)){
	for(rep in 1:num_reps){
		p = ks_size[[size]]
		counter = counter+1
		Lambda = matrix(0,p,k)
		numeff = sample((p/30):(p/4),k,replace=T)
		numeff[numeff<3]=3
		for(h in 1:k){
			Lambda[sample(1:p,numeff[h]),h] = rnorm(numeff[h])
		}
		Lambda = Lambda[,order(-diag(t(Lambda) %*% Lambda))]
		cols=1:k
		factor_h2s = rep(0,k)
		factor_h2s[1:k_G] = factor_h2
		# factor_h2s[sample(1:k,ks_size[[size]][1])] = factor_h2
		# g_cols = factor_h2s>0
		# Lambda =Lambda[do.call("order", unname(split(-abs(Lambda[,cols[factor_h2s>0]]), col(Lambda[,cols[factor_h2s>0]])))),]

		E_a_prec = rep(1/0.2,p)
		resid_Y_prec = rep(1/0.2,p)

		gen_factor_Lambda = (Lambda %*% diag(sqrt(factor_h2s)))[,factor_h2s>0]
		error_factor_Lambda = Lambda

		name = paste('p',size,rep,sep='_')
		dir = paste('Sim',counter,sep="_")
		try(dir.create(dir))
		setwd(dir)
		prepare_simulation_factorModel(A_chol,A,Lambda,factor_h2s,E_a_prec,resid_Y_prec,name)
		setwd('..')
	}
}

