new_halfSib_spline_simulation = function(name, nSire,nRep,p, Time, k, k_G, factor_h2s = rep(0.5,k_G),resid_h2 = rep(0.6,p)){
  require(MCMCglmm)
  require(pedantics)
  # build pedigree
  pedigree = data.frame(ind=nSire*nRep + nSire + 1:(nSire*nRep),dam=1:(nSire*nRep) + nSire, sire = gl(nSire,nRep))
  pedigree<-fixPedigree(pedigree)
  children = !is.na(pedigree[,3])

  #generate A matrix as 2* kinship matrix from whole pedigree
  Kinv = forceSymmetric(inverseA(pedigree)$Ainv)
  K = forceSymmetric(solve(Kinv))
  rownames(K) = rownames(Kinv)
  K = K[children,children]
  # K[K==K[1,2]] = 0.9999

  K_chol = chol(K)

  # Lambda matrix
  # factor_h2s = rep(0,k)
  # factor_h2s[2+1:k_G] = runif(k_G)
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
  # resid_h2 = factor_h2s

  # G, R matrices
  G = tcrossprod(sweep(Lambda,2,sqrt(factor_h2s),'*')) + diag(resid_h2/tot_Y_prec)
  R = tcrossprod(sweep(Lambda,2,sqrt(1-factor_h2s),'*')) + diag((1-resid_h2)/tot_Y_prec)

  # fixed effect design
  n = nrow(K)
  X = cbind(1,gl(2,n/2))
  colnames_X = c('intercept','group')
  X_F = X[,-1]

  # RE design
  data = droplevels(data.frame(X,Sire=factor(pedigree$sire[children]),animal = as.factor(pedigree$id[children])))
  Z = diag(1,n)
  # Z = model.matrix(~0+Sire,data)
  # r = ncol(Z)

  b = 2
  B_F = matrix(rnorm((b-1)*k),nrow = (b-1),ncol=k)
  B = rbind(rnorm(p),matrix(0,(b-1),p))

  U_F = t(K_chol) %*% matrix(rnorm(n*k,0,sqrt(factor_h2s)),n,k,byrow=T)

  F = X_F %*% B_F + Z %*% U_F + matrix(rnorm(n*k,0,sqrt(1-factor_h2s)),n,k,byrow=T)

  U_R = t(K_chol) %*% matrix(rnorm(n*p,0,sqrt(resid_h2/tot_Y_prec)),n,p,byrow=T)

  # Eta = X %*% B + F %*% t(Lambda) + U_R + matrix(rnorm(n*p,0,sqrt((1-resid_h2)/tot_Y_prec)),n,p,byrow=T)
  Eta = X %*% B + U_R + matrix(rnorm(n*p,0,sqrt((1-resid_h2)/tot_Y_prec)),n,p,byrow=T)

  # coefficients = splines::bs(Time,df = p)
  coefficients = poly(Time,degree = p)
  observations = c()
  for(i in 1:nrow(Eta)){
    observations = rbind(observations,data.frame(ID = i, covariate = Time, Y = coefficients %*% Eta[i,]))
  }
  setup = list(
    observations = observations,
    Eta = Eta,
    data = data,
    K = K,
    B = B,
    B_F = B_F,
    error_factor_Lambda = Lambda,
    h2 = resid_h2,
    factor_h2s = factor_h2s,
    G = G,
    R = R,
    X = X,
    X_F = X_F,
    U_F = U_F,
    U_R = U_R,
    F = F,
    name = name
  )
  save(setup,file='setup.RData')
}
