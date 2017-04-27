new_halfSib_simulation_eQTL = function(name, nSire,nRep,p, b, factor_h2s, Va = 0.2, Ve = 0.2,Vb = 0,V_cis,nSNP,bSNP = 1){
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
  # K[K>0 & K<1] = 0.5

  K_chol = chol(K)

  # Lambda matrix
  # factor_h2s = rep(0,k)
  # factor_h2s[2+1:k_G] = runif(k_G)
  k = length(factor_h2s)
  Lambda = matrix(0,p,k)
  numeff = sample((p/30):(p/4),k,replace=T)
  for(h in 1:k){
    Lambda[sample(1:p,numeff[h]),h] = rnorm(numeff[h])
  }
  Lambda = Lambda[,order(-diag(t(Lambda) %*% Lambda))]
  cols=1:k
  g_cols = factor_h2s>0
  if(sum(g_cols) == 0) g_cols = 1:k
  Lambda = Lambda[do.call("order", unname(split(-abs(Lambda[,cols[g_cols]]), col(Lambda[,cols[g_cols]])))),]

  # resid variances
  # tot_Y_prec = rep(0.5,p)
  # resid_h2 = runif(p,0,0.6)
  resid_Va = rep(Va,p)
  resid_Ve = rep(Ve,p)
  resid_h2 = resid_Va/(resid_Va+resid_Ve)
  tot_Y_prec = 1/(resid_Va+resid_Ve)

  # G, R matrices
  G = tcrossprod(sweep(Lambda,2,sqrt(factor_h2s),'*')) + diag(resid_h2/tot_Y_prec)
  R = tcrossprod(sweep(Lambda,2,sqrt(1-factor_h2s),'*')) + diag((1-resid_h2)/tot_Y_prec)

  # fixed effect design
  n = nrow(K)
  X = cbind(1,matrix(rnorm(n*(b-1)),nrow=n))
  colnames_X = 'intercept'
  if(b > 1) colnames_X = c(colnames_X,paste0('Fixed',1:max(1,(b-1))))
  colnames(X) = colnames_X
  X_F = X[,-1]

  # cis_genotypes
  X_cis = matrix(sample(c(0,1),n*p,replace=T),n,p)
  b_cis = rnorm(p,0,sqrt(V_cis))

  # eQTL - one per factor
  X_SNP = matrix(sample(c(0,1),n*nSNP,replace=T),n,nSNP)
  B_SNP = matrix(0,nSNP,k)
  diag(B_SNP) = bSNP


  # RE design
  data = droplevels(data.frame(X,Sire=factor(pedigree$sire[children]),animal = as.factor(pedigree$id[children])))
  Z = diag(1,n)
  # Z = model.matrix(~0+Sire,data)
  # r = ncol(Z)

  B_F = matrix(rnorm((b-1)*k),nrow = (b-1),ncol=k) * sqrt(Vb)
  B = rbind(rnorm(p),matrix(0,(b-1),p)) * sqrt(Vb)

  U_F = t(K_chol) %*% matrix(rnorm(n*k,0,sqrt(factor_h2s)),n,k,byrow=T)
  E_F = matrix(rnorm(n*k,0,sqrt(1-factor_h2s)),n,k,byrow=T)
  F = X_F %*% B_F + X_SNP %*% B_SNP + Z %*% U_F + E_F

  U_R = t(K_chol) %*% matrix(rnorm(n*p,0,sqrt(resid_h2/tot_Y_prec)),n,p,byrow=T)
  E_R = matrix(rnorm(n*p,0,sqrt((1-resid_h2)/tot_Y_prec)),n,p,byrow=T)

  Y = X %*% B + F %*% t(Lambda) + U_R + E_R + sapply(1:p,function(j) X_cis[,j] * b_cis[j])

  # recover()
  Y = as.matrix(Y)
  colnames(Y) = paste0('gene',1:p)


  setup = list(
    Y = Y,
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
    X_cis = X_cis,
    b_cis = b_cis,
    X_SNP = X_SNP,
    B_SNP = B_SNP,
    # E_F = E_F,
    # E_R = E_R,
    name = name
  )
  save(setup,file=sprintf('setup_%s.RData',name))
}
