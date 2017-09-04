new_halfSib_simulation_2env = function(name, nSire,nRep,p, b, factor_h2s, Va = 0.2, Ve = 0.2,Vb = 0, numeff = NULL){
  library(Matrix)
  Sire = gl(nSire,nRep)
  K = .25*tcrossprod(Matrix(model.matrix(~0+Sire))) + diag(.75,length(Sire))
  rownames(K) = 1:nrow(K)

    K_chol = chol(K+diag(1e-6,nrow(K)))

  # Lambda matrix
    # idea here: traits 1:p in env 1, p+1:2p in env2. So no phenotype correlation, but still genetic correlation
    # will have within-environment covariance identical



  # factor_h2s = rep(0,k)
  # factor_h2s[2+1:k_G] = runif(k_G)
  k = length(factor_h2s)
  Lambda = matrix(0,p,k)
  if(is.null(numeff)) numeff = sample((p/30):(p/4),k,replace=T)
  for(h in 1:k){
    Lambda[sample(1:p,numeff[h]),h] = rnorm(numeff[h])
  }
  Lambda = Lambda[,order(-diag(t(Lambda) %*% Lambda))]
  cols=1:k
  g_cols = factor_h2s>0
  if(sum(g_cols) == 0) g_cols = 1:k
  Lambda = Lambda[do.call("order", unname(split(-abs(Lambda[,cols[g_cols]]), col(Lambda[,cols[g_cols]])))),]

  Lambda = rbind(Lambda,Lambda)
  p = 2*p

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

  # RE design
  data = droplevels(data.frame(Rep = 1:nRep, Sire=Sire,animal = rownames(K)))
  data$Env = c(1,2)[(data$Rep > nRep/2)+1]
  Z = diag(1,n)
  # Z = model.matrix(~0+Sire,data)
  # r = ncol(Z)

  B_F = matrix(rnorm((b-1)*k),nrow = (b-1),ncol=k) * sqrt(Vb)
  B = rbind(rnorm(p),matrix(0,(b-1),p)) * sqrt(Vb)

  U_F = t(K_chol) %*% matrix(rnorm(n*k,0,sqrt(factor_h2s)),n,k,byrow=T)
  E_F = matrix(rnorm(n*k,0,sqrt(1-factor_h2s)),n,k,byrow=T)
  F = X_F %*% B_F + Z %*% U_F + E_F

  U_R = t(K_chol) %*% matrix(rnorm(n*p,0,sqrt(resid_h2/tot_Y_prec)),n,p,byrow=T)
  E_R = matrix(rnorm(n*p,0,sqrt((1-resid_h2)/tot_Y_prec)),n,p,byrow=T)
  Y = X %*% B + F %*% t(Lambda) + U_R + E_R

  # recover()
  Y = as.matrix(Y)
  colnames(Y) = paste0('gene',1:p)
  Y[data$Env == 1, p/2+1:(p/2)] = NA
  Y[data$Env == 2, 1:(p/2)] = NA

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
    F = F,
    # E_F = E_F,
    # E_R = E_R,
    name = name
  )
  save(setup,file='setup.RData')
}
