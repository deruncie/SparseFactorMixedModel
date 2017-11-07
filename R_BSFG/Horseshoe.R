a0 = 1
b0 = 1

get_truncated_exp = function(mn,trunc_point,v){
  p = length(mn)
  eps = .Machine$double.eps
  r = mn * trunc_point
  sml = abs(r) < .Machine$double.eps
  x = rep(0,p)
  tmp = rep(0,p)
  if(any(sml)){
    print(paste('r',sum(sml)))
    tmp[sml] = expm1(-r[sml])
  }
  tmp[!sml] = (exp(-r[!sml])-1)
  tmp = tmp * v

  sml = abs(tmp) < eps
  if(any(sml)){
    print(paste('tmp',sum(sml)))
    x[sml] = -log1p(temp[sml])
  }
  x[!sml] = -log(1+tmp[!sml])
  x = x/m
  return(x)
}

horseshoe2 = function(y,X,a0,b0,nIter,burn,verbose=F){
  #asdf
  n = length(y)
  p = ncol(X)

  beta = rep(0,p)


  s2 = 1/rgamma(1,shape=a0,rate=b0)

  xi = 1/rgamma(1,shape=1/2,rate=1)
  eta = 1/rgamma(p,shape=1/2,rate=1)

  # xi = 1
  # eta = rep(1,p)

  tau2 = 1
  lambda2 = rep(1,p)

  chol_R = as(diag(1,n),'dgCMatrix')
  prior_mean = rep(0,p)

  XtX = t(X) %*% X

  samples = matrix(0,nIter,p+p+1+1)
  for(i in (-burn):nIter){
    if(verbose && i %% (nIter/10) == 0) print(i)
    xi_new = exp(rnorm(1,log(xi),.8))

    if(p > n) {
      XDXt = sweep(X,2,eta,'/') %*% t(X)
      M_new = diag(1,n) + XDXt/xi_new
      chol_M_new = chol(M_new)
      log_det_M_new = 2*sum(log(diag(chol_M_new)))
      scaled_resid2_new = solve(t(chol_M_new),y)^2

      M = diag(1,n) + XDXt/xi
      chol_M = chol(M)
      log_det_M = 2*sum(log(diag(chol_M)))
      scaled_resid2 = solve(t(chol_M),y)^2
    } else{
      Xty = t(X) %*% y
      inner = XtX
      diag(inner) = diag(inner) + c(eta)*xi_new
      Miy = y - X %*% solve(inner,Xty)
      scaled_resid2_new =  y*Miy
      log_det_M_new = -p*log(xi_new) + sum(-log(eta)) + determinant(inner,log=T)$modulus

      inner = XtX
      diag(inner) = diag(inner) + c(eta)*xi
      Miy = y - X %*% solve(inner,Xty)
      scaled_resid2 =  y*Miy
      log_det_M = -p*log(xi) + sum(-log(eta)) + determinant(inner,log=T)$modulus
    }

    log_L_new = -1/2*log_det_M_new - (n+a0)/2*log(b0/2 + sum(scaled_resid2_new/2))
    log_prior_new = -log(sqrt(xi_new)*(1+xi_new))

    log_L = -1/2*log_det_M - (n+a0)/2*log(b0/2 + sum(scaled_resid2/2))
    log_prior = -log(sqrt(xi)*(1+xi))

    log_MH_ratio = log_L_new + log_prior_new + log(xi_new) - (log_L + log_prior + log(xi))

    if(log(runif(1)) < log_MH_ratio){
      xi = xi_new
      scaled_resid2 = scaled_resid2_new
    }

    s2 = 1/rgamma(1,shape = (a0+n)/2,rate = (sum(scaled_resid2) + b0)/2)

    prior_prec = eta*xi/(s2)
    randn_theta = rnorm(p)
    randn_e = matrix(0,0,0)
    if(p>n) randn_e = rnorm(n)
    beta = sample_MME_single_diagK(y,X,prior_mean,prior_prec,chol_R,1/s2,randn_theta,randn_e)


    u = runif(p,0,1/(eta+1))
    r = (1-u)/u
    v = runif(p,0,1)
    m = beta^2*xi/(2*s2)
    eta = get_truncated_exp(m,r,v)
    # eta = -log(1-(1-exp(-m*r))*v)/m

    tau2 = 1/xi
    lambda2 = 1/eta

    if(i>0) samples[i,] = c(beta,lambda2,tau2,s2)

  }
  return(samples)
}

# This doesn't shrink as well as horseshoe, but easier and faster and works pretty well. Maybe horeshoe better for genotypes
ard = function(y,X,a0,b0,a1,b1,a2,b2,nIter,burn,verbose=F){
  n = length(y)
  p = ncol(X)

  beta = rep(0,p)

  s2 = 1/rgamma(1,shape=a0,rate=b0)

  lambda = rep(1,p)
  tau = 1

  chol_R = as(diag(1,n),'dgCMatrix')
  prior_mean = rep(0,p)


  samples = matrix(0,nIter,p+p+1+1)
  for(i in (-burn):nIter){
    if(verbose && i %% (nIter/10) == 0) print(i)

    prior_prec = lambda*tau
    randn_theta = rnorm(p)
    randn_e = matrix(0,0,0)
    if(p>n) randn_e = rnorm(n)
    beta = sample_MME_single_diagK(y,X,prior_mean,prior_prec,chol_R,1/s2,randn_theta,randn_e)

    resids = y - X %*% beta
    s2 = 1/rgamma(1,shape = a0+(n)/2,rate = b0+sum(resids^2)/2)

    lambda = rgamma(p,shape = a1+1/2,rate = b1+beta^2*tau)
    tau = rgamma(1,shape = a2+p/2,rate = b2+beta^2*lambda)

    if(i>0) samples[i,] = c(beta,lambda,tau,s2)

  }
  return(samples)
}







# This doesn't mix as well as Johndrow's horeshoe1 above, but approximatly as well if sample hyperparameters multipel times
horseshoe1 = function(y,X,a0,b0,nIter,burn,thin=1,verbose = F){
  n = length(y)
  p = ncol(X)

  beta = rep(0,p)

  s2 = 1/rgamma(1,shape=a0,rate=b0)
  xi = 1/rgamma(1,shape=1/2,rate=1)
  nu = 1/rgamma(p,shape=1/2,rate=1)

  xi =1
  nu = rep(1,p)

  tau2 = 1
  lambda2 = rep(1,p)

  chol_R = as(diag(1,n),'dgCMatrix')
  prior_mean = rep(0,p)

  samples = matrix(0,nIter/thin,p+p+1+1)#+p+1
  for(i in (-burn):nIter){
    if(verbose && i %% (nIter/10) == 0) print(i)

    prior_prec = 1/(lambda2*tau2*s2)
    randn_theta = rnorm(p)
    randn_e = matrix(0,0,0)
    if(p>n) randn_e = rnorm(n)
    beta = sample_MME_single_diagK(y,X,prior_mean,prior_prec,chol_R,1/s2,randn_theta,randn_e)

    resids = y - X %*% beta
    s2 = 1/rgamma(1,shape = (n+p)/2,rate = sum(resids^2)/2 + sum(beta^2/(lambda2*tau2))/2)
for(j in 1:10){
    lambda2 = 1/rgamma(p,shape = 1,rate = 1/nu + beta^2/(2*tau2*s2))
    tau2 = 1/rgamma(1,shape=(p+1)/2,rate = 1/xi + 1/(2*s2)*sum(beta^2/lambda2))

    nu = 1/rgamma(p,shape=1,rate=1 + 1/lambda2)

    xi = 1/rgamma(1,shape=1,rate=1+1/tau2)
}

    if(i>0 && i %% thin == 0) samples[ceiling(i/thin),] = c(beta,lambda2,tau2,s2)#,nu,xi

  }
  return(samples)
}
