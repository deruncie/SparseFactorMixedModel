fast_BSFG_sampler_init_comb = function(priors,run_parameters){
  # require(PEIP)
  require(Matrix)
  
  # data_matrices,run_parameters,priors,current_state,Posterior,simulation = F)
  #   Preps everything for an analysis with the function MA_sampler. 
  #       Loads data from setup.mat
  #       Decides if the data is from a simulation
  #       Uses prior to find initial values for all parameters
  #       Initializes Posterior struct
  #       Saves run parameters
  #
  #   Input file setup.mat should contain:
  #     Y      gene expression data data: n x p
  #     X      fixed-effect design matrix (optional)
  #     Z_1    random effects 1 incidence matrix n x r1
  #     Z_2    random effects 2 incidence matrix n x r2 (optional)
  #     A      kinship matrix of lines: r1 x r1
  #     optional:
  #         U_act
  #         gen_factor_Lambda
  #         error_factor_Lambda
  #         h2_act
  #         G
  #         R
  #         B
  #         factor_h2s
  #         name
  
  # ----------------------- #
  # ------read data-------- #
  # ----------------------- #
  if(file.exists('../setup.RData')) {
    load('../setup.RData')
    Y = lapply(setup, function(x) x$Y)
    Mean_Y = colMeans(do.call(rbind,Y),na.rm = TRUE)
  }
  
  k = priors$k_init
  
  simulation = F
  if('gen_factor_Lambda' %in% names(setup)){
    simulation = T
    print('simulated data')
    run_parameters$setup = setup
  }
  run_parameters$simulation = simulation
  
  pops = c("LConG1","LConG2","LConG3")
  # do a loop for each Y in the list
  Y =list();U_act=list();E_act=list();B=list();Z_1=list();A=list()
  X=list();n=list();r=list();r2=list();b=list();traitnames=list()
  Z_2=list();priors$b_X_prec = list();resid_Y_prec=list()
  E_a_prec=list();E_a=list();F_h2=list();F_a=list();F=list();W_prec=list();W=list()
  Ainv=list();A_2_inv=list();invert_aI_bZAZ=list();invert_aZZt_Ainv=list();invert_aPXA_bDesignDesignT=list()
  invert_aPXA_bDesignDesignT_rand2=list()
  
  for (pop in pops){
  Y[[pop]]       = setup[[pop]]$Y
  U_act[[pop]]   = setup[[pop]]$U_act
  E_act[[pop]]   = setup[[pop]]$E_act
  B[[pop]]       = setup[[pop]]$B_act   #Is B = B_act? There is no use of B_act in this code.
  Z_1[[pop]]     = setup[[pop]]$Z_1
  A[[pop]]       = setup[[pop]]$A	
  X[[pop]]       = setup[[pop]]$X	
  n[[pop]]       = setup[[pop]]$n	
  r[[pop]]       = setup[[pop]]$r	
  traitnames[[pop]] = setup[[pop]]$traitnames
  
  
  #Determine if 'setup.mat' contains output of a simulation, based on if
  #known factor loadings are included. Affects plotting functions
  
  # Factors:
  #  initial number of factors
  
  resid_Y_prec_shape  = priors$resid_Y_prec_shape
  resid_Y_prec_rate   = priors$resid_Y_prec_rate
  E_a_prec_shape = priors$E_a_prec_shape
  E_a_prec_rate  = priors$E_a_prec_rate
  W_prec_shape = priors$W_prec_shape
  W_prec_rate  = priors$W_prec_rate
  p = ncol(Y[[1]]) 
  #k = p

  #normalize Y to have zero mean and unit variances among observed values,
  #allowing for NaNs.
  n[[pop]] = nrow(Y[[pop]])
  r[[pop]] = ncol(Z_1[[pop]])
  
  Y_missing = is.na(Y[[pop]])        # matrix recording missing data in Y
  #Mean_Y = colMeans(Y[[pop]],na.rm=T)
  VY = apply(Y[[pop]],2,var,na.rm=T)
  # Don't remove the mean and standardize the variance if it's a
  # simulation because this makes comparing to the simulated values more
  # difficult. 
  # Do we want to do this for real data, or let the user do it?
  if(simulation) {
    VY = rep(1,p)
  } else{
    Y[[pop]] = sweep(Y[[pop]],2,Mean_Y,'-')
    Y[[pop]] = sweep(Y[[pop]],2,sqrt(VY),'/')
  }
  
  
  #determine if a design matrix (X) exists (is loaded from setup.mat). If
  #not, make a dummy X-matrix with no columns.
  if(is.null(X[[pop]])) {
    X[[pop]]=matrix(0,nr = n[[pop]],nc = 0)
    B[[pop]]=matrix(0,nr = 0,nc = p)
  }
  if(nrow(X[[pop]]) == 1) X[[pop]] = t(X[[pop]])       # because I used to give X transposed
  stopifnot(nrow(X[[pop]]) == n[[pop]])
  b[[pop]] = ncol(X[[pop]])
  
  
  #Determine if a second random effects design matrix exists. If not, make a
  #dummy matrix
  # if(! 'Z_2' %in% ls()) {
  #   Z_2[[pop]]=matrix(0,nr = n,nc = 0)
  # }
    if(is.null(Z_2[[pop]])){
      Z_2[[pop]]=matrix(0,nr = n[[pop]],nc = 0)
    }
    stopifnot(nrow(Z_2[[pop]]) == n[[pop]])
    r2[[pop]] = ncol(Z_2[[pop]])
  

  
  #fixed effect priors
  
  priors$b_X_prec[[pop]]            =   1e-10*diag(1,b[[pop]])   
  
  # ----------------------------- #
  # -----Initialize variables---- #
  # ----------------------------- # 
  
  # --- transcript-level model
  # p-vector of probe residual precisions. 
  #  Prior: Gamma distribution for each element
  #       shape = resid_Y_prec_shape
  #       rate = resid_Y_prec_rate
  resid_Y_prec[[pop]]        = rgamma(p,shape = resid_Y_prec_shape,rate = resid_Y_prec_rate)
  
  
  # g-vector of specific precisions of genetic effects. 
  #  Prior: Gamma distribution for each element
  #       shape = E_a_prec_shape
  #       rate = E_a_prec_rate
  E_a_prec[[pop]]       = rgamma(p,shape = E_a_prec_shape,rate = E_a_prec_rate)
  
  # Genetic effects not accounted for by factors.
  #   Prior: Normal distribution on each element.
  #       mean = 0
  #       sd = 1./sqrt(E_a_prec)' on each row
  E_a[[pop]] = matrix(rnorm(p*r[[pop]],0,sqrt(1/E_a_prec[[pop]])),nr = r[[pop]],nc = p, byrow = T)
  
  # Latent factor heritabilties. h2 can take h2_divisions values
  #   between 0 and 1.
  #   Prior: 0.5: h2=0, .05: h2 > 0. 
  
  # change k back  to p
  F_h2[[pop]] = runif(k)
  
  # Genetic effects on the factor scores.
  #  Prior: Normal distribution for each element
  #       mean = 0
  #       sd = sqrt(F_h2') for each row.
  #F_a[[pop]] = matrix(rnorm(k*r,0,sqrt(F_h2[[pop]])),nr = r,nc = k, byrow = T)
  F_a[[pop]] = matrix(rnorm(k*r[[pop]],0,sqrt(F_h2[[pop]])),nr = r[[pop]],nc = k, byrow = T)
  
  # Full Factor scores. Combination of genetic and residual variation on
  # each factor.
  #  Prior: Normal distribution for each element
  #       mean = Z_1 * F_a
  #       sd = sqrt(1-F_h2') for each row.
  F[[pop]] = Z_1[[pop]] %*% F_a[[pop]] + matrix(rnorm(k*n[[pop]],0,sqrt(1-F_h2[[pop]])),nr = n[[pop]],nc = k, byrow = T)    
  
  # g-vector of specific precisions of genetic effects. 
  #  Prior: Gamma distribution for each element
  #       shape = E_a_prec_shape
  #       rate = E_a_prec_rate
  W_prec[[pop]]       = rgamma(p,shape = W_prec_shape,rate = W_prec_rate)
  
  # Genetic effects not accounted for by factors.
  #   Prior: Normal distribution on each element.
  #       mean = 0
  #       sd = 1./sqrt(E_a_prec)' on each row
  W[[pop]] = matrix(rnorm(p*r2[[pop]],0,sqrt(1/W_prec[[pop]])),nr = r2[[pop]], nc = p, byrow = T)
  
  # Fixed effect coefficients.
  #  Prior: Normal distribution for each element
  #       mean = 0
  #       sd = sqrt(1/fixed_effects_prec)
  B[[pop]] = matrix(rnorm(b[[pop]]*p),nr = b[[pop]], nc = p)
  
  # ------------------------------------ #
  # ----Precalculate some matrices------ #
  # ------------------------------------ #
  
  # recover()
  #invert the random effect covariance matrices
  Ainv[[pop]] = solve(A[[pop]])
  A_2_inv[[pop]] = diag(1,r2[[pop]]) #Z_2 random effects are assumed to have covariance proportional to the identity. Can be modified.
  
  #pre-calculate transformation parameters to diagonalize aI + bZAZ for fast
  #inversion: inv(aI + bZAZ) = 1/b*U*diag(1./(s+a/b))*U'
  #uses singular value decomposition of ZAZ for stability when ZAZ is low
  #rank
  #     XZ = [X_f Z_1]
  #     [U,S,~]          = svd(XZ*blkdiag(1e6*eye(b_f),A)*XZ')
  # recover()
  # r = rgamma(nrow(Z_1),1,1)
  # result = GSVD_2_c(cholcov(diag(r)),cholcov(Z_1 %*% A %*% t(Z_1)))
  # r2 = rgamma(nrow(Z_1),1,1)
  # result2 = GSVD_2_c(cholcov(diag(r2)),cholcov(Z_1 %*% A %*% t(Z_1)))
  
  result = svd(Z_1[[pop]] %*% A[[pop]] %*% t(Z_1[[pop]]))
  invert_aI_bZAZ[[pop]] = list(
    U = result$u,
    s = result$d
  )
  
  #fixed effects + random effects 1
  #diagonalize mixed model equations for fast inversion: 
  #inv(a*bdiag(priors$b_X_prec,Ainv) + b*t(cbind(X,Z_1)) %*% cbind(X,Z_1)) = U %*% diag(1/(a*s1+b*s2)) %*% t(U)
  Design= cbind(X[[pop]],Z_1[[pop]])
  Design2 = t(Design) %*% Design
  result = GSVD_2_c(cholcov(bdiag(priors$b_X_prec[[pop]],Ainv[[pop]])),cholcov(Design2))
  invert_aPXA_bDesignDesignT[[pop]] = list(
    U = t(solve(result$X)),
    s1 = diag(result$C)^2,
    s2 = diag(result$S)^2
  )
  invert_aPXA_bDesignDesignT[[pop]]$Design_U = Design %*% invert_aPXA_bDesignDesignT[[pop]]$U
  
  #random effects 2
  #diagonalize mixed model equations for fast inversion: 
  #inv(a*A_2_inv + b*Z_2'Z_2]) = U*diag(1./(a.*s1+b.*s2))*U'
  if(r2[[pop]] > 0) {
    Design = Z_2[[pop]]
    Design2 = t(Design) %*% Design
    result = GSVD_2_c(cholcov(A_2_inv[[pop]]),cholcov(Design2))
    invert_aPXA_bDesignDesignT_rand2[[pop]] = list(
      U = t(solve(result$X)),
      s1 = diag(result$C)^2,
      s2 = diag(result$S)^2
    )
    invert_aPXA_bDesignDesignT_rand2[[pop]]$Design_U = Design %*% invert_aPXA_bDesignDesignT_rand2[[pop]]$U
  } else{
    invert_aPXA_bDesignDesignT_rand2[[pop]] = list()
  }
  
  #genetic effect variances of factor traits
  # diagonalizing a*Z_1'*Z_1 + b*Ainv for fast inversion
  #diagonalize mixed model equations for fast inversion: 
  # inv(a*Z_1'*Z_1 + b*Ainv) = U*diag(1./(a.*s1+b.*s2))*U'
  #similar to fixed effects + random effects 1 above, but no fixed effects.
  ZZt = t(Z_1[[pop]]) %*% Z_1[[pop]]
  result = GSVD_2_c(cholcov(ZZt),cholcov(Ainv[[pop]]))
  invert_aZZt_Ainv[[pop]] = list(
    U = t(solve(result$X)),
    s1 = diag(result$C)^2,
    s2 = diag(result$S)^2
  )
  # the end of loop for pops
  }
  
  
  data_matrices = list(
    Y         = Y,
    Z_1       = Z_1,
    Z_2       = Z_2,
    X         = X,
    Y_missing = Y_missing
  )
  
  # ----------------------------- #
  # -----Initialize lambda   ---- #
  # ----------------------------- # 
  # Factor loading precisions (except column penalty tauh).
  #  Prior: Gamma distribution for each element. 
  #       shape = Lambda_df/2
  #       rate = Lambda_df/2
  #    Marginilizes to t-distribution with Lambda_df degrees of freedom
  #    on each factor loading, conditional on tauh
  Lambda_df   = priors$Lambda_df
  Lambda_prec = matrix(rgamma(p*k,shape = Lambda_df/2,rate = Lambda_df/2),nr = p,nc = k)
  
  # Factor penalty. tauh(h) = \prod_{i=1}^h \delta_i
  #  Prior: Gamma distribution for each element of delta
  #     delta_1:
  #       shape = delta_1_shape
  #       rate = delta_1_rate
  #     delta_2 ... delta_m:
  #       shape = delta_2_shape
  #       rate = delta_2_rate
  delta_1_shape  = priors$delta_1_shape
  delta_1_rate   = priors$delta_1_rate
  delta_2_shape  = priors$delta_2_shape
  delta_2_rate   = priors$delta_2_rate
  delta          = c(rgamma(1,shape = delta_1_shape,rate = delta_1_rate),rgamma(k-1,shape = delta_2_shape,rate = delta_2_rate))
  tauh           = cumprod(delta)
  
  # Total Factor loading precisions Lambda_prec * tauh
  Plam = sweep(Lambda_prec,2,tauh,'*')
  
  # Lambda - factor loadings
  #   Prior: Normal distribution for each element.
  #       mu = 0
  #       sd = sqrt(1/Plam)
  Lambda = matrix(rnorm(p*k,0,sqrt(1/Plam)),nr = p,nc = k)
  
  # ----------------------- #
  # -Initialize Posterior-- #
  # ----------------------- #
  Posterior = list(
    Lambda        = matrix(0,nr=0,nc=0),
    E_a           = matrix(0,nr =nrow(do.call(rbind,F_a)) ,nc = p),
    F_a           = matrix(0,nr=0,nc=0),
    F             = matrix(0,nr=0,nc=0),
    delta         = matrix(0,nr=0,nc=0),
    F_h2          = matrix(0,nr=0,nc=0),
    resid_Y_prec  = matrix(0,nr = p,nc = 0),
    E_a_prec      = matrix(0,nr = p,nc = 0),
    W_prec        = matrix(0,nr = p,nc = 0),
    B             = matrix(0,nr = b[[1]],nc = p),
    W             = matrix(0,nr = r2[[1]],nc = p)
    #E_a           = matrix(0,nr = r[[1]],nc = p)
  )
  # ----------------------- #
  # ---Save initial values- #
  # ----------------------- #
  current_state = list(
    resid_Y_prec  = resid_Y_prec,
    Lambda_prec   = Lambda_prec,
    delta         = delta,
    tauh          = tauh,
    Plam          = Plam,
    Lambda        = Lambda,
    F_h2          = F_h2,
    E_a_prec      = E_a_prec,
    W_prec        = W_prec,
    F_a           = F_a,
    F             = F,
    E_a           = E_a,
    B             = B,         
    W             = W,
    nrun 		  = 0
  )
  # ----------------------------- #
  # ----Save run parameters------ #
  # ----------------------------- #
  
  run_variables = list(
    p       = p,
    n       = n,
    r       = r,
    r2      = r2,
    b       = b,
    Mean_Y  = Mean_Y,
    VY      = VY,
    Ainv    = Ainv,
    A_2_inv = A_2_inv,
    invert_aI_bZAZ                   = invert_aI_bZAZ,
    invert_aPXA_bDesignDesignT       = invert_aPXA_bDesignDesignT,
    invert_aZZt_Ainv                 = invert_aZZt_Ainv,
    invert_aPXA_bDesignDesignT_rand2 = invert_aPXA_bDesignDesignT_rand2
  )
  
  RNG = list(
    Random.seed = .Random.seed,
    RNGkind = RNGkind()
  )
  
  return(list(
    data_matrices  = data_matrices,
    run_parameters = run_parameters,
    run_variables  = run_variables,
    priors         = priors,
    current_state  = current_state,
    Posterior      = Posterior,
    simulation     = simulation,
    RNG            = RNG,
    traitnames     = traitnames
  ))
}
