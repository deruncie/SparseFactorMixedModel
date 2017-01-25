
fast_BSFG_ipx_sampler_init = function(priors,run_parameters){
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
        names(setup) = sub('.','_',names(setup),fixed=T)
    }
    else{
        require(R.matlab)
        setup = readMat('../setup.mat')
        for(i in 1:10) names(setup) = sub('.','_',names(setup),fixed=T)
    }
	Y       = setup$Y
	U_act   = setup$U_act
	E_act   = setup$E_act
	B       = setup$B_act   #Is B = B_act? There is no use of B_act in this code.
	Z_1     = Matrix(setup$Z_1,sparse=T)
	A       = Matrix(setup$A,sparse=T)
	X       = setup$X	
	n       = setup$n	
	r       = setup$r	
	traitnames = setup$traitnames

    #Determine if 'setup.mat' contains output of a simulation, based on if
    #known factor loadings are included. Affects plotting functions
    simulation = F
    if('gen_factor_Lambda' %in% names(setup)){
    	simulation = T
    	print('simulated data')
    	run_parameters$setup = setup
        run_parameters$name = setup$name
    }
    run_parameters$simulation = simulation

    #normalize Y to have zero mean and unit variances among observed values,
    #allowing for NaNs.
    n = nrow(Y)
    p = ncol(Y)
    r = ncol(Z_1)

    Y_missing = is.na(Y)        # matrix recording missing data in Y
    Mean_Y = colMeans(Y,na.rm=T)
    VY = apply(Y,2,var,na.rm=T)
    # Don't remove the mean and standardize the variance if it's a
    # simulation because this makes comparing to the simulated values more
    # difficult. 
    # Do we want to do this for real data, or let the user do it?
    if(simulation) {
        VY = rep(1,p)
    } else{
    	Y = sweep(Y,2,Mean_Y,'-')
    	Y = sweep(Y,2,sqrt(VY),'/')
    }
       
    
    #determine if a design matrix (X) exists (is loaded from setup.mat). If
    #not, make a dummy X-matrix with no columns.
    if(is.null(X)) {
        X=matrix(0,nr = n,nc = 0)
        B=matrix(0,nr = 0,nc = p)
    }
    if(nrow(X) == 1) X = t(X)       # because I used to give X transposed
    stopifnot(nrow(X) == n)
    b = ncol(X)
    
    Z_1_sparse = Matrix(Z_1,sparse=T)


    data_matrices = list(
		    Y           = Y,
            Z_1         = Z_1,
            Z_1_sparse  = Z_1_sparse,
		    X           = X,
            Y_missing   = Y_missing
    		)
    
    #fixed effect priors
    
    priors$b_X_prec            =   1e-10*diag(1,b)
        
# ----------------------------- #
# -----Initialize variables---- #
# ----------------------------- # 

   # --- transcript-level model
    # p-vector of gene precisions after taking removing effect of factors. 
	#  Prior: Gamma distribution for each element
    #       shape = tot_Y_prec_shape
    #       rate = tot_Y_prec_rate
    tot_Y_prec = with(priors,rgamma(p,shape = tot_Y_prec_shape,rate = tot_Y_prec_rate))
   
    # Factors:
    #  initial number of factors
    k = priors$k_init

    # Factor loading precisions (except column penalty tauh).
	 #  Prior: Gamma distribution for each element. 
     #       shape = Lambda_df/2
     #       rate = Lambda_df/2
     #    Marginilizes to t-distribution with Lambda_df degrees of freedom
     #    on each factor loading, conditional on tauh
    Lambda_prec = with(priors,matrix(rgamma(p*k,shape = Lambda_df/2,rate = Lambda_df/2),nr = p,nc = k))
    
    # Factor penalty. tauh(h) = \prod_{i=1}^h \delta_i
	 #  Prior: Gamma distribution for each element of delta
     #     delta_1:
     #       shape = delta_1_shape
     #       rate = delta_1_rate
     #     delta_2 ... delta_m:
     #       shape = delta_2_shape
     #       rate = delta_2_rate
    delta          = with(priors,c(rgamma(1,shape = delta_1_shape,rate = delta_1_rate),rgamma(k-1,shape = delta_2_shape,rate = delta_2_rate)))
    tauh           = cumprod(delta)
    
    # Total Factor loading precisions Lambda_prec * tauh
    Plam = sweep(Lambda_prec,2,tauh,'*')
    
    # Lambda - factor loadings
     #   Prior: Normal distribution for each element.
     #       mu = 0
     #       sd = sqrt(1/Plam)
    Lambda = matrix(rnorm(p*k,0,sqrt(1/Plam)),nr = p,nc = k)

    # p-vector of heritabilities for residuals on factors. 
	#  Prior: discrete prior on interval [0,1)
    h2s = with(run_parameters,seq(0,1,length = h2_divisions+1)[1:h2_divisions])
    resid_h2 = with(priors,sample(h2s,p,prob = h2_priors_resids, replace = T))
    
    # Genetic effects not accounted for by factors.
    #   Prior: Normal distribution on each element.
    #       mean = 0
    #       sd = 1./sqrt(E_a_prec)' on each row
    E_a = matrix(rnorm(p*r,0,sqrt(resid_h2 * tot_Y_prec)),nr = r,nc = p, byrow = T)
        
    # Latent factor variances
    tot_F_prec = with(priors,rgamma(k,shape = tot_F_prec_shape,rate = tot_F_prec_rate))
    F_h2 = with(priors,sample(h2s,k,prob = h2_priors_factors, replace = T))
        
    # Genetic effects on the factor scores.
	#  Prior: Normal distribution for each element
    #       mean = 0
    #       sd = sqrt(F_h2') for each row.
    F_a = matrix(rnorm(k*r,0,sqrt(F_h2)),nr = r,nc = k, byrow = T)
    
    # Full Factor scores. Combination of genetic and residual variation on
    # each factor.
	#  Prior: Normal distribution for each element
    #       mean = Z_1 * F_a
    #       sd = sqrt(1-F_h2') for each row.
    F = as.matrix(Z_1 %*% F_a + matrix(rnorm(k*n,0,sqrt(1-F_h2)),nr = n,nc = k, byrow = T))
    
    # Fixed effect coefficients.
	#  Prior: Normal distribution for each element
    #       mean = 0
    #       sd = sqrt(1/fixed_effects_prec)
    B = matrix(rnorm(b*p),nr = b, nc = p)

    
# ----------------------- #
# -Initialize Posterior-- #
# ----------------------- #
    Posterior = list(
		    Lambda        = matrix(0,nr=0,nc=0),
		    F_a           = matrix(0,nr=0,nc=0),
		    F             = matrix(0,nr=0,nc=0),
		    delta         = matrix(0,nr=0,nc=0),
            tot_F_prec    = matrix(0,nr=0,nc=0),
            F_h2          = matrix(0,nr=0,nc=0),
            tot_Y_prec    = matrix(0,nr = p,nc = 0),
            resid_h2      = matrix(0,nr = p,nc = 0),
            resid_Y_prec  = matrix(0,nr = p,nc = 0),
            E_a_prec      = matrix(0,nr = p,nc = 0),
		    B             = matrix(0,nr = b,nc = p),
		    E_a           = matrix(0,nr = r,nc = p)
    	)
# ----------------------- #
# ---Save initial values- #
# ----------------------- #
    current_state = list(
    		Lambda_prec   = Lambda_prec,
    		delta         = delta,
    		tauh          = tauh,
    		Plam          = Plam,
    		Lambda        = Lambda,
            tot_Y_prec    = tot_Y_prec,
            resid_h2      = resid_h2,
            tot_F_prec    = tot_F_prec,
            F_h2          = F_h2,
    		F_a           = F_a,
    		F             = F,
    		E_a           = E_a,
    		B             = B,
    		nrun 		  = 0
    	)


# ------------------------------------ #
# ----Precalculate some matrices------ #
# ------------------------------------ #

    # recover()
    #invert the random effect covariance matrices
    Ainv = solve(A)
    chol_Ainv = chol(Ainv)

    #pre-calculate transformation parameters to diagonalize aI + bZAZ for fast
    #inversion: inv(aI + bZAZ) = 1/b*U*diag(1./(s+a/b))*U'
    #uses singular value decomposition of ZAZ for stability when ZAZ is low
    #rank
#     XZ = [X_f Z_1]
#     [U,S,~]          = svd(XZ*blkdiag(1e6*eye(b_f),A)*XZ')
    
    result = svd(Z_1 %*% A %*% t(Z_1))
    invert_aI_bZAZ = list(
        U = result$u,
        s = result$d
    )
    invert_aI_bZAZ$U_sparse = Matrix(invert_aI_bZAZ$U,sparse=T)
   
    #genetic effect variances of factor traits
    # diagonalizing a*Z_1'*Z_1 + b*Ainv for fast inversion
    #diagonalize mixed model equations for fast inversion: 
    # inv(a*Z_1'*Z_1 + b*Ainv) = U*diag(1./(a.*s1+b.*s2))*U'
    #similar to fixed effects + random effects 1 above, but no fixed effects.
    ZZt = t(Z_1) %*% Z_1
    result = GSVD_2_c(cholcov(ZZt),cholcov(Ainv))
	invert_aZZt_Ainv = list(
		U = t(solve(result$X)),
			s1 = diag(result$C)^2,
			s2 = diag(result$S)^2
		)

# ----------------------------- #
# ----Save run parameters------ #
# ----------------------------- #

	run_variables = list(
			p                = p,
			n                = n,
			r                = r,
			b                = b,
			Mean_Y           = Mean_Y,
			VY               = VY,
            Ainv             = Ainv,
            chol_Ainv        = chol_Ainv,
			invert_aI_bZAZ   = invert_aI_bZAZ,
			invert_aZZt_Ainv = invert_aZZt_Ainv
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
