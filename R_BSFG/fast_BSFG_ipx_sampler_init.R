load_simulation_data = function(){
    require(Matrix)
    if(file.exists('../setup.RData')) {
        load('../setup.RData')
        names(setup) = sub('.','_',names(setup),fixed=T)
    }
    else{
        require(R.matlab)
        setup = readMat('../setup.mat')
        for(i in 1:10) names(setup) = sub('.','_',names(setup),fixed=T)
    }
    r = dim(setup$A)[1]
    n = nrow(setup$Y)
    data = data.frame(Group = gl(r,n/r))
    rownames(setup$A) = data$Group
    return(list(Y = setup$Y, data = data, A_mats = list(Group = setup$A),setup = setup))
}



fast_BSFG_ipx_sampler_init = function(Y, fixed, random, data, priors, run_parameters, A_mat = NULL, A_inv_mat = NULL, 
                                    fixed_Factors = NULL, scaleY = TRUE,
                                    simulation = F,setup = NULL,verbose=T){
# function(Y, fixed, random, data, priors, run_parameters, scaleY = TRUE,simulation = FALSE,setup = NULL){
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
#     Z      random effects 1 incidence matrix n x r
#     A      kinship matrix of lines: r x r
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

    run_parameters$setup = setup
    run_parameters$name = setup$name
    run_parameters$simulation = simulation

    # model dimensions
    n = nrow(Y)
    p = ncol(Y)
    traitnames = colnames(Y)

    # missing data
    Y_missing = is.na(Y)

    # scale Y
    if(scaleY){
        Mean_Y = colMeans(Y,na.rm=T)
        VY = apply(Y,2,var,na.rm=T)
        Y = sweep(Y,2,Mean_Y,'-')
        Y = sweep(Y,2,sqrt(VY),'/')
    } else {
        Mean_Y = rep(0,p)
        VY = rep(1,p)
    }

    # build X from fixed model
    X = model.matrix(fixed,data)
    b = ncol(X)

    # build Z for random effect
    RE_name = rownames(attr(terms(random),'factors'))[1]
    Z = Matrix(model.matrix(formula(sprintf('~0 + %s',RE_name)),data),sparse = TRUE)
    Z = Z[,paste0(RE_name,levels(data[[RE_name]]))]    
    r = ncol(Z)

    if(is.null(A_mat)) {
        if(is.null(A_inv_mat)){
            A = Diagonal(ncol(Z))
        } else{
            A = Matrix(forceSymmetric(solve(A_inv_mat)))
        }
    } else{
        A = forceSymmetric(A_mat)
    }

    data_matrices = list(
            Y         = Y,
            Z         = Z,
            X         = X,
            Y_missing = Y_missing,
            A         = A
    		)        

# ----------------------------- #
# ----- re-formulate priors --- #
# ----------------------------- # 
    priors$tot_Y_prec_shape = with(priors$tot_Y_var,V * nu)
    priors$tot_Y_prec_rate  = with(priors$tot_Y_var,nu - 2)
    priors$tot_F_prec_shape = with(priors$tot_F_var,V * nu)
    priors$tot_F_prec_rate  = with(priors$tot_F_var,nu - 2)
    priors$delta_1_shape    = with(priors$delta_1,V * nu)
    priors$delta_1_rate     = with(priors$delta_1,nu - 2)
    priors$delta_2_shape    = with(priors$delta_2,V * nu)
    priors$delta_2_rate     = with(priors$delta_2,nu - 2)
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
    k = run_parameters$k_init

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
    #       mean = Z   * F_a
    #       sd = sqrt(1-F_h2') for each row.
    F = as.matrix(Z   %*% F_a + matrix(rnorm(k*n,0,sqrt(1-F_h2)),nr = n,nc = k, byrow = T))
    
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
#     XZ = [X_f Z  ]
#     [U,S,~]          = svd(XZ*blkdiag(1e6*eye(b_f),A)*XZ')
    
    result = svd(Z %*% A %*% t(Z))
    invert_aI_bZAZ = list(
        U = Matrix(result$u),
        s = result$d
    )
   
    #genetic effect variances of factor traits
    # diagonalizing a*Z  '*Z   + b*Ainv for fast inversion
    #diagonalize mixed model equations for fast inversion: 
    # inv(a*Z  '*Z   + b*Ainv) = U*diag(1./(a.*s1+b.*s2))*U'
    #similar to fixed effects + random effects 1 above, but no fixed effects.
    ZZt = crossprod(Z)

    # Z1 = matrix(rnorm(10*4),10)
    # A1 = matrix(rnorm(16),4);A1 = crossprod(A1) + diag(1,4)
    # ZZt1 = crossprod(Z1)
    # A1inv = solve(A1)
    # result = GSVD_2_c(cholcov(ZZt1),cholcov(A1inv))
    # r = svd(t(solve(t(cholcov(A1inv)),t(cholcov(ZZt1)))))
    # norm_factor = sqrt(1 + r$d^2)
    # c = r$d / norm_factor
    # s = 1/norm_factor
    # X = sweep(t(cholcov(A1inv)) %*% r$v,2,norm_factor,'*')

    # a=3;b=5;
    # solve(a * ZZt1 + b * A1inv)
    # t(solve(result$X)) %*% diag(1/(a*c^2 + b*s^2)) %*% solve(result$X)

    # chol_base = Cholesky(forceSymmetric(Matrix(ZZt1 + A1inv,sparse=T)),super=T,perm=F)

    # r = svd(solve(A1inv,ZZt1))
    # r2 = svd(ZZt1 %*% A1)

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
