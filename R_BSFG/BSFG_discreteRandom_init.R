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

BSFG_discreteRandom_init = function(Y, fixed, random, data, priors, run_parameters, A_mats = NULL, A_inv_mats = NULL, 
									fixed_Factors = NULL, scaleY = TRUE,
									ncores = detectCores(),simulation = F,setup = NULL,verbose=T){
	require(Matrix)
	require(parallel)

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

	# build Z matrices from random model
	RE_names = rownames(attr(terms(random),'factors'))
	n_RE = length(RE_names)
	Z_matrices = lapply(RE_names,function(re) {
		Z = Matrix(model.matrix(formula(sprintf('~0 + %s',re)),data),sparse = TRUE)
		Z[,paste0(re,levels(data[[re]]))]
	})
	names(Z_matrices) = RE_names
	r_RE = sapply(Z_matrices,function(x) ncol(x))

	Z_all = do.call(cbind,Z_matrices)

	h2_divisions = expand.grid(lapply(RE_names,function(re) 0:run_parameters$discrete_divisions)) / run_parameters$discrete_divisions
	colnames(h2_divisions) = RE_names
	h2_divisions = t(h2_divisions[rowSums(h2_divisions) < 1,,drop=FALSE])

	data_matrices = list(
		Y             = Y,
		Y_missing     = Y_missing,
		X             = X,
		Z_matrices    = Z_matrices,
		Z_all         = Z_all,
		h2_divisions  = h2_divisions
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

# ----------------------------- #
# -----Initialize variables---- #
# ----------------------------- # 

  # Factors loadings:
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

  # Factor scores:
     # p-vector of factor precisions. Note - this is a 'redundant' parameter designed to give the Gibbs sampler more flexibility 
	 #  Prior: Gamma distribution for each element
     #       shape = tot_F_prec_shape
     #       rate = tot_F_prec_rate
    tot_F_prec = with(priors,rgamma(k,shape = tot_F_prec_shape,rate = tot_F_prec_rate))

    # Factor discrete variances
     # k-matrix of n_RE x k with 
    F_h2_index = sample(1:ncol(h2_divisions),k,replace=T)
    F_h2 = h2_divisions[,F_h2_index,drop=FALSE]

    F_a = lapply(RE_names,function(effect){
    	matrix(rnorm(r_RE[effect] * k, 0, sqrt(F_h2[effect,] / tot_F_prec)),ncol = k, byrow = T)
    })
    names(F_a) = RE_names

    F = matrix(rnorm(n * k, 0, sqrt(rowSums(t(F_h2) / tot_F_prec))),ncol = k, byrow = T)
    for(effect in RE_names) {
    	F = F + Z_matrices[[effect]] %*% F_a[[effect]]
    }
    F = as.matrix(F)
    F_a = do.call(rbind,F_a)

  # residuals
     # p-vector of factor precisions. Note - this is a 'redundant' parameter designed to give the Gibbs sampler more flexibility 
	 #  Prior: Gamma distribution for each element
     #       shape = tot_Y_prec_shape
     #       rate = tot_Y_prec_rate
    tot_Y_prec = with(priors,rgamma(p,shape = tot_Y_prec_shape,rate = tot_Y_prec_rate))

    # Resid discrete variances
     # p-matrix of n_RE x p with 
    resid_h2_index = sample(1:ncol(h2_divisions),p,replace=T)
    resid_h2 = h2_divisions[,resid_h2_index,drop=FALSE]

    E_a = do.call(rbind,lapply(RE_names,function(effect){
    	matrix(rnorm(r_RE[effect] * p, 0, sqrt(resid_h2[effect,] / tot_Y_prec)),ncol = p, byrow = T)
    }))

  # Fixed effects
    B = matrix(rnorm(b*p), ncol = p)

# ----------------------- #
# ---Save initial values- #
# ----------------------- #
    current_state = list(
		Lambda_prec      = Lambda_prec,
		delta            = delta,
		tauh             = tauh,
		Plam             = Plam,
		Lambda           = Lambda,
		tot_F_prec       = tot_F_prec,
		F_h2_index       = F_h2_index,
		F_h2             = F_h2,
		F_a              = F_a,
		F                = F,
		tot_Y_prec       = tot_Y_prec,
		resid_h2_index   = resid_h2_index,
		resid_h2         = resid_h2,
		E_a              = E_a,
		B                = B,
		nrun             = 0
    	)

    
# ----------------------- #
# -Initialize Posterior-- #
# ----------------------- #

    Posterior = clear_Posterior(list(current_state = current_state))$Posterior

# ------------------------------------ #
# ----Precalculate ZAZts, chol_As ---- #
# ------------------------------------ #
    if(verbose) print('Pre-calculating matrices')

	fix_A = function(x) forceSymmetric(drop0(x,tol = 1e-10))

	A_mats = lapply(RE_names,function(re) {
		if(re %in% names(A_mats)) {
			A = A_mats[[re]]
		} else if(re %in% names(A_inv_mats)){
 			A = solve(A_inv_mats[[re]])
 			rownames(A) = rownames(A_inv_mats[[re]])
		} else{
			A = Diagonal(ncol(Z_matrices[[re]]))
			rownames(A) = levels(data[[re]])
		}
		index = match(sub(re,'',colnames(Z_matrices[[re]])),rownames(A))  # A must have rownames
		fix_A(A[index,index])
	})
	names(A_mats) = RE_names

	Ai_mats = lapply(RE_names,function(re) {
		if(re %in% names(A_inv_mats)){
			Ai = A_inv_mats[[re]]
			index = match(sub(re,'',colnames(Z_matrices[[re]])),rownames(Ai))  # Ai must have rownames
			Ai = Ai[index,index]
		} else {
			Ai = solve(A_mats[[re]])
		}
		Ai = fix_A(Ai)
		Ai
	})
	chol_Ai_mats = lapply(Ai_mats,chol)

	#C
	ZtZ = crossprod(Z_all)
	make_Chol_Ai = function(chol_Ai_mats,h2s){
		do.call(bdiag,lapply(1:length(h2s),function(i) {
			if(h2s[i] == 0) return(Diagonal(nrow(chol_Ai_mats[[i]]),Inf))  # if h2==0, then we want a Diagonal matrix with Inf diagonal.
			chol_Ai = chol_Ai_mats[[i]]
			chol_Ai@x = chol_Ai@x / sqrt(h2s[i])
			chol_Ai
		}))
	}
	# setup of symbolic Cholesky of C
	Ai = forceSymmetric(crossprod(make_Chol_Ai(chol_Ai_mats,rep(1,n_RE)/(n_RE+1))))
	Cholesky_C = Cholesky(ZtZ + Ai)

	randomEffect_C_Choleskys = mclapply(1:ncol(h2_divisions),function(i) {    	
		if(i %% 100 == 0 && verbose) print(sprintf('randomEffects_C %d of %d',i,ncol(h2_divisions)))
		h2s = h2_divisions[,i]
		chol_A_inv = make_Chol_Ai(chol_Ai_mats,h2s)
		Ai = crossprod(chol_A_inv)
		C = ZtZ/(1-sum(h2s))
		C = C + Ai
		Cholesky_C_i = update(Cholesky_C,forceSymmetric(C))

		return(list(Cholesky_C = Cholesky_C_i, chol_A_inv = chol_A_inv))
	},mc.cores = ncores)

	# Sigma
	ZAZts = list()
	for(i in 1:n_RE){
		ZAZts[[i]] = forceSymmetric(Z_matrices[[i]] %*% A_mats[[i]] %*% t(Z_matrices[[i]]))
	}


	make_Sigma = function(ZAZts,h2s){
		R = 0
		for(i in 1:length(h2s)){
			R = R + h2s[i]*ZAZts[[i]]
		}
		forceSymmetric(R + (1-sum(h2s)) * Diagonal(nrow(R)))
	}

	# setup of symbolic Cholesky of Sigma
	Sigma = make_Sigma(ZAZts,h2_divisions[,2])
	Cholesky_Sigma_base = Cholesky(Sigma,perm=T,super=T)
	stopifnot(!isLDL(Cholesky_Sigma_base))
	Sigma_Perm = expand(Cholesky_Sigma_base)$P
	if(all(diag(Sigma_Perm))) Sigma_Perm = NULL

	Sigma_Choleskys = mclapply(1:ncol(h2_divisions),function(i) {
		if(i %% 100 == 0 && verbose) print(sprintf('Sigma_Choleskys %d of %d',i,ncol(h2_divisions)))
		Sigma = forceSymmetric(make_Sigma(ZAZts,h2_divisions[,i]))
		stopifnot(class(Sigma) == 'dsCMatrix')
		Cholesky_Sigma = update(Cholesky_Sigma_base,Sigma)
		log_det = 2*determinant(Cholesky_Sigma,logarithm=T)$modulus
		if(is.null(Sigma_Perm)) {
			chol_Sigma = expand(Cholesky_Sigma)$L
		} else{
			chol_Sigma = t(Sigma_Perm) %*% expand(Cholesky_Sigma)$L
		}
		list(log_det = log_det,Cholesky_Sigma = Cholesky_Sigma,chol_Sigma=chol_Sigma,Sigma = Sigma)
	},mc.cores = ncores)

# ----------------------------- #
# ----Save run parameters------ #
# ----------------------------- #

	run_variables = list(
			p                        = p,
			n                        = n,
			r_RE                     = r_RE,
			b                        = b,
			Mean_Y                   = Mean_Y,
			VY                       = VY,
			Sigma_Choleskys          = Sigma_Choleskys,
			Sigma_Perm               = Sigma_Perm,
			randomEffect_C_Choleskys = randomEffect_C_Choleskys,
			chol_Ai_mats             = chol_Ai_mats
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
