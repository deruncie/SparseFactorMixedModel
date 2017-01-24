BSFG_discreteRandom_loadData = function(){
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
	return(list(Y = setup$Y, data = data, randomEffects = list(Group = setup$A),setup = setup))
}

BSFG_discreteRandom_init = function(Y, fixed, randomEffects, data, priors, run_parameters, scaleY = TRUE,simulation = F,setup = NULL){
	require(Matrix)

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
	RE_names = names(randomEffects)
	n_RE = length(randomEffects)
	r_RE = sapply(RE_names,function(x) dim(randomEffects[[x]])[1])
	Z_matrices = lapply(RE_names,function(re) {
		Matrix(model.matrix(formula(sprintf('~0 + %s',re)),data),sparse = TRUE)
	})
	names(Z_matrices) = RE_names
	Z_all = do.call(cbind,Z_matrices)

	# check Z_matrix dimensions 
	for(re in RE_names){
		stopifnot(dim(Z_matrices[[re]])[2] == r_RE[[re]])
	}

	h2_divisions = expand.grid(lapply(RE_names,function(re) 0:run_parameters$discrete_divisions)) / run_parameters$discrete_divisions
	colnames(h2_divisions) = RE_names
	h2_divisions = t(h2_divisions[rowSums(h2_divisions) < 1,,drop=FALSE])

	data_matrices = list(
		Y             = Y,
		Y_missing     = Y_missing,
		X             = X,
		Z_matrices    = Z_matrices,
		Z_all         = Z_all,
		randomEffects = randomEffects,
		h2_divisions  = h2_divisions
		)


# ----------------------------- #
# -----Initialize variables---- #
# ----------------------------- # 

  # Factors loadings:
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
    	matrix(rnorm(r_RE[effect] * k, 0, sqrt(resid_h2[effect,] / tot_Y_prec)),ncol = k, byrow = T)
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

    # Posterior = clear_Posterior(current_state, c('Lambda','delta','F_h2','F_a','F','tot_Y_prec','resid_h2','E_a','B'))
# r2 = 0
#     Posterior = list(
# 		    Lambda        = matrix(0,nr=0,nc=0),
# 		    F_a           = matrix(0,nr=0,nc=0),
# 		    F             = matrix(0,nr=0,nc=0),
# 		    delta         = matrix(0,nr=0,nc=0),
#             tot_F_prec    = matrix(0,nr=0,nc=0),
#             F_h2          = matrix(0,nr=0,nc=0),
#             tot_Y_prec    = matrix(0,nr = p,nc = 0),
#             resid_h2      = matrix(0,nr = p,nc = 0),
#             resid_Y_prec  = matrix(0,nr = p,nc = 0),
#             E_a_prec      = matrix(0,nr = p,nc = 0),
# 		    B             = matrix(0,nr = b,nc = p),
# 		    E_a           = matrix(0,nr = sum(r_RE),nc = p),
#             W_prec        = matrix(0,nr = p,nc = 0),
# 		    W             = matrix(0,nr = r2,nc = p)
#     	)
    Posterior = clear_Posterior(list(current_state = current_state))$Posterior

# ------------------------------------ #
# ----Precalculate ZAZts, chol_As ---- #
# ------------------------------------ #
    randomEffects = lapply(randomEffects,function(x) Matrix(x))
    names(randomEffects) = RE_names

    ZAZts = list()
    for(i in 1:n_RE){
    	ZAZts[[i]] = tcrossprod(Z_matrices[[i]] %*% randomEffects[[i]],Z_matrices[[i]])
    }
    # inverse_Sigmas = list()
    # for(i in 1:nrow(h2_divisions)){
    # 	Sigma = 0
    # 	for(re in RE_names){
    # 		Sigma = Sigma + h2_divisions[i,re] * randomEffects[[re]]
    # 	}
    # 	Sigma = Matrix(Sigma,sparse=T)
    # 	chol_Sigma = chol(Sigma)
    # 	chol_Sigma_inv = solve(chol_Sigma)
    # 	inverse_Sigmas[[i]]$sqrt_det = prod(diag(chol_Sigma_inv))
    # 	inverse_Sigmas[[i]]$Sigma_inv = chol_Sigma_inv %*% t(chol_Sigma_inv)
    # }

    chol_As = lapply(1:n_RE,function(i) chol(randomEffects[[i]]))
    Ai_mats = lapply(1:n_RE,function(i) chol2inv(chol_As[[i]]))

    ZtZ = crossprod(Z_all)
    randomEffect_Cs = lapply(1:ncol(h2_divisions),function(i) {    	
    	h2s = h2_divisions[,i]
    	h2s = pmax(1e-10,h2s)
		Ai = do.call(bdiag,lapply(1:length(h2s),function(i) Ai_mats[[i]]/h2s[i]))
		C = ZtZ/(1-sum(h2s))
		C = C + Ai
		C
    })

	make_Sigma = function(ZAZts,h2s){
		R = 0
		for(i in 1:length(h2s)){
			R = R + h2s[i]*ZAZts[[i]]
		}
		R + (1-sum(h2s)) * Diagonal(nrow(R))
	}
	
	Sigmas = lapply(1:ncol(h2_divisions),function(i) {
		Sigma = make_Sigma(ZAZts,h2_divisions[,i])
		det = det(Sigma)
		chol = chol(Sigma)
		list(Sigma = Sigma, det = det,chol = chol)
	})


# ----------------------------- #
# ----Save run parameters------ #
# ----------------------------- #

	run_variables = list(
			p               = p,
			n               = n,
			r_RE            = r_RE,
			b               = b,
			Mean_Y          = Mean_Y,
			VY              = VY,
			ZAZts           = ZAZts,
			chol_As         = chol_As,
			Ai_mats         = Ai_mats,
			Sigmas          = Sigmas,
			randomEffect_Cs = randomEffect_Cs
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
