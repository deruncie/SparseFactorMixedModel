initialize_BSFG.general_BSFG = function(BSFG_state, K_mats = NULL, chol_Ki_mats = NULL,
									ncores = detectCores(),verbose=T,...){

    Y          = BSFG_state$data_matrices$Y
    X_F        = BSFG_state$data_matrices$X_F
    Z_matrices = BSFG_state$data_matrices$Z_matrices
    Z          = BSFG_state$data_matrices$Z
    h2s_matrix = BSFG_state$data_matrices$h2s_matrix
    cis_effects_index = BSFG_state$data_matrices$cis_effects_index

    RE_names   = rownames(h2s_matrix)
    n_RE       = length(RE_names)

    traitnames     = BSFG_state$traitnames
    priors         = BSFG_state$priors
    run_parameters = BSFG_state$run_parameters
    run_variables  = BSFG_state$run_variables

    p    = BSFG_state$run_variables$p
    n    = BSFG_state$run_variables$n
    r_RE = BSFG_state$run_variables$r_RE
    b    = BSFG_state$run_variables$b
    b_F  = BSFG_state$run_variables$b_F

    # cis effects
    cis_effects = matrix(rnorm(length(cis_effects_index),0,1),nrow=1)


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
    delta = with(priors,matrix(c(rgamma(1,shape = delta_1_shape,rate = delta_1_rate),rgamma(k-1,shape = delta_2_shape,rate = delta_2_rate)),nrow=1))
    tauh  = matrix(cumprod(delta),nrow=1)

    # Total Factor loading precisions Lambda_prec * tauh
    Plam = sweep(Lambda_prec,2,tauh,'*')

    # Lambda - factor loadings
     #   Prior: Normal distribution for each element.
     #       mu = 0
     #       sd = sqrt(1/Plam)
    Lambda = matrix(rnorm(p*k,0,sqrt(1/Plam)),nr = p,nc = k)
    rownames(Lambda) = traitnames

  # Factor scores:
     # p-vector of factor precisions. Note - this is a 'redundant' parameter designed to give the Gibbs sampler more flexibility
	 #  Prior: Gamma distribution for each element
     #       shape = tot_F_prec_shape
     #       rate = tot_F_prec_rate
    tot_F_prec = with(priors,matrix(rgamma(k,shape = tot_F_prec_shape,rate = tot_F_prec_rate),nrow=1))

    # Factor discrete variances
     # k-matrix of n_RE x k with
    F_h2_index = sample(1:ncol(h2s_matrix),k,replace=T)
    F_h2 = h2s_matrix[,F_h2_index,drop=FALSE]

    U_F = lapply(RE_names,function(effect){
    	matrix(rnorm(r_RE[effect] * k, 0, sqrt(F_h2[effect,] / tot_F_prec)),ncol = k, byrow = T)
    })
    names(U_F) = RE_names

    # Factor fixed effects
    B_F = matrix(rnorm(b_F * k),b_F,k)

    F = X_F %*% B_F + matrix(rnorm(n * k, 0, sqrt((1-colSums(F_h2)) / tot_F_prec)),ncol = k, byrow = T)
    for(effect in RE_names) {
    	F = F + Z_matrices[[effect]] %*% U_F[[effect]]
    }
    F = as.matrix(F)
    U_F = do.call(rbind,U_F)
    rownames(U_F) = colnames(Z)

  # residuals
     # p-vector of factor precisions. Note - this is a 'redundant' parameter designed to give the Gibbs sampler more flexibility
	 #  Prior: Gamma distribution for each element
     #       shape = tot_Eta_prec_shape
     #       rate = tot_Eta_prec_rate
    tot_Eta_prec = with(priors,matrix(rgamma(p,shape = tot_Eta_prec_shape,rate = tot_Eta_prec_rate),nrow = 1))
    colnames(tot_Eta_prec) = traitnames

    # Resid discrete variances
     # p-matrix of n_RE x p with
    resid_h2_index = sample(1:ncol(h2s_matrix),p,replace=T)
    resid_h2 = h2s_matrix[,resid_h2_index,drop=FALSE]

    U_R = do.call(rbind,lapply(RE_names,function(effect){
    	matrix(rnorm(r_RE[effect] * p, 0, sqrt(resid_h2[effect,] / tot_Eta_prec)),ncol = p, byrow = T)
    }))
    colnames(U_R) = traitnames
    rownames(U_R) = colnames(Z)

  # Fixed effects
    B = matrix(rnorm(b*p), ncol = p)
    colnames(B) = traitnames

    if(b > 0) {
      tau_B = matrix(c(1e-10,rgamma(b-1,shape = priors$fixed_prec_shape, rate = priors$fixed_prec_rate)),nrow=1)
      tau_B_F = tau_B[1,-1,drop=FALSE]
    } else{
      tau_B = matrix(0,ncol=0,nrow=1)
      tau_B_F = tau_B
    }
    prec_B = matrix(tau_B,nrow = b, ncol = p)
    prec_B_F = matrix(tau_B_F,nrow = b_F, ncol = k)

# ----------------------- #
# ---Save initial values- #
# ----------------------- #
    current_state = list(
        k              = k,
    		Lambda_prec    = Lambda_prec,
    		delta          = delta,
    		tauh           = tauh,
    		Plam           = Plam,
    		Lambda         = Lambda,
    		tot_F_prec     = tot_F_prec,
    		F_h2_index     = F_h2_index,
    		F_h2           = F_h2,
    		U_F            = U_F,
    		F              = F,
    		tot_Eta_prec     = tot_Eta_prec,
    		resid_h2_index = resid_h2_index,
    		resid_h2       = resid_h2,
    		U_R            = U_R,
    		B              = B,
    		B_F            = B_F,
    		tau_B          = tau_B,
    		tau_B_F        = tau_B_F,
    		prec_B         = prec_B,
    		prec_B_F       = prec_B_F,
    		cis_effects    = cis_effects,
    		traitnames     = traitnames,
    		nrun           = 0,
    		total_time     = 0
    )
    BSFG_state$current_state = current_state

# ------------------------------------ #
# ----Precalculate ZKZts, chol_Ks ---- #
# ------------------------------------ #
  if(verbose) print('Pre-calculating matrices')

	#C
	# recover()
	ZtZ = crossprod(Z)
	make_Chol_Ki = function(chol_Ki_mats,h2s){
		do.call(bdiag,lapply(1:length(h2s),function(i) {
			if(h2s[i] == 0) return(Diagonal(nrow(chol_Ki_mats[[i]]),Inf))  # if h2==0, then we want a Diagonal matrix with Inf diagonal.
			chol_Ki = chol_Ki_mats[[i]]
			chol_Ki@x = chol_Ki@x / sqrt(h2s[i])
			chol_Ki
		}))
	}

	# setup of symbolic Cholesky of C
	print('creating randomEffects_C')

	Ki = forceSymmetric(crossprod(make_Chol_Ki(chol_Ki_mats,rep(1,n_RE)/(n_RE+1))))
	Cholesky_C = Cholesky(ZtZ + Ki)

	randomEffect_C_Choleskys = mclapply(1:ncol(h2s_matrix),function(i) {
		if(i %% 100 == 0 && verbose) print(sprintf('randomEffects_C %d of %d',i,ncol(h2s_matrix)))
		h2s = h2s_matrix[,i]
		chol_K_inv = make_Chol_Ki(chol_Ki_mats,h2s)
		Ki = crossprod(chol_K_inv)
		C = ZtZ/(1-sum(h2s))
		C = C + Ki
		Cholesky_C_i = update(Cholesky_C,forceSymmetric(C))

		return(list(Cholesky_C = Cholesky_C_i, chol_K_inv = chol_K_inv))
	},mc.cores = ncores)

	# Sigma
	ZKZts = list()
	for(i in 1:n_RE){
		ZKZts[[i]] = forceSymmetric(Z_matrices[[i]] %*% K_mats[[i]] %*% t(Z_matrices[[i]]))
	}


	make_Sigma = function(ZKZts,h2s){
		R = 0
		for(i in 1:length(h2s)){
			R = R + h2s[i]*ZKZts[[i]]
		}
		forceSymmetric(R + (1-sum(h2s)) * Diagonal(nrow(R)))
	}

	# setup of symbolic Cholesky of Sigma
	print('creating Sigma_Choleskys')
	Sigma = make_Sigma(ZKZts,h2s_matrix[,2])
	Cholesky_Sigma_base = Cholesky(Sigma,perm=T,super=T)
	stopifnot(!isLDL(Cholesky_Sigma_base))
	Sigma_Perm = expand(Cholesky_Sigma_base)$P
	if(all(diag(Sigma_Perm))) Sigma_Perm = NULL

	Sigma_Choleskys = mclapply(1:ncol(h2s_matrix),function(i) {
		if(i %% 100 == 0 && verbose) print(sprintf('Sigma_Choleskys %d of %d',i,ncol(h2s_matrix)))
		Sigma = forceSymmetric(make_Sigma(ZKZts,h2s_matrix[,i]))
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

	run_variables = c(run_variables,list(
			Sigma_Choleskys          = Sigma_Choleskys,
			Sigma_Perm               = Sigma_Perm,
			randomEffect_C_Choleskys = randomEffect_C_Choleskys,
			chol_Ki_mats             = chol_Ki_mats
    ))

    RNG = list(
    	Random.seed = .Random.seed,
    	RNGkind = RNGkind()
    )

    BSFG_state$run_variables = run_variables
    BSFG_state$RNG = RNG

    return(BSFG_state)
}
