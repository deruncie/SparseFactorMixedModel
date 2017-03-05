initialize_BSFG.general_BSFG = function(BSFG_state, A_mats = NULL, chol_Ai_mats = NULL,
									ncores = detectCores(),verbose=T,...){

    Y          = BSFG_state$data_matrices$Y
    Y_missing  = BSFG_state$data_matrices$Y_missing
    X_F        = BSFG_state$data_matrices$X_F
    Z_matrices = BSFG_state$data_matrices$Z_matrices
    Z          = BSFG_state$data_matrices$Z
    h2s_matrix = BSFG_state$data_matrices$h2s_matrix

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


# ----------------------------- #
# -----Initialize variables---- #
# ----------------------------- #

  # Initialize Eta
    #    Here, Eta = Y.
    #    With missing data, Eta is complete data
    #    Eta could be parameters of a Y-level model (independent across individuals)
    Eta = Y
    if(sum(Y_missing)>0) {
      Eta[Y_missing] = rnorm(sum(Y_missing))
    }

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

    F_a = lapply(RE_names,function(effect){
    	matrix(rnorm(r_RE[effect] * k, 0, sqrt(F_h2[effect,] / tot_F_prec)),ncol = k, byrow = T)
    })
    names(F_a) = RE_names

    # Factor fixed effects
    B_F = matrix(rnorm(b_F * k),b_F,k)

    F = X_F %*% B_F + matrix(rnorm(n * k, 0, sqrt((1-colSums(F_h2)) / tot_F_prec)),ncol = k, byrow = T)
    for(effect in RE_names) {
    	F = F + Z_matrices[[effect]] %*% F_a[[effect]]
    }
    F = as.matrix(F)
    F_a = do.call(rbind,F_a)

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

    E_a = do.call(rbind,lapply(RE_names,function(effect){
    	matrix(rnorm(r_RE[effect] * p, 0, sqrt(resid_h2[effect,] / tot_Eta_prec)),ncol = p, byrow = T)
    }))
    colnames(E_a) = traitnames

  # Fixed effects
    B = matrix(rnorm(b*p), ncol = p)
    colnames(B) = traitnames

    if(b > 0) {
      prec_B = matrix(c(1e-10,rgamma(b-1,shape = priors$fixed_prec_shape, rate = priors$fixed_prec_rate)),nrow=1)
    } else{
      prec_B = matrix(1e-10,ncol=1,nrow=1)
    }
# ----------------------- #
# ---Save initial values- #
# ----------------------- #
    current_state = list(
        Eta            = Eta,
    		Lambda_prec    = Lambda_prec,
    		delta          = delta,
    		tauh           = tauh,
    		Plam           = Plam,
    		Lambda         = Lambda,
    		tot_F_prec     = tot_F_prec,
    		F_h2_index     = F_h2_index,
    		F_h2           = F_h2,
    		F_a            = F_a,
    		F              = F,
    		tot_Eta_prec     = tot_Eta_prec,
    		resid_h2_index = resid_h2_index,
    		resid_h2       = resid_h2,
    		E_a            = E_a,
    		B              = B,
    		B_F            = B_F,
    		prec_B         = prec_B,
    		traitnames     = traitnames,
    		nrun           = 0,
    		total_time     = 0
    )


# ----------------------- #
# -Initialize Posterior-- #
# ----------------------- #
    Posterior = list(
        sample_params = c('Lambda','F_a','F','delta','tot_F_prec','F_h2','tot_Eta_prec','resid_h2', 'B', 'B_F', 'prec_B'),
        posteriorMean_params = c('E_a'),
        per_trait_params = c('tot_Eta_prec','resid_h2','B'),
        total_samples = 0,
        folder = run_parameters$Posterior_folder,
        files = c()
        # per_trait_params = c('tot_Eta_prec','resid_h2','B','E_a')
    )
    if(run_parameters$save_Eta) Posterior$sample_params = c(Posterior$sample_params,'Eta')
    Posterior = reset_Posterior(Posterior,current_state)

# ------------------------------------ #
# ----Precalculate ZAZts, chol_As ---- #
# ------------------------------------ #
  if(verbose) print('Pre-calculating matrices')

	#C
	# recover()
	ZtZ = crossprod(Z)
	make_Chol_Ai = function(chol_Ai_mats,h2s){
		do.call(bdiag,lapply(1:length(h2s),function(i) {
			if(h2s[i] == 0) return(Diagonal(nrow(chol_Ai_mats[[i]]),Inf))  # if h2==0, then we want a Diagonal matrix with Inf diagonal.
			chol_Ai = chol_Ai_mats[[i]]
			chol_Ai@x = chol_Ai@x / sqrt(h2s[i])
			chol_Ai
		}))
	}

	# setup of symbolic Cholesky of C
	print('creating randomEffects_C')

	Ai = forceSymmetric(crossprod(make_Chol_Ai(chol_Ai_mats,rep(1,n_RE)/(n_RE+1))))
	Cholesky_C = Cholesky(ZtZ + Ai)

	randomEffect_C_Choleskys = mclapply(1:ncol(h2s_matrix),function(i) {
		if(i %% 100 == 0 && verbose) print(sprintf('randomEffects_C %d of %d',i,ncol(h2s_matrix)))
		h2s = h2s_matrix[,i]
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
	print('creating Sigma_Choleskys')
	Sigma = make_Sigma(ZAZts,h2s_matrix[,2])
	Cholesky_Sigma_base = Cholesky(Sigma,perm=T,super=T)
	stopifnot(!isLDL(Cholesky_Sigma_base))
	Sigma_Perm = expand(Cholesky_Sigma_base)$P
	if(all(diag(Sigma_Perm))) Sigma_Perm = NULL

	Sigma_Choleskys = mclapply(1:ncol(h2s_matrix),function(i) {
		if(i %% 100 == 0 && verbose) print(sprintf('Sigma_Choleskys %d of %d',i,ncol(h2s_matrix)))
		Sigma = forceSymmetric(make_Sigma(ZAZts,h2s_matrix[,i]))
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
			chol_Ai_mats             = chol_Ai_mats
    ))

    RNG = list(
    	Random.seed = .Random.seed,
    	RNGkind = RNGkind()
    )

    BSFG_state$run_variables = run_variables
    BSFG_state$RNG = RNG
    BSFG_state$Posterior = Posterior
    BSFG_state$current_state = current_state

    return(BSFG_state)
}
