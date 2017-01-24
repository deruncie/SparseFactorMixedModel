require(parallel)
BSFG_discreteRandom_sampler = function(BSFG_state,n_samples,ncores = detectCores()) {

	data_matrices  = BSFG_state$data_matrices
	params         = BSFG_state$params
	priors         = BSFG_state$priors
	Posterior      = BSFG_state$Posterior
	current_state  = BSFG_state$current_state
	run_parameters = BSFG_state$run_parameters
	run_variables  = BSFG_state$run_variables


	# ----------------------------------------------- #
	# -----------Reset Global Random Number Stream--- #
	# ----------------------------------------------- #
	do.call("RNGkind",as.list(BSFG_state$RNG$RNGkind))  ## must be first!
	assign(".Random.seed", BSFG_state$RNG$Random.seed, .GlobalEnv)

	# ----------------------------------------------- #
	# ----------------Set up run--------------------- #
	# ----------------------------------------------- #
	save_freq    = run_parameters$save_freq
	burn         = run_parameters$burn
	thin         = run_parameters$thin
	start_i      = current_state$nrun

	# ----------------------------------------------- #
	# ---Extend posterior matrices for new samples--- #
	# ----------------------------------------------- #

	sp = (start_i + n_samples - burn)/thin - ncol(Posterior$Lambda)
	if(sp > 0){
		for(param in names(Posterior)){
			if(param %in% c('B','E_a')) next
			Posterior[[param]] = cbind(Posterior[[param]],matrix(0,nr = nrow(Posterior[[param]]),nc = sp))
		}
		# Posterior$Lambda        = cbind(Posterior$Lambda,matrix(0,nr = nrow(Posterior$Lambda),nc = sp))
		# Posterior$F             = cbind(Posterior$F,matrix(0,nr = nrow(Posterior$F),nc = sp))
		# Posterior$F_a           = cbind(Posterior$F_a,matrix(0,nr = nrow(Posterior$F_a),nc = sp))
		# Posterior$delta         = cbind(Posterior$delta,matrix(0,nr = nrow(Posterior$delta),nc = sp))
		# Posterior$tot_Y_prec    = cbind(Posterior$tot_Y_prec,matrix(0,nr = nrow(Posterior$tot_Y_prec),nc = sp))
		# Posterior$resid_h2      = cbind(Posterior$resid_h2  ,matrix(0,nr = nrow(Posterior$resid_h2  ),nc = sp))
		# Posterior$tot_F_prec    = cbind(Posterior$tot_F_prec,matrix(0,nr = nrow(Posterior$tot_F_prec),nc = sp))
		# Posterior$F_h2          = cbind(Posterior$F_h2,matrix(0,nr = nrow(Posterior$F_h2),nc = sp))
		# Posterior$resid_Y_prec  = cbind(Posterior$resid_Y_prec,matrix(0,nr = nrow(Posterior$tot_Y_prec),nc = sp))
		# Posterior$E_a_prec      = cbind(Posterior$E_a_prec,matrix(0,nr = nrow(Posterior$tot_Y_prec),nc = sp))
	}

	# ----------------------------------------------- #
	# --------------start gibbs sampling------------- #
	# ----------------------------------------------- #

	current_state$E_a_prec = with(current_state,tot_Y_prec / rowSums(resid_h2))

	start_time = Sys.time()
	for(i in start_i+(1:n_samples)){
		current_state$nrun = i
		current_state_names = names(current_state)
		current_state = within(c(current_state,priors,run_parameters, run_variables,data_matrices), {
			k = ncol(Lambda)

		 # -----fill in missing phenotypes----- #
			#conditioning on everything else
			# if(sum(Y_missing)>0) {
			# 	meanTraits = X %*% B + F %*% t(Lambda) + Z_1_sparse %*% E_a + Z_2_sparse %*% W
			# 	resids = matrix(rnorm(p*n,0,sqrt(1/resid_Y_prec)),nr = n,nc = p,byrow=T)
			# 	Y[Y_missing] = meanTraits[Y_missing] + resids[Y_missing]
			# }
			# recover()
			  
		 # -----Sample Lambda and B ------------------ #
			#conditioning on W, F, marginalizing over random effects (conditional on resid_h2)
			Design = cbind(X,F)
			rows = b + k
			prior_mean = matrix(0,rows,p)
			prior_prec = rbind(1e-6,t(Plam))
			coefs = sample_MME_fixedEffects(Y,Design,Sigmas, resid_h2_index, tot_Y_prec, prior_mean, prior_prec,ncores)
			# coefs = sample_coefs_DR_parallel_sparse_c( T, Design, tot_Y_prec, prior_mean, prior_prec, inverse_Sigmas, resid_h2_index, 1)

			if(b > 0){
				B = coefs[1:b,]
			}
			Lambda = t(coefs[b + 1:k,])

	 	# # -----Sample Lambda_prec------------- #
			Lambda2 = Lambda^2
			Lambda_prec = matrix(rgamma(p*k,shape = (Lambda_df + 1)/2,rate = (Lambda_df + sweep(Lambda2,2,tauh,'*'))/2),nr = p,nc = k)

		 # # -----Sample delta, update tauh------ #
			delta = sample_delta_c( delta,tauh,Lambda_prec,delta_1_shape,delta_1_rate,delta_2_shape,delta_2_rate,Lambda2,times = 100)
			tauh  = cumprod(delta)
			
		 # # -----Update Plam-------------------- #
			Plam = sweep(Lambda_prec,2,tauh,'*')

		 # -----Sample tot_Y_prec, resid_h2, E_a ---------------- #
			#conditioning on B, F, Lambda, resid_h2, tot_Y_prec
			Y_tilde = Y - X %*% B - F %*% t(Lambda)
			tot_Y_prec = sample_tot_prec(Y_tilde, tot_Y_prec_shape, tot_Y_prec_rate, Sigmas, resid_h2_index,ncores)
			resid_h2_index = sample_h2s_discrete(Y_tilde,tot_Y_prec, Sigmas,Resid_discrete_priors,ncores)
			resid_h2 = h2_divisions[,resid_h2_index,drop=FALSE]
			E_a_prec = tot_Y_prec / colSums(resid_h2)

			# Y_tilde = as.matrix(Y - X %*% B - F %*% t(Lambda))
			# tot_Y_prec = sample_tot_prec_sparse_c(Y_tilde, tot_Y_prec_shape, tot_Y_prec_rate, inverse_Sigmas, resid_h2_index)

			# resid_h2_index = sample_discrete_h2s_given_p_sparse_c( Y_tilde, tot_Y_prec, Resid_discrete_priors, inverse_Sigmas)
			# resid_h2 = h2_divisions[resid_h2_index,]

			prior_mean = matrix(0,sum(r_RE),p)
			E_a = sample_MME_ZAZts(Y_tilde, Z_all, tot_Y_prec, prior_mean, randomEffect_Cs, Ai_mats, resid_h2, resid_h2_index,chol_As,ncores)
			
		 # -----Sample tot_F_prec, F_h2, F_a ---------------- #
			#conditioning on B, F, Lambda, F_h2, tot_F_prec

			tot_F_prec = sample_tot_prec(F, tot_F_prec_shape, tot_F_prec_rate, Sigmas, F_h2_index,ncores)
			F_h2_index = sample_h2s_discrete(F,tot_F_prec, Sigmas,F_discrete_priors,ncores)
			F_h2 = h2_divisions[,F_h2_index,drop=FALSE]

			# tot_F_prec = sample_tot_prec_sparse_c(F, tot_F_prec_shape, tot_F_prec_rate, inverse_Sigmas, F_h2_index)

			# F_h2_index = sample_discrete_h2s_given_p_sparse_c( F, tot_F_prec, F_discrete_priors, inverse_Sigmas)
			# F_h2 = h2_divisions[F_h2_index,]

			prior_mean = matrix(0,sum(r_RE),k)
			F_a = sample_MME_ZAZts(F, Z_all, tot_F_prec, prior_mean, randomEffect_Cs, Ai_mats, F_h2, F_h2_index,chol_As,ncores)

		 # -----Sample F----------------------- #
			#conditioning on B, F_a,E_a,Lambda, F_h2
			Y_tilde = as.matrix(Y - X %*% B - Z_all %*% E_a)
			F_e_prec = tot_F_prec / (1-colSums(F_h2))
			resid_Y_prec = tot_Y_prec / (1-colSums(resid_h2))
			F = sample_factors_scores_ipx_sparse_c( Y_tilde, Z_all,Lambda,resid_Y_prec,F_a,F_e_prec )
	
	})
	current_state = current_state[current_state_names]

	 # -- adapt number of factors to samples ---#
		current_state = update_k( current_state, priors, run_parameters, data_matrices)

	 # -- save sampled values (after thinning) -- #
		if( (i-burn) %% thin == 0 && i > burn) {
				
			sp_num = (i-burn)/thin
			
			Posterior = save_posterior_samples( sp_num,current_state, Posterior)
			
			if(sp_num %% save_freq == 0) save(Posterior,file='Posterior.RData')
		}
	}
	end_time = Sys.time()
	print(end_time - start_time)
	save(Posterior,file = 'Posterior.RData')


	# ----------------------------------------------- #
	# ------------Save state for restart------------- #
	# ----------------------------------------------- #


	save(current_state,file='current_state.RData')


	BSFG_state$current_state = current_state
	BSFG_state$Posterior = Posterior
	BSFG_state$RNG = list(
		Random.seed = .Random.seed,
		RNGkind = RNGkind()
	)
	return(BSFG_state)
}