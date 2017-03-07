sample_BSFG.general_BSFG = function(BSFG_state,n_samples,ncores = detectCores(),...) {
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

	sp = (start_i + n_samples - burn)/thin - Posterior$total_samples
	Posterior = expand_Posterior(Posterior,max(0,sp))

	# ----------------------------------------------- #
	# --------------start gibbs sampling------------- #
	# ----------------------------------------------- #

	current_state$E_a_prec = with(current_state,tot_Eta_prec / rowSums(resid_h2))

	start_time = Sys.time()
	for(i in start_i+(1:n_samples)){
		current_state$nrun = i
		current_state_names = names(current_state)
		current_state = within(c(current_state,priors,run_parameters, run_variables,data_matrices), {
			k = ncol(Lambda)

		 # -----Sample Lambda and B ------------------ #
			#conditioning on W, F, marginalizing over random effects (conditional on resid_h2)
			Design = cbind(X,F)
			rows = b + k
			prior_mean = matrix(0,rows,p)
			if(b > 0) {
			  prior_prec = rbind(matrix(prec_B,b,p),t(Plam))
			} else{ # b == 0
			  prior_prec = t(Plam)
			}
			coefs = sample_MME_fixedEffects(Eta,Design,Sigma_Choleskys, Sigma_Perm,  resid_h2_index, tot_Eta_prec, prior_mean, prior_prec,ncores)
			if(b > 0){
				B[] = coefs[1:b,,drop=FALSE]
			}
			Lambda[] = t(coefs[b + 1:k,,drop=FALSE])

	 	# # -----Sample Lambda_prec------------- #
			Lambda2 = Lambda^2
			Lambda_prec = matrix(rgamma(p*k,shape = (Lambda_df + 1)/2,rate = (Lambda_df + sweep(Lambda2,2,tauh,'*'))/2),nr = p,nc = k)

		 # # -----Sample delta, update tauh------ #
			delta[] = sample_delta_c( delta,tauh,Lambda_prec,delta_1_shape,delta_1_rate,delta_2_shape,delta_2_rate,Lambda2,times = 100)
			tauh[]  = matrix(cumprod(delta),nrow=1)

		 # # -----Update Plam-------------------- #
			Plam[] = sweep(Lambda_prec,2,tauh,'*')

		 # -----Sample tot_Eta_prec, resid_h2, E_a ---------------- #
			#conditioning on B, F, Lambda, resid_h2, tot_Eta_prec
			Eta_tilde = Eta - X %*% B - F %*% t(Lambda)
			tot_Eta_prec[] = sample_tot_prec(Eta_tilde, tot_Eta_prec_shape, tot_Eta_prec_rate, Sigma_Choleskys, Sigma_Perm, resid_h2_index,ncores)
			resid_h2_index = sample_h2s_discrete(Eta_tilde,tot_Eta_prec, Sigma_Choleskys, Sigma_Perm, Resid_discrete_priors,ncores)
			resid_h2[] = h2s_matrix[,resid_h2_index,drop=FALSE]
			E_a_prec = tot_Eta_prec / colSums(resid_h2)

			E_a[] = sample_MME_ZAZts(Eta_tilde, Z, tot_Eta_prec, randomEffect_C_Choleskys, resid_h2, resid_h2_index,chol_Ai_mats,ncores)


			# -----Sample Lambda and B_F ------------------ #
			# F, marginalizing over random effects (conditional on F_h2, tot_F_prec)
			if(b_F > 0){
  			prior_mean = matrix(0,b_F,p)
  			prior_prec = matrix(prec_B[-1],b_F,k)
  			B_F = sample_MME_fixedEffects(F,X_F,Sigma_Choleskys, Sigma_Perm, F_h2_index, tot_F_prec, prior_mean, prior_prec,ncores)
  			F_tilde = F - X_F %*% B_F
			} else{
			  F_tilde = F
			}

		 # -----Sample tot_F_prec, F_h2, F_a ---------------- #
			#conditioning on B, F, Lambda, F_h2, tot_F_prec

			tot_F_prec[] = sample_tot_prec(F_tilde, tot_F_prec_shape, tot_F_prec_rate, Sigma_Choleskys, Sigma_Perm, F_h2_index,ncores)
			F_h2_index = sample_h2s_discrete(F_tilde,tot_F_prec, Sigma_Choleskys, Sigma_Perm, F_discrete_priors,ncores)
			F_h2[] = h2s_matrix[,F_h2_index,drop=FALSE]

			F_a[] = sample_MME_ZAZts(F_tilde, Z, tot_F_prec, randomEffect_C_Choleskys, F_h2, F_h2_index,chol_Ai_mats,ncores)

		 # -----Sample F----------------------- #
			#conditioning on B, F_a,E_a,Lambda, F_h2
			Eta_tilde = as.matrix(Eta - X %*% B - Z %*% E_a)
			F_e_prec = tot_F_prec / (1-colSums(F_h2))
			resid_Eta_prec = tot_Eta_prec / (1-colSums(resid_h2))
			if(b_F > 0) {
			  prior_mean = X_F %*% B_F + Z %*% F_a
			} else {
			  prior_mean = Z %*% F_a
			}
			F[] = sample_factors_scores_sparse_c( Eta_tilde, as.matrix(prior_mean),Lambda,resid_Eta_prec,F_e_prec )

		 # -----Sample prec_B------------- #
			if(b > 1) {
			  if(b_F > 0){
			    B2 = cbind(B[-1,],B_F)^2
			  } else{
			    B2 = B[-1,,drop=FALSE]^2
			  }
			  prec_B[1,-1] = rgamma(b-1, shape = fixed_prec_shape + ncol(B2)/2, rate = fixed_prec_rate + rowSums(B2)/2)
			}
	  })
	  current_state = current_state[current_state_names]

	  # ----- sample Eta ----- #
	  data_model_state = run_parameters$data_model(data_matrices$Y,run_parameters$data_model_parameters,current_state,data_matrices)
	  current_state[names(data_model_state)] = data_model_state

	 # -- adapt number of factors to samples ---#
		current_state = update_k( current_state, priors, run_parameters, data_matrices)

	 # -- save sampled values (after thinning) -- #
		if( (i-burn) %% thin == 0 && i > burn) {
			Posterior = save_posterior_sample(current_state, Posterior)
		}
	}
	end_time = Sys.time()
	print(end_time - start_time)
	current_state$total_time = current_state$total_time + end_time - start_time


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
