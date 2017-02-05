fast_BSFG_ipx_sampler = function(BSFG_state,n_samples) {
	# -- Daniel Runcie -- #
	#
	# Gibbs sampler for genetic covariance estimation based on mixed effects
	# model, with missing data
	#
	# Code for:
	# 
	# Runcie and Mukherjee (2013) Dissecting high-dimensional traits
	# with Bayesian sparse factor analysis of genetic covariance matrices.
	# GENETICS.
	#
	# (c) July 30, 2013
	#
	# code based on original provided by Anirban Bhattacharya
	#
	#         
	# This function implements the BSF-G partially collapsed Gibbs sampler.
	# Variable notation follows the Runcie and Mukherjee manuscript as closely
	# as possible.
	#
	# All input and initialization functions are carried out by the function
	# fast_BSFG_sampler_init See the documentation of that function for details.
	# 
	# The sampler is designed to do a short-medium length run and then return
	# the state of the chain. After assessing the progress of the chain,
	# (optionally), the Posterior matrices can be reset, and the chain
	# re-started where it left off. Right now, the random seed is not saved. Probably should
	# be if I can figure out how.
	# 
	# This function takes the following inputs:
	#     data_matrices: struct holding the (should be imutable) data, design and incidence matrices:
	#          Y:  		Full phenotype data. n x p
	# 		   X: 		Fixed effect design matrix. n x b
	#          Z_1:     Random effect 1 incidence matrix. n x r
	#          Y_missing: incidence matrix of missing data in Y
	#     start_i: iteration number of the end of the last run.
	#     draw_iter: frequency of updating diagnostic plots
	#     burn: number of burnin samples
	#     sp: total number of samples to collect
	#     thin: thinning rate of chain
	#     simulation: boolean. Is this a simulation?
	#     params: struct with chain parameters.
	#     priors: struct with all relevant prior hyperparameters
	#     Posterior: struct with posterior matrices, or running posterior means. 
	#            Note: not sure how well posterior means work after re-starting chain. Probably not well.
	#     current_state: current (initial) conditions of all model parameters
	# 
	# Several diagnostic plots are produced during the run. 
	#     Their interpretation is described within the source codes:
	#         draw_simulation_diagnostics.m: For simulated data with known true values
	#         draw_results_diagnostics.m: Otherwise
	#         

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
		Posterior$Lambda        = cbind(Posterior$Lambda,matrix(0,nr = nrow(Posterior$Lambda),nc = sp))
		Posterior$F             = cbind(Posterior$F,matrix(0,nr = nrow(Posterior$F),nc = sp))
		Posterior$F_a           = cbind(Posterior$F_a,matrix(0,nr = nrow(Posterior$F_a),nc = sp))
		Posterior$delta         = cbind(Posterior$delta,matrix(0,nr = nrow(Posterior$delta),nc = sp))
		Posterior$tot_Y_prec    = cbind(Posterior$tot_Y_prec,matrix(0,nr = nrow(Posterior$tot_Y_prec),nc = sp))
		Posterior$resid_h2      = cbind(Posterior$resid_h2  ,matrix(0,nr = nrow(Posterior$resid_h2  ),nc = sp))
		Posterior$tot_F_prec    = cbind(Posterior$tot_F_prec,matrix(0,nr = nrow(Posterior$tot_F_prec),nc = sp))
		Posterior$F_h2          = cbind(Posterior$F_h2,matrix(0,nr = nrow(Posterior$F_h2),nc = sp))
		Posterior$resid_Y_prec  = cbind(Posterior$resid_Y_prec,matrix(0,nr = nrow(Posterior$tot_Y_prec),nc = sp))
		Posterior$E_a_prec      = cbind(Posterior$E_a_prec,matrix(0,nr = nrow(Posterior$tot_Y_prec),nc = sp))
	}

	# ----------------------------------------------- #
	# --------------start gibbs sampling------------- #
	# ----------------------------------------------- #

	start_time = Sys.time()
	for(i in start_i+(1:n_samples)){
		current_state$nrun = i
		current_state_names = names(current_state)
		current_state = within(c(current_state,priors,run_parameters, run_variables,data_matrices), {
			k = ncol(Lambda)

		 # -----fill in missing phenotypes----- #
			#conditioning on everything else
			if(sum(Y_missing)>0) {
				meanTraits = X %*% B + F %*% t(Lambda) + Z %*% E_a
				resids = matrix(rnorm(p*n,0,sqrt(1/resid_Y_prec)),nr = n,nc = p,byrow=T)
				Y[Y_missing] = meanTraits[Y_missing] + resids[Y_missing]
			}
			  
		 # -----Sample Lambda and B ------------------ #
			#conditioning on F, marginalizing over E_a

			Design = as.matrix(cbind(X,F))
			rows = b + k
			prior_mean = matrix(0,rows,p)
			prior_prec = rbind(0,t(Plam))
			coefs = sample_coefs_parallel_sparse_c( Y,Design,resid_h2, tot_Y_prec,prior_mean,prior_prec,invert_aI_bZAZ,1)
			
			if(b > 0){
				B = coefs[1:b,]
			}
			Lambda = t(coefs[b + 1:k,])

	 	# # -----Sample Lambda_prec------------- #
			Lambda2 = Lambda^2
			Lambda_prec = matrix(rgamma(p*k,shape = (Lambda_df + 1)/2,rate = (Lambda_df + sweep(Lambda2,2,tauh,'*'))/2),nr = p,nc = k)

		 # # -----Sample delta, update tauh------ #
			# delta = sample_delta_ipx( delta,tauh,Lambda_prec,delta_1_shape,delta_1_rate,delta_2_shape,delta_2_rate,Lambda2,times = 100)
			delta = sample_delta_c( delta,tauh,Lambda_prec,delta_1_shape,delta_1_rate,delta_2_shape,delta_2_rate,Lambda2,times = 100)
			tauh  = cumprod(delta)
			
		 # # -----Update Plam-------------------- #
			Plam = sweep(Lambda_prec,2,tauh,'*')

		 # -----Sample resid_h2, tot_Y_prec, E_a ---------------- #
			#conditioning on W, B, F, Lambda, marginalizing over E_a 
			Y_tilde = as.matrix(Y - X %*% B - F %*% t(Lambda))
			tot_Y_prec = sample_tot_prec_sparse_c(Y_tilde,resid_h2,tot_Y_prec_shape,tot_Y_prec_rate,invert_aI_bZAZ)
			
			resid_h2_index = sample_h2s_discrete_given_p_sparse_c(Y_tilde,h2_divisions,h2_priors_resids,tot_Y_prec,invert_aI_bZAZ)
			resid_h2 = h2s_matrix[,resid_h2_index,drop=FALSE]

			E_a = sample_randomEffects_parallel_sparse_c( Y_tilde, Z, tot_Y_prec, resid_h2, invert_aZZt_Ainv, 1)

			resid_Y_prec = tot_Y_prec / (1-resid_h2)
			
		 # -----Sample F_h2 and tot_F_prec, F_a -------------------- #
			#conditioning on F, F_a
			tot_F_prec = sample_tot_prec_sparse_c(F,F_h2,tot_F_prec_shape,tot_F_prec_rate,invert_aI_bZAZ)

			F_h2_index = sample_h2s_discrete_given_p_sparse_c(F,h2_divisions,h2_priors_factors,tot_F_prec,invert_aI_bZAZ)
			F_h2 = h2s_matrix[,F_h2_index,drop=FALSE]

	    	F_a = sample_randomEffects_parallel_sparse_c(F,Z,tot_F_prec, F_h2, invert_aZZt_Ainv, 1)
			
		 # -----Sample F----------------------- #
			#conditioning on B, F_a,E_a,W,Lambda, F_h2
			Y_tilde = as.matrix(Y - X %*% B - Z %*% E_a)
			F_e_prec = tot_F_prec / (1-F_h2)
			F = sample_factors_scores_ipx_sparse_c( Y_tilde, Z,Lambda,resid_Y_prec,F_a,F_e_prec )
						
	})
	current_state = current_state[current_state_names]

	 # -- adapt number of factors to samples ---#
		current_state = update_k( current_state, priors, run_parameters, data_matrices)

	 # -- save sampled values (after thinning) -- #
		if( (i-burn) %% thin == 0 && i > burn) {
				
			sp_num = (i-burn)/thin
			
			Posterior = save_posterior_samples( sp_num,current_state, Posterior)
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