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
	#          Z_1:     Random effect 1 incidence matrix. n x r1
	#          Z_2:     Random effect 2 incidence matrix. n x r2
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
	# ----------------Load data matrices------------- #
	# ----------------------------------------------- #

	Y  			= data_matrices$Y
	Z_1     	= data_matrices$Z_1
	Z_2     	= data_matrices$Z_2
	X       	= data_matrices$X
	Y_missing 	= data_matrices$Y_missing

	p   = run_variables$p
	n   = run_variables$n
	r   = run_variables$r
	r2  = run_variables$r2
	b   = run_variables$b


	# ----------------------------------------------- #
	# ----------------Load priors-------------------- #
	# ----------------------------------------------- #

	resid_Y_prec_shape =   priors$resid_Y_prec_shape
	resid_Y_prec_rate  =   priors$resid_Y_prec_rate
	E_a_prec_shape     =   priors$E_a_prec_shape
	E_a_prec_rate      =   priors$E_a_prec_rate
	F_a_prec_shape     =   priors$F_a_prec_shape
	F_a_prec_rate      =   priors$F_a_prec_rate
	F_e_prec_shape     =   priors$F_e_prec_shape
	F_e_prec_rate      =   priors$F_e_prec_rate
	tot_Y_prec_shape     =   priors$tot_Y_prec_shape
	tot_Y_prec_rate      =   priors$tot_Y_prec_rate
	tot_F_prec_shape     =   priors$tot_F_prec_shape
	tot_F_prec_rate      =   priors$tot_F_prec_rate
	W_prec_shape       =   priors$W_prec_shape
	W_prec_rate        =   priors$W_prec_rate
	Lambda_df          =   priors$Lambda_df
	delta_1_shape      =   priors$delta_1_shape
	delta_1_rate       =   priors$delta_1_rate
	delta_2_shape      =   priors$delta_2_shape
	delta_2_rate       =   priors$delta_2_rate
	h2_priors_factors  =   priors$h2_priors_factors
	h2_priors_resids   =   priors$h2_priors_resids


	# ----------------------------------------------- #
	# ----------------Load current state------------- #
	# ----------------------------------------------- #
	resid_Y_prec =   current_state$resid_Y_prec
	F            =   current_state$F
	Lambda       =   current_state$Lambda
	E_a          =   current_state$E_a
	W            =   current_state$W
	Lambda_prec  =   current_state$Lambda_prec
	delta        =   current_state$delta
	tauh         =   current_state$tauh
	E_a_prec     =   current_state$E_a_prec
	W_prec       =   current_state$W_prec
	Plam         =   current_state$Plam
	F_a_prec     =   current_state$F_a_prec
	F_e_prec     =   current_state$F_e_prec
	F_h2    	 =   current_state$F_h2
	F_a          =   current_state$F_a
	B            =   current_state$B
	start_i      =   current_state$nrun
	k = ncol(F)

	# ----------------------------------------------- #
	# -----------Reset Global Random Number Stream--- #
	# ----------------------------------------------- #
	do.call("RNGkind",as.list(BSFG_state$RNG$RNGkind))  ## must be first!
	assign(".Random.seed", BSFG_state$RNG$Random.seed, .GlobalEnv)

	# ----------------------------------------------- #
	# -----------Load pre-calcualted matrices-------- #
	# ----------------------------------------------- #
	Ainv                             = run_variables$Ainv
	A_2_inv                          = run_variables$A_2_inv
	chol_Ainv                        = run_variables$chol_Ainv
	chol_A_2_inv                     = run_variables$A_2_inv
	invert_aI_bZAZ                   = run_variables$invert_aI_bZAZ  
	invert_aPXA_bDesignDesignT       = run_variables$invert_aPXA_bDesignDesignT 
	invert_aZZt_Ainv                 = run_variables$invert_aZZt_Ainv
	invert_aPXA_bDesignDesignT_rand2 = run_variables$invert_aPXA_bDesignDesignT_rand2  


	# ----------------------------------------------- #
	# ----------------Set up run--------------------- #
	# ----------------------------------------------- #
	#     b0,b1: parameters controlling rate of adaptation of factor model size
	#     h2_divisions: number of discrete steps for each factor heritability parameter
	#     epsilon: truncation point for factor loadings during adaptation
	b0           = run_parameters$b0
	b1           = run_parameters$b1
	h2_divisions = run_parameters$h2_divisions
	epsilon      = run_parameters$epsilon
	prop         = run_parameters$prop
	save_freq    = run_parameters$save_freq
	burn         = run_parameters$burn
	thin         = run_parameters$thin


	# ----------------------------------------------- #
	# ---Extend posterior matrices for new samples--- #
	# ----------------------------------------------- #

	sp = (start_i + n_samples - burn)/thin - ncol(Posterior$Lambda)
	if(sp > 0){
		Posterior$Lambda        = cbind(Posterior$Lambda,matrix(0,nr = nrow(Posterior$Lambda),nc = sp))
		Posterior$F             = cbind(Posterior$F,matrix(0,nr = nrow(Posterior$F),nc = sp))
		Posterior$F_a           = cbind(Posterior$F_a,matrix(0,nr = nrow(Posterior$F_a),nc = sp))
		Posterior$delta         = cbind(Posterior$delta,matrix(0,nr = nrow(Posterior$delta),nc = sp))
		Posterior$F_h2          = cbind(Posterior$F_h2,matrix(0,nr = nrow(Posterior$F_h2),nc = sp))
		Posterior$resid_Y_prec  = cbind(Posterior$resid_Y_prec,matrix(0,nr = nrow(Posterior$resid_Y_prec),nc = sp))
		Posterior$E_a_prec      = cbind(Posterior$E_a_prec,matrix(0,nr = nrow(Posterior$E_a_prec),nc = sp))
		Posterior$W_prec        = cbind(Posterior$W_prec,matrix(0,nr = nrow(Posterior$W_prec),nc = sp))
	}

	# ----------------------------------------------- #
	# --------------start gibbs sampling------------- #
	# ----------------------------------------------- #

	start_time = Sys.time()
	for(i in start_i+(1:n_samples)){
		current_state$nrun = i
		current_state_names = names(current_state)
		current_state = within(current_state, {
			k = ncol(Lambda)

		 # -----fill in missing phenotypes----- #
			#conditioning on everything else
			if(sum(Y_missing)>0) {
				meanTraits = X %*% B + F %*% t(Lambda) + Z_1 %*% E_a + Z_2 %*% W
				resids = matrix(rnorm(p*n,0,sqrt(1/resid_Y_prec)),nr = n,nc = p,byrow=T)
				Y[Y_missing] = meanTraits[Y_missing] + resids[Y_missing]
			}
			# recover()
			  
		 # -----Sample Lambda------------------ #
			#conditioning on W, B, F, marginalizing over E_a
			Y_tilde = Y - X %*% B - Z_2 %*% W
			# Lambda = sample_Lambda( Y_tilde,F,resid_Y_prec, E_a_prec,Plam,invert_aI_bZAZ )
			# Lambda = sample_Lambda_c( Y_tilde,F,resid_Y_prec, E_a_prec,Plam,invert_aI_bZAZ )
			Lambda = sample_Lambda_parallel_c( Y_tilde,F,resid_Y_prec, E_a_prec,Plam,invert_aI_bZAZ,1)

	 	# # -----Sample Lambda_prec------------- #
			Lambda2 = Lambda^2
			Lambda_prec = matrix(rgamma(p*k,shape = (Lambda_df + 1)/2,rate = (Lambda_df + sweep(Lambda2,2,tauh,'*'))/2),nr = p,nc = k)

		 # # -----Sample delta, update tauh------ #
			delta = sample_delta( delta,tauh,Lambda_prec,delta_1_shape,delta_1_rate,delta_2_shape,delta_2_rate,Lambda2 )
			tauh  = cumprod(delta)
			
		 # # -----Update Plam-------------------- #
			Plam = sweep(Lambda_prec,2,tauh,'*')

		 # -----Sample E_a_prec, tot_Y_prec ---------------- #
			#conditioning on W, B, F, marginalizing over E_a 
			Y_tilde = Y - X %*% B - F %*% t(Lambda)  - Z_2 %*% W

			resid_h2 = resid_Y_prec/(resid_Y_prec + E_a_prec)
			U = invert_aI_bZAZ$U
			s = invert_aI_bZAZ$s
			UtY_tilde = t(U) %*% Y_tilde
			n = nrow(Y)
			tot_Y_prec = sapply(1:p,function(j) {
				Sigma_sqrt = sqrt(resid_h2[j]*s + (1-resid_h2[j]))
				SiUtY_tilde_j = UtY_tilde[,j] / Sigma_sqrt
				rgamma(1,shape = tot_Y_prec_shape + n/2, rate = tot_Y_prec_rate+1/2*sum(SiUtY_tilde_j^2))
			})

			resid_h2 = sample_h2s_discrete_given_p_c(Y_tilde,h2_divisions,h2_priors_resids,tot_Y_prec,invert_aI_bZAZ)
			E_a_prec = tot_Y_prec / resid_h2
			resid_Y_prec = tot_Y_prec / (1-resid_h2)

		 # -----Sample B and E_a--------------- #
			#conditioning on W, F, Lambda
			Y_tilde = Y - F %*% t(Lambda) - Z_2 %*% W
			# location_sample = sample_means( Y_tilde, resid_Y_prec, E_a_prec, invert_aPXA_bDesignDesignT )
			# location_sample = sample_means_c( Y_tilde, resid_Y_prec, E_a_prec, invert_aPXA_bDesignDesignT)
			location_sample = sample_means_parallel_c( Y_tilde, resid_Y_prec, E_a_prec, invert_aPXA_bDesignDesignT ,1)
			B   = location_sample[1:b,]
			E_a = location_sample[b+(1:r),]

		 # -----Sample W ---------------------- #
			#conditioning on B, E_a, F, Lambda
			if(ncol(Z_2) > 0) {
				Y_tilde = Y - X %*% B - Z_1 %*% E_a - F %*% t(Lambda)
				# location_sample = sample_means( Y_tilde, resid_Y_prec, W_prec, invert_aPXA_bDesignDesignT_rand2 )
				# location_sample = sample_means_c( Y_tilde, resid_Y_prec, W_prec, invert_aPXA_bDesignDesignT_rand2 )
				location_sample = sample_means_parallel_c( Y_tilde, resid_Y_prec, W_prec, invert_aPXA_bDesignDesignT_rand2,1)
				W = location_sample
			}
		 # -----Sample F----------------------- #
			#conditioning on B,F_a,E_a,W,Lambda, F_h2
			
		 # -----Sample F_a_prec and F_e_prec -------------------- #
			#conditioning on F, F_a
			U = invert_aI_bZAZ$U
			s = invert_aI_bZAZ$s
			UtF = t(U) %*% F
			n = nrow(F)
			tot_F_prec = sapply(1:k,function(j) {
				Sigma_sqrt = sqrt(s*F_h2[j] + (1-F_h2[j]))
				SiUtF_j = UtF[,j] / Sigma_sqrt
				rgamma(1,shape = tot_F_prec_shape + n/2, rate = tot_F_prec_rate+1/2*sum(SiUtF_j^2))
			})

			F_h2 = sample_h2s_discrete_given_p_c(F,h2_divisions,h2_priors_factors,tot_F_prec,invert_aI_bZAZ)
			F_a_prec = tot_F_prec / F_h2
			F_e_prec = tot_F_prec / (1-F_h2)
			
		 # -----Sample F_a--------------------- #
			#conditioning on F, F_h2
	    	F_a = sample_F_a_ipx_c(F,Z_1,F_a_prec,F_e_prec,invert_aZZt_Ainv);

		 # -----Sample F----------------------- #
			#conditioning on B,F_a,E_a,W,Lambda, F_h2
			Y_tilde = Y - X %*% B - Z_1 %*% E_a - Z_2 %*% W
			# F1 = sample_factors_scores( Y_tilde, Z_1,Lambda,resid_Y_prec,F_a,F_h2 )
			F = sample_factors_scores_ipx_c( Y_tilde, Z_1,Lambda,resid_Y_prec,F_a,F_e_prec )
					
		 # -----Sample W_prec------------------ #
			if(ncol(Z_2) > 0) {
				W_prec =  rgamma(p, W_prec_shape + r2/2,rate = W_prec_rate + 1/2*colSums((chol_A_2_inv %*% W)^2))
			}
		 
	
	})
	current_state = current_state[current_state_names]

	 # -- adapt number of factors to samples ---#
		current_state = update_k( current_state, priors, run_parameters, data_matrices)
		
	 # -- save sampled values (after thinning) -- #
		if( (i-burn) %% thin == 0 && i > burn) {
				
			sp_num = (i-burn)/thin
			
			Posterior = save_posterior_samples_ipx( sp_num,current_state, Posterior)
			
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