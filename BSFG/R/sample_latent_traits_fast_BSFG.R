sample_latent_traits.fast_BSFG = function(BSFG_state,grainSize,...) {
  data_matrices  = BSFG_state$data_matrices
  priors         = BSFG_state$priors
  run_parameters = BSFG_state$run_parameters
  run_variables  = BSFG_state$run_variables
  current_state  = BSFG_state$current_state
  invert_aI_bZKZ = BSFG_state$run_variables$invert_aI_bZKZ
  Ut = t(invert_aI_bZKZ$U)
  s = invert_aI_bZKZ$s

	current_state_names = names(current_state)
	current_state = with(c(priors,run_parameters, run_variables,data_matrices),within(current_state, {
		k = ncol(Lambda)

		if(is.null(cis_genotypes)){
		  XB = X %*% B
		} else{
		  XB = matrix(0,ncol = p, nrow = n)
		  for(j in 1:p){
		    cis_X_j = cis_genotypes[[j]]
		    XB[,j] = X %*% B[,j] + cis_X_j %*% cis_effects[cis_effects_index[j]]
		  }
		}
		if(class(XB) == 'Matrix') XB = as.matrix(XB)

	 # -----Sample resid_h2, tot_Eta_prec, U_R ---------------- #
		#conditioning on W, B, F, Lambda, marginalizing over U_R
		Eta_tilde = Eta - XB - F %*% t(Lambda)
		UtEta_tilde = as.matrix(Ut %*% Eta_tilde)
		# randg_draws = rgamma(p,shape = tot_Eta_prec_shape + n/2,rate = 1)
		# tot_Eta_prec[] = sample_tot_prec_sparse_c(UtEta_tilde,resid_h2,s,tot_Eta_prec_rate,randg_draws)
		scores = tot_prec_scores_c(UtEta_tilde,resid_h2,s)
		tot_Eta_prec[] = rgamma(p,shape = tot_Eta_prec_shape + n/2,rate = tot_Eta_prec_rate + 0.5*scores)

		# resid_h2_index = sample_h2s_discrete_given_p_sparse_c(UtEta_tilde,h2_divisions,h2_priors_resids,tot_Eta_prec,s)
		resid_h2_index = sample_h2s_discrete_fast(UtEta_tilde, tot_Eta_prec, h2_priors_resids,s,grainSize)
		resid_h2[] = h2s_matrix[,resid_h2_index,drop=FALSE]

		# U_R[] = sample_randomEffects_parallel_sparse_c(Eta_tilde, Z, tot_Eta_prec, resid_h2, invert_aZZt_Kinv, grainSize)
		randn = matrix(rnorm(ncol(Z)*p),ncol(Z))
		U_R[] = sample_randomEffects_parallel_sparse_c_Eigen(Eta_tilde, Z, tot_Eta_prec, resid_h2, invert_aZZt_Kinv, randn,grainSize)

		resid_Eta_prec = tot_Eta_prec / (1-resid_h2)

	# -----Sample Lambda and B_F ------------------ #
		# tot_F_prec[] = 1
		# marginalizing over random effects (conditional on F, F_h2, tot_F_prec, prec_B)
		if(b_F > 0){
		  prior_mean = matrix(0,b_F,k)
		  prior_prec = sweep(prec_B_F,2,tot_F_prec,'*')  # prior for B_F includes tot_F_prec
		  if(b_F > 100){
		    n_sets = ceiling(b_F/100)
		    sets = gl(n_sets,b_F/n_sets)
		    for(set in unique(sets)){
		      index = sets==set
		      if(sum(!index) > 0){
		        X_F_set = X_F[,index,drop=FALSE]
		        F_tilde = F - X_F[,!index,drop=FALSE] %*% B_F[!index,,drop=FALSE]
		      } else{
		        X_F_set = X_F
		        F_tilde = F
		      }
		      # B_F[index,] = sample_coefs_parallel_sparse_c(as.matrix(Ut %*% F_tilde),
		      #                                              as.matrix(Ut %*% X_F_set),
		      #                                              F_h2, tot_F_prec,s,
		      #                                              prior_mean[index,],
		      #                                              prior_prec[index,],
		      #                                              grainSize)
		      randn_theta = matrix(rnorm(ncol(X_F_set)*k),ncol(X_F_set))
		      randn_e = matrix(rnorm(n*k),n)
		      # recover()
		      B_F[index,] = sample_coefs_parallel_sparse_c_Eigen(Ut,F_tilde,X_F_set,
		                                                  F_h2, tot_F_prec,s,
		                                                  prior_mean[index,],
		                                                  prior_prec[index,],
		                                                  randn_theta,randn_e,
		                                                  grainSize)
		    }
		  } else{
		    # B_F = sample_coefs_parallel_sparse_c(as.matrix(Ut %*% F),as.matrix(Ut %*% X_F),F_h2, tot_F_prec,s, prior_mean,prior_prec,grainSize)
		    randn_theta = matrix(rnorm(b_F*k),b_F)
		    randn_e = matrix(rnorm(n*k),n)
		    B_F = sample_coefs_parallel_sparse_c_Eigen(Ut,F,X_F,F_h2, tot_F_prec,s, prior_mean,prior_prec,randn_theta,randn_e,grainSize)
		  }
		  XFBF = X_F %*% B_F
		  if(class(XFBF) == 'Matrix') XFBF = as.matrix(XFBF)
		  F_tilde = F - XFBF # not sparse.
		} else{
		  F_tilde = F
		}

	 # -----Sample F_h2 and tot_F_prec, U_F -------------------- #
		#conditioning on F, U_F
		UtF_tilde = as.matrix(Ut %*% F_tilde)

		if(nrun > 0) {
  		if(b_F == 0) {
    		# randg_draws = rgamma(k,shape = tot_F_prec_shape + n/2,rate = 1)
    		# tot_F_prec[] = sample_tot_prec_sparse_c(UtF_tilde,F_h2,s,tot_F_prec_rate,randg_draws)
    		scores = tot_prec_scores_c(UtF_tilde,F_h2,s)
    		tot_F_prec[] = rgamma(k,shape = tot_F_prec_shape + n/2,rate = tot_F_prec_rate + 0.5*scores)
  		} else{
  		  # randg_draws = rgamma(k,shape = tot_F_prec_shape + n/2 + b_F/2,rate = 1)
  		  # tot_F_prec[] = sample_tot_prec_sparse_withX_c(UtF_tilde,B_F,F_h2,s,prec_B_F,tot_F_prec_rate,randg_draws)
  		  scores = tot_prec_scores_withX_c(UtF_tilde,B_F,F_h2,s,prec_B_F)
  		  tot_F_prec[] = rgamma(k,shape = tot_F_prec_shape + n/2+ b_F/2,rate = tot_F_prec_rate + 0.5*scores)
  		}
		} else{
		  tot_F_prec[] = 1
		}

		# F_h2_index = sample_h2s_discrete_given_p_sparse_c(UtF_tilde,h2_divisions,h2_priors_factors,tot_F_prec,s)
		F_h2_index = sample_h2s_discrete_fast(UtF_tilde, tot_F_prec, h2_priors_factors,s,grainSize)
		F_h2[] = h2s_matrix[,F_h2_index,drop=FALSE]

    # U_F[] = sample_randomEffects_parallel_sparse_c(F_tilde,Z,tot_F_prec, F_h2, invert_aZZt_Kinv, grainSize)
    randn = matrix(rnorm(ncol(Z)*k),ncol(Z))
    U_F[] = sample_randomEffects_parallel_sparse_c_Eigen(F_tilde, Z, tot_F_prec, F_h2, invert_aZZt_Kinv, randn,grainSize)

	 # -----Sample F----------------------- #
		#conditioning on B, U_F,U_R,W,Lambda, F_h2
		Eta_tilde = Eta - XB - as.matrix(Z %*% U_R)
		F_e_prec = tot_F_prec / (1-F_h2)
		prior_mean = as.matrix(Z %*% U_F)
		if(b_F > 0) {
		  prior_mean = prior_mean + XFBF
		}
		# F[] = sample_factors_scores_sparse_c( Eta_tilde, prior_mean,Lambda,resid_Eta_prec,F_e_prec )
		randn = matrix(rnorm(n*k),n)
		F[] = sample_factors_scores_sparse_c_Eigen( Eta_tilde, prior_mean,Lambda,resid_Eta_prec,F_e_prec,randn )
		# tot_F_prec[] = temp
  }))
	current_state = current_state[current_state_names]

	return(current_state)
}
