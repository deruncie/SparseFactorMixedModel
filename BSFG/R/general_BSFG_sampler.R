sample_factor_model.general_BSFG = function(BSFG_state,ncores = detectCores(),...) {
	data_matrices  = BSFG_state$data_matrices
	priors         = BSFG_state$priors
	run_parameters = BSFG_state$run_parameters
	run_variables  = BSFG_state$run_variables
	current_state  = BSFG_state$current_state

	current_state_names = names(current_state)
	current_state = with(c(priors,run_parameters, run_variables,data_matrices),within(current_state, {
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


	 # -----Sample tot_Eta_prec, resid_h2, E_a ---------------- #
		#conditioning on B, F, Lambda, resid_h2, tot_Eta_prec
		Eta_tilde = Eta - X %*% B - F %*% t(Lambda)
		tot_Eta_prec[] = sample_tot_prec(Eta_tilde, tot_Eta_prec_shape, tot_Eta_prec_rate, Sigma_Choleskys, Sigma_Perm, resid_h2_index,ncores)
		resid_h2_index = sample_h2s_discrete(Eta_tilde,tot_Eta_prec, Sigma_Choleskys, Sigma_Perm, Resid_discrete_priors,ncores)
		resid_h2[] = h2s_matrix[,resid_h2_index,drop=FALSE]
		# E_a_prec = tot_Eta_prec / colSums(resid_h2)

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
  }))
	current_state = current_state[current_state_names]

	return(current_state)
}

