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
		  prior_prec = rbind(prec_B,t(Plam))
		} else{ # b == 0
		  prior_prec = t(Plam)
		}
		if(is.null(cis_genotypes)){
		  coefs = sample_MME_fixedEffects(Eta,Design,Sigma_Choleskys, Sigma_Perm,  resid_h2_index, tot_Eta_prec, prior_mean, prior_prec,ncores)
		  if(b > 0){
		    B[] = coefs[1:b,,drop=FALSE]
		  }
		  Lambda[] = t(coefs[b + 1:k,,drop=FALSE])
		  XB = X %*% B
		} else{
		  XB = matrix(0,ncol = p, nrow = n)
  		for(j in 1:p){
  		  cis_X_j = cis_genotypes[[j]]
  		  Design_j = cbind(Design,cis_X_j)
  		  prior_mean_j = rbind(prior_mean[,j,drop=FALSE],0)
  		  prior_prec_j = rbind(prior_prec[,j,drop=FALSE],1e-10)
  		  coefs_j = sample_MME_fixedEffects(Eta[,j,drop=FALSE],Design_j,Sigma_Choleskys, Sigma_Perm,  resid_h2_index[j], tot_Eta_prec[,j,drop=FALSE], prior_mean_j, prior_prec_j,ncores)
  		  if(b > 0){
  		    B[,j] = coefs_j[1:b]
  		  }
  		  Lambda[j,] = coefs_j[b+1:k]
  		  cis_effects[,cis_effects_index[j]] = coefs_j[-c(1:(b+k))]
  		  XB[,j] = X %*% B[,j] + cis_X_j %*% cis_effects[cis_effects_index[j]]
  		}
		}


	 # -----Sample tot_Eta_prec, resid_h2, U_R ---------------- #
		#conditioning on B, F, Lambda, resid_h2, tot_Eta_prec
		Eta_tilde = Eta - XB - F %*% t(Lambda)
		tot_Eta_prec[] = sample_tot_prec(Eta_tilde, tot_Eta_prec_shape, tot_Eta_prec_rate, Sigma_Choleskys, Sigma_Perm, resid_h2_index,ncores)

		# big_Sigma_Choleskys = rep(Sigma_Choleskys,5)
		# big_Resid_discrete_priors = rep(Resid_discrete_priors,5)
		# big_candidate_states_h2 = rep(candidate_states_h2,5)
		# recover()
		# resid_h2_index = sample_h2s_discrete(Eta_tilde,tot_Eta_prec, Sigma_Choleskys, Sigma_Perm, Resid_discrete_priors,ncores)
		resid_h2_index = sample_h2s_discrete_MH(Eta_tilde,tot_Eta_prec, Sigma_Choleskys,Resid_discrete_priors,h2s_matrix,resid_h2_index,step_size = 0.2,ncores)
		# resid_h2_index2b = sample_h2s_discrete_MH(Eta_tilde,tot_Eta_prec, big_Sigma_Choleskys,big_Resid_discrete_priors,h2s_matrix,resid_h2_index,step_size = 0.2,ncores)
		# resid_h2_index3 = sample_h2s_discrete_MH2(Eta_tilde,tot_Eta_prec, Sigma_Choleskys,Resid_discrete_priors,h2s_matrix,resid_h2_index,candidate_states_h2,ncores)
		#
		# microbenchmark(
		#   sample_h2s_discrete(Eta_tilde,tot_Eta_prec, Sigma_Choleskys, Sigma_Perm, Resid_discrete_priors,ncores),
		#   sample_h2s_discrete_MH(Eta_tilde,tot_Eta_prec, Sigma_Choleskys,Resid_discrete_priors,h2s_matrix,resid_h2_index,step_size = 0.2,ncores),
		#   sample_h2s_discrete_MH(Eta_tilde,tot_Eta_prec, big_Sigma_Choleskys,big_Resid_discrete_priors,h2s_matrix,resid_h2_index,step_size = 0.2,ncores),
		#   sample_h2s_discrete_MH2(Eta_tilde,tot_Eta_prec, Sigma_Choleskys,Resid_discrete_priors,h2s_matrix,resid_h2_index,candidate_states_h2,ncores),
		#   sample_h2s_discrete_MH2(Eta_tilde,tot_Eta_prec, big_Sigma_Choleskys,big_Resid_discrete_priors,h2s_matrix,resid_h2_index,big_candidate_states_h2,ncores),
		#   times=10
		# )

		resid_h2[] = h2s_matrix[,resid_h2_index,drop=FALSE]

		U_R[] = sample_MME_ZKZts(Eta_tilde, Z, tot_Eta_prec, randomEffect_C_Choleskys, resid_h2, resid_h2_index,chol_Ki_mats,ncores)

		# -----Sample Lambda and B_F ------------------ #
		# F, marginalizing over random effects (conditional on F_h2, tot_F_prec)
		if(b_F > 0){
			prior_mean = matrix(0,b_F,p)
			prior_prec = prec_B_F
			B_F = sample_MME_fixedEffects(F,X_F,Sigma_Choleskys, Sigma_Perm, F_h2_index, tot_F_prec, prior_mean, prior_prec,ncores)
			XFBF = X_F %*% B_F
			F_tilde = F - XFBF
		} else{
		  F_tilde = F
		}
	 # -----Sample tot_F_prec, F_h2, U_F ---------------- #
		#conditioning on B, F, Lambda, F_h2, tot_F_prec

		tot_F_prec[] = sample_tot_prec(F_tilde, tot_F_prec_shape, tot_F_prec_rate, Sigma_Choleskys, Sigma_Perm, F_h2_index,ncores)
		F_h2_index = sample_h2s_discrete(F_tilde,tot_F_prec, Sigma_Choleskys, Sigma_Perm, F_discrete_priors,ncores)
		F_h2[] = h2s_matrix[,F_h2_index,drop=FALSE]

		U_F[] = sample_MME_ZKZts(F_tilde, Z, tot_F_prec, randomEffect_C_Choleskys, F_h2, F_h2_index,chol_Ki_mats,ncores)

	 # -----Sample F----------------------- #
		#conditioning on B, U_F,U_R,Lambda, F_h2
		Eta_tilde = as.matrix(Eta - XB - Z %*% U_R)
		F_e_prec = tot_F_prec / (1-colSums(F_h2))
		resid_Eta_prec = tot_Eta_prec / (1-colSums(resid_h2))
		if(b_F > 0) {
		  prior_mean = XFBF + Z %*% U_F
		} else {
		  prior_mean = Z %*% U_F
		}
		F[] = sample_factors_scores_sparse_c( Eta_tilde, as.matrix(prior_mean),Lambda,resid_Eta_prec,F_e_prec )

  }))
	current_state = current_state[current_state_names]

	return(current_state)
}

