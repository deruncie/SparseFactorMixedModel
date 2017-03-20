sample_factor_model.fast_BSFG = function(BSFG_state,...) {
  data_matrices  = BSFG_state$data_matrices
  priors         = BSFG_state$priors
  run_parameters = BSFG_state$run_parameters
  run_variables  = BSFG_state$run_variables
  current_state  = BSFG_state$current_state

	current_state_names = names(current_state)
	current_state = with(c(priors,run_parameters, run_variables,data_matrices),within(current_state, {
		k = ncol(Lambda)

	 # -----Sample Lambda and B ------------------ #
		#conditioning on F, marginalizing over E_a

		Design = as.matrix(cbind(X,F))
		rows = b + k
		prior_mean = matrix(0,rows,p)
		if(b > 0) {
		  prior_prec = rbind(matrix(prec_B,b,p),t(Plam))
		} else{ # b == 0
		  prior_prec = t(Plam)
		}
		coefs = sample_coefs_parallel_sparse_c( Eta,Design,resid_h2, tot_Eta_prec,prior_mean,prior_prec,invert_aI_bZKZ,1)

		if(b > 0){
			B[] = coefs[1:b,,drop=FALSE]
		}
		Lambda[] = t(coefs[b + 1:k,,drop=FALSE])


	 # -----Sample resid_h2, tot_Eta_prec, E_a ---------------- #
		#conditioning on W, B, F, Lambda, marginalizing over E_a
		Eta_tilde = as.matrix(Eta - X %*% B - F %*% t(Lambda))
		tot_Eta_prec[] = sample_tot_prec_sparse_c(Eta_tilde,resid_h2,tot_Eta_prec_shape,tot_Eta_prec_rate,invert_aI_bZKZ)

		resid_h2_index = sample_h2s_discrete_given_p_sparse_c(Eta_tilde,h2_divisions,h2_priors_resids,tot_Eta_prec,invert_aI_bZKZ)
		resid_h2[] = h2s_matrix[,resid_h2_index,drop=FALSE]

		E_a[] = sample_randomEffects_parallel_sparse_c( Eta_tilde, Z, tot_Eta_prec, resid_h2, invert_aZZt_Kinv, 1)

		resid_Eta_prec = tot_Eta_prec / (1-resid_h2)

	# -----Sample Lambda and B_F ------------------ #
		# marginalizing over random effects (conditional on F, F_h2, tot_F_prec, prec_B)
		if(b_F > 0){
		  prior_mean = matrix(0,b_F,p)
		  prior_prec = matrix(prec_B[-1],b_F,k)
		  B_F = sample_coefs_parallel_sparse_c(F,X_F,F_h2, tot_F_prec,prior_mean,prior_prec,invert_aI_bZKZ,1)
		  F_tilde = F - X_F %*% B_F # not sparse.
		} else{
		  F_tilde = F
		}

	 # -----Sample F_h2 and tot_F_prec, F_a -------------------- #
		#conditioning on F, F_a
		tot_F_prec[] = sample_tot_prec_sparse_c(F_tilde,F_h2,tot_F_prec_shape,tot_F_prec_rate,invert_aI_bZKZ)

		F_h2_index = sample_h2s_discrete_given_p_sparse_c(F_tilde,h2_divisions,h2_priors_factors,tot_F_prec,invert_aI_bZKZ)
		F_h2[] = h2s_matrix[,F_h2_index,drop=FALSE]

    F_a[] = sample_randomEffects_parallel_sparse_c(F_tilde,Z,tot_F_prec, F_h2, invert_aZZt_Kinv, 1)

	 # -----Sample F----------------------- #
		#conditioning on B, F_a,E_a,W,Lambda, F_h2
		Eta_tilde = as.matrix(Eta - X %*% B - Z %*% E_a)
		F_e_prec = tot_F_prec / (1-F_h2)
		if(b_F > 0) {
		  prior_mean = as.matrix(X_F %*% B_F + Z %*% F_a)
		} else {
		  prior_mean = as.matrix(Z %*% F_a)
		}
		F[] = sample_factors_scores_sparse_c( Eta_tilde, prior_mean,Lambda,resid_Eta_prec,F_e_prec )

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
