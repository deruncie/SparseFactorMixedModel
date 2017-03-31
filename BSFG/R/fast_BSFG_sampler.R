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
		#conditioning on F, marginalizing over U_R

		Design = as.matrix(cbind(X,F))
		rows = b + k
		prior_mean = matrix(0,rows,p)
		if(b > 0) {
		  prior_prec = rbind(prec_B,t(Plam))
		} else{ # b == 0
		  prior_prec = t(Plam)
		}
		if(is.null(cis_genotypes)){
		  coefs = sample_coefs_parallel_sparse_c( Eta,Design,resid_h2, tot_Eta_prec,prior_mean,prior_prec,invert_aI_bZKZ,1)
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
		    coefs_j = sample_coefs_parallel_sparse_c(Eta[,j,drop=FALSE],Design_j,resid_h2[,j,drop=FALSE], tot_Eta_prec[,j,drop=FALSE],prior_mean_j,prior_prec_j,invert_aI_bZKZ,1)
		    # j2 = rep(j,100)
		    # prior_mean_j2 = rbind(prior_mean[,j2,drop=FALSE],0)
		    # prior_prec_j2 = rbind(prior_prec[,j2,drop=FALSE],1e-10)
		    # coefs_j = sample_coefs_parallel_sparse_c(Eta[,j2,drop=FALSE],Design_j,resid_h2[,j2,drop=FALSE], tot_Eta_prec[,j2,drop=FALSE],prior_mean_j2,prior_prec_j2,invert_aI_bZKZ,1)
		    # coefs_j2 = sample_coefs_parallel_sparse_c2(Eta[,j2,drop=FALSE],Design_j,resid_h2[,j2,drop=FALSE], tot_Eta_prec[,j2,drop=FALSE],prior_mean_j2,prior_prec_j2,invert_aI_bZKZ,1)
		    # recover()
		    if(b > 0){
		      B[,j] = coefs_j[1:b]
		    }
		    Lambda[j,] = coefs_j[b+1:k]
		    cis_effects[,cis_effects_index[j]] = coefs_j[-c(1:(b+k))]
		    XB[,j] = X %*% B[,j] + cis_X_j %*% cis_effects[cis_effects_index[j]]
		  }
		}


	 # -----Sample resid_h2, tot_Eta_prec, U_R ---------------- #
		#conditioning on W, B, F, Lambda, marginalizing over U_R
		Eta_tilde = as.matrix(Eta - XB - F %*% t(Lambda))
		tot_Eta_prec[] = sample_tot_prec_sparse_c(Eta_tilde,resid_h2,tot_Eta_prec_shape,tot_Eta_prec_rate,invert_aI_bZKZ)

		resid_h2_index = sample_h2s_discrete_given_p_sparse_c(Eta_tilde,h2_divisions,h2_priors_resids,tot_Eta_prec,invert_aI_bZKZ)
		resid_h2[] = h2s_matrix[,resid_h2_index,drop=FALSE]

		U_R[] = sample_randomEffects_parallel_sparse_c( Eta_tilde, Z, tot_Eta_prec, resid_h2, invert_aZZt_Kinv, 1)

		resid_Eta_prec = tot_Eta_prec / (1-resid_h2)

	# -----Sample Lambda and B_F ------------------ #
		# marginalizing over random effects (conditional on F, F_h2, tot_F_prec, prec_B)
		if(b_F > 0){
		  prior_mean = matrix(0,b_F,p)
		  prior_prec = prec_B_F
		  B_F = sample_coefs_parallel_sparse_c(F,X_F,F_h2, tot_F_prec,prior_mean,prior_prec,invert_aI_bZKZ,1)
		  XFBF = X_F %*% B_F
		  F_tilde = F - XFBF # not sparse.
		} else{
		  F_tilde = F
		}

	 # -----Sample F_h2 and tot_F_prec, U_F -------------------- #
		#conditioning on F, U_F
		tot_F_prec[] = sample_tot_prec_sparse_c(F_tilde,F_h2,tot_F_prec_shape,tot_F_prec_rate,invert_aI_bZKZ)

		F_h2_index = sample_h2s_discrete_given_p_sparse_c(F_tilde,h2_divisions,h2_priors_factors,tot_F_prec,invert_aI_bZKZ)
		F_h2[] = h2s_matrix[,F_h2_index,drop=FALSE]

    U_F[] = sample_randomEffects_parallel_sparse_c(F_tilde,Z,tot_F_prec, F_h2, invert_aZZt_Kinv, 1)

	 # -----Sample F----------------------- #
		#conditioning on B, U_F,U_R,W,Lambda, F_h2
		Eta_tilde = as.matrix(Eta - XB - Z %*% U_R)
		F_e_prec = tot_F_prec / (1-F_h2)
		if(b_F > 0) {
		  prior_mean = as.matrix(XFBF + Z %*% U_F)
		} else {
		  prior_mean = as.matrix(Z %*% U_F)
		}
		F[] = sample_factors_scores_sparse_c( Eta_tilde, prior_mean,Lambda,resid_Eta_prec,F_e_prec )

  }))
	current_state = current_state[current_state_names]

	return(current_state)
}
