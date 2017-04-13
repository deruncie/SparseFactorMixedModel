sample_factor_model.fast_BSFG = function(BSFG_state,grainSize,...) {
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
		UtEta = as.matrix(Ut %*% Eta)

	 # -----Sample Lambda and B ------------------ #
		#conditioning on F, marginalizing over U_R

		Design = as.matrix(cbind(X,F))
		UtDesign = as.matrix(Ut %*% Design)
		rows = b + k
		prior_mean = matrix(0,rows,p)
		if(b > 0) {
		  prior_prec = rbind(prec_B,t(Plam))
		} else{ # b == 0
		  prior_prec = t(Plam)
		}
		if(is.null(cis_genotypes)){
		  coefs = sample_coefs_parallel_sparse_c( UtEta,UtDesign,resid_h2, tot_Eta_prec,s, prior_mean,prior_prec,grainSize)
		  if(b > 0){
		    B[] = coefs[1:b,,drop=FALSE]
		  }
		  Lambda[] = t(coefs[b + 1:k,,drop=FALSE])
		  XB = X %*% B
		} else{
		  XB = matrix(0,ncol = p, nrow = n)
		  for(j in 1:p){
		    cis_X_j = cis_genotypes[[j]]
		    UtDesign_j = Ut %*% cis_X_j
		    if(is(UtDesign_j,'Matrix')) UtDesign_j = UtDesign_j@x
		    prior_mean_j = rbind(prior_mean[,j,drop=FALSE],0)
		    prior_prec_j = rbind(prior_prec[,j,drop=FALSE],1e-10)
		    coefs_j = sample_coefs_parallel_sparse_c(UtEta[,j,drop=FALSE],cbind(UtDesign,UtDesign_j),resid_h2[,j,drop=FALSE], tot_Eta_prec[,j,drop=FALSE],s,prior_mean_j,prior_prec_j,grainSize)
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
		Eta_tilde = Eta - XB - F %*% t(Lambda)
		UtEta_tilde = as.matrix(Ut %*% Eta_tilde)
		tot_Eta_prec[] = sample_tot_prec_sparse_c(UtEta_tilde,resid_h2,s,tot_Eta_prec_shape,tot_Eta_prec_rate)

		resid_h2_index = sample_h2s_discrete_given_p_sparse_c(UtEta_tilde,h2_divisions,h2_priors_resids,tot_Eta_prec,s)
		resid_h2[] = h2s_matrix[,resid_h2_index,drop=FALSE]

		U_R[] = sample_randomEffects_parallel_sparse_c(Eta_tilde, Z, tot_Eta_prec, resid_h2, invert_aZZt_Kinv, grainSize)

		resid_Eta_prec = tot_Eta_prec / (1-resid_h2)

	# -----Sample Lambda and B_F ------------------ #
		# marginalizing over random effects (conditional on F, F_h2, tot_F_prec, prec_B)
		if(b_F > 0){
		  prior_mean = matrix(0,b_F,k)
		  prior_prec = prec_B_F
		  B_F = sample_coefs_parallel_sparse_c(as.matrix(Ut %*% F),as.matrix(Ut %*% X_F),F_h2, tot_F_prec,s, prior_mean,prior_prec,grainSize)
		  XFBF = X_F %*% B_F
		  F_tilde = F - XFBF # not sparse.
		} else{
		  F_tilde = F
		}

	 # -----Sample F_h2 and tot_F_prec, U_F -------------------- #
		#conditioning on F, U_F
		UtF_tilde = as.matrix(Ut %*% F_tilde)
		tot_F_prec[] = sample_tot_prec_sparse_c(UtF_tilde,F_h2,s,tot_F_prec_shape,tot_F_prec_rate)

		F_h2_index = sample_h2s_discrete_given_p_sparse_c(UtF_tilde,h2_divisions,h2_priors_factors,tot_F_prec,s)
		F_h2[] = h2s_matrix[,F_h2_index,drop=FALSE]

    U_F[] = sample_randomEffects_parallel_sparse_c(F_tilde,Z,tot_F_prec, F_h2, invert_aZZt_Kinv, grainSize)

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
