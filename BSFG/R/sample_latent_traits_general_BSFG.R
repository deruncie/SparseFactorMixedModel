sample_latent_traits.general_BSFG = function(BSFG_state,grainSize = 1,...) {
  data_matrices  = BSFG_state$data_matrices
  priors         = BSFG_state$priors
  run_parameters = BSFG_state$run_parameters
  run_variables  = BSFG_state$run_variables
  current_state  = BSFG_state$current_state

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

    # -----Sample tot_Eta_prec, resid_h2, U_R ---------------- #
    #conditioning on B, F, Lambda, resid_h2, tot_Eta_prec
    Eta_tilde = Eta - XB - F %*% t(Lambda)
    tot_Eta_prec[] = sample_tot_prec(Eta_tilde, tot_Eta_prec_shape, tot_Eta_prec_rate, Sigma_Choleskys, resid_h2_index,grainSize)

    if(is.null(h2_step_size)) {
      resid_h2_index = sample_h2s_discrete(Eta_tilde,tot_Eta_prec, Sigma_Choleskys, h2_priors_resids,grainSize)
    } else{
      resid_h2_index = sample_h2s_discrete_MH(Eta_tilde,tot_Eta_prec, Sigma_Choleskys,h2_priors_resids,h2s_matrix,resid_h2_index,h2_step_size,grainSize)
    }

    resid_h2[] = h2s_matrix[,resid_h2_index,drop=FALSE]

    U_R[] = sample_MME_ZKZts(Eta_tilde, Z, tot_Eta_prec, randomEffect_C_Choleskys, resid_h2, resid_h2_index,grainSize)

    resid_Eta_prec = tot_Eta_prec / (1-colSums(resid_h2))

    # -----Sample Lambda and B_F ------------------ #
    # F, marginalizing over random effects (conditional on F_h2, tot_F_prec)
    if(b_F > 0){
      prior_mean = matrix(0,b_F,k)
      prior_prec = prec_B_F
      B_F = sample_MME_fixedEffects(F,X_F,Sigma_Choleskys, F_h2_index, tot_F_prec, prior_mean, prior_prec,grainSize)
      XFBF = X_F %*% B_F
      if(class(XFBF) == 'Matrix') XFBF = as.matrix(XFBF)
      F_tilde = F - XFBF
    } else{
      F_tilde = F
    }
    # -----Sample tot_F_prec, F_h2, U_F ---------------- #
    #conditioning on B, F, Lambda, F_h2, tot_F_prec

    tot_F_prec[] = sample_tot_prec(F_tilde, tot_F_prec_shape, tot_F_prec_rate, Sigma_Choleskys, F_h2_index,grainSize)

    if(is.null(h2_step_size)) {
      F_h2_index = sample_h2s_discrete(F_tilde,tot_F_prec, Sigma_Choleskys,h2_priors_factors,grainSize)
    } else{
      F_h2_index = sample_h2s_discrete_MH(F_tilde,tot_F_prec, Sigma_Choleskys,h2_priors_factors,h2s_matrix,F_h2_index,h2_step_size,grainSize)
    }

    F_h2[] = h2s_matrix[,F_h2_index,drop=FALSE]

    U_F[] = sample_MME_ZKZts(F_tilde, Z, tot_F_prec, randomEffect_C_Choleskys, F_h2, F_h2_index,grainSize)

    # -----Sample F----------------------- #
    #conditioning on B, U_F,U_R,Lambda, F_h2
    Eta_tilde = Eta - XB - as.matrix(Z %*% U_R)
    F_e_prec = tot_F_prec / (1-colSums(F_h2))
    prior_mean = as.matrix(Z %*% U_F)
    if(b_F > 0) {
      prior_mean = prior_mean + XFBF
    }
    randn = matrix(rnorm(n*k),n)
    F[] = sample_factors_scores_sparse_c_Eigen( Eta_tilde, prior_mean,Lambda,resid_Eta_prec,F_e_prec,randn )

  }))
  current_state = current_state[current_state_names]

  return(current_state)
}

