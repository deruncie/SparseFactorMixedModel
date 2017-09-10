sample_latent_traits.fast_missing_BSFG = function(BSFG_state,grainSize,...) {
  data_matrices  = BSFG_state$data_matrices
  priors         = BSFG_state$priors
  run_parameters = BSFG_state$run_parameters
  run_variables  = BSFG_state$run_variables
  current_state  = BSFG_state$current_state
  Y_row_obs_index = BSFG_state$run_variables$Y_row_obs_index
  Y_row_obs_sets = BSFG_state$run_variables$Y_row_obs_sets
  invert_aI_bZKZ = BSFG_state$run_variables$invert_aI_bZKZ
  Y_missing = BSFG_state$run_parameters$observation_model_parameters$Y_missing
  Y_col_obs_index = BSFG_state$run_variables$Y_col_obs_index
  Qt = BSFG_state$run_variables$invert_aI_bZKZ[[1]]$Qt
  s = BSFG_state$run_variables$invert_aI_bZKZ[[1]]$s

  current_state_names = names(current_state)
  current_state = with(c(priors,run_parameters, run_variables,data_matrices),within(current_state, {
    k = ncol(Lambda)

    XB = X %**% B
    if(!is.null(cis_genotypes)){
      for(j in 1:p){
        cis_X_j = cis_genotypes[[j]]
        XB[,j] = XB[,j] + cis_X_j %*% cis_effects[cis_effects_index[j]:(cis_effects_index[j+1]-1)]
      }
    }

    # -----Sample resid_h2, tot_Eta_prec, U_R ---------------- #
    #conditioning on W, B, F, Lambda, marginalizing over U_R
    Eta_tilde = Eta - XB - F %*% t(Lambda)
    scores = tot_prec_scores_missing_c(Eta_tilde,resid_h2,invert_aI_bZKZ)
    n_obs = sapply(Y_col_obs_index,function(x) length(invert_aI_bZKZ[[x]]$Y_obs))
    tot_Eta_prec[] = rgamma(p,shape = tot_Eta_prec_shape + n_obs/2,rate = tot_Eta_prec_rate + 0.5*scores)

    if(!length(h2_priors_resids) == ncol(h2s_matrix)) stop('wrong length of h2_priors_resids')
    if(is.null(h2_step_size)) {
      log_ps = log_p_h2s_fast_missing(Eta_tilde, tot_Eta_prec, h2_priors_resids,invert_aI_bZKZ,grainSize)
      resid_h2_index = sample_h2s(log_ps,grainSize)
    } else{
      resid_h2_index = sample_h2s_discrete_MH_fast_missing_c(Eta_tilde,tot_Eta_prec,h2_priors_resids,resid_h2_index,h2s_matrix,h2_step_size,invert_aI_bZKZ,grainSize)+1
    }
    resid_h2[] = h2s_matrix[,resid_h2_index,drop=FALSE]

    U_R[] = sample_randomEffects_parallel_sparse_missing_c_Eigen(Eta_tilde,tot_Eta_prec,resid_h2,invert_aZZt_Kinv,ncol(Z),grainSize)

    resid_Eta_prec = tot_Eta_prec / (1-resid_h2)

    # -----Sample Lambda and B_F ------------------ #
    # tot_F_prec[] = 1
    # marginalizing over random effects (conditional on F, F_h2, tot_F_prec, prec_B)

    if(b_F > 0){
      resid_F_prec = uncorrelated_prec_mat(F_h2,tot_F_prec,s)
      # non-QTL fixed effects
      X_F1 = X_F
      b_F1 = ncol(X_F1)
      F_tilde = F
      if(length(QTL_columns_factors) > 0) {
        X_F1 = X_F[,-QTL_columns_factors,drop=FALSE]
        b_F1 = ncol(X_F1)
        F_tilde = F - toDense(QTL_factors_Z %*% QTL_factors_X %*% B_F[-c(1:b_F1),,drop=FALSE])
      }
      prior_mean = matrix(0,b_F1,k)
      prior_prec = B_F_prec[1:b_F1,,drop=FALSE] * tot_F_prec[rep(1,b_F1),,drop=FALSE]  # prior for B_F includes tot_F_prec
      B_F[1:b_F1,] = sample_coefs_parallel_sparse_c_Eigen( Qt %**% F_tilde,Qt %**% X_F1,resid_F_prec, prior_mean,prior_prec,grainSize)

      # QTL fixed effects
      if(length(QTL_columns_factors) > 0){
        F_tilde = F - toDense(X_F1 %*% B_F[1:b_F1,,drop=FALSE])
        b_F_QTL = ncol(QTL_factors_X)
        prior_mean = matrix(0,b_F_QTL,k)
        prior_prec = B_F_prec[QTL_columns_factors,] * tot_F_prec[rep(1,b_F_QTL),]  # prior for B_F includes tot_F_prec
        randn_theta = matrix(rnorm(b_F_QTL*k),b_F_QTL)
        randn_e = matrix(rnorm(n*k),n)
        B_F[QTL_columns_factors,] = sample_coefs_hierarchical_parallel_sparse_c_Eigen(Qt,F_tilde,QTL_factors_Z,QTL_factors_X,F_h2, tot_F_prec,s, prior_mean,prior_prec,randn_theta,randn_e,grainSize)
      }

      XFBF = X_F %**% B_F
      F_tilde = F - XFBF # not sparse.
    } else{
      F_tilde = F
    }

    # -----Sample F_h2 and tot_F_prec, U_F -------------------- #
    #conditioning on F, U_F
    QtF_tilde = Qt %**% F_tilde

    resid_F_prec = uncorrelated_prec_mat(F_h2,rep(1,p),s)
    scores = tot_prec_scores_c(QtF_tilde,resid_F_prec)
    if(b_F > 0) {
      scores = scores + colSums((B_F^2*B_F_prec)[!X_F_zero_variance,,drop=FALSE])   # add this if tot_F_prec part of the prior for B_F
    }
    tot_F_prec[] = rgamma(k,shape = tot_F_prec_shape + n/2 + sum(!X_F_zero_variance)/2,rate = tot_F_prec_rate + scores/2)
    # tot_F_prec[] = rgamma(k,shape = tot_F_prec_shape + n/2,rate = tot_F_prec_rate + scores/2)

    if(!length(h2_priors_factors) == ncol(h2s_matrix)) stop('wrong length of h2_priors_factors')
    if(is.null(h2_step_size)) {
      log_ps = log_p_h2s_fast(QtF_tilde, tot_F_prec, h2_priors_factors,s,grainSize)
      F_h2_index = sample_h2s(log_ps,grainSize)
    } else{
      F_h2_index = sample_h2s_discrete_MH_fast_c(QtF_tilde,tot_F_prec,h2_priors_factors,F_h2_index,h2s_matrix,s,h2_step_size,grainSize)+1
    }
    F_h2[] = h2s_matrix[,F_h2_index,drop=FALSE]

    U_F[] = sample_randomEffects_parallel_sparse_c_Eigen(F_tilde, Z, tot_F_prec, F_h2, invert_aZZt_Kinv[[1]], grainSize)

    # -----Sample F----------------------- #
    #conditioning on B, U_F,U_R,W,Lambda, F_h2
    Eta_tilde = Eta - XB - Z %**% U_R
    F_e_prec = tot_F_prec / (1-F_h2)
    prior_mean = Z %**% U_F
    if(b_F > 0) {
      prior_mean = prior_mean + XFBF
    }
    F[] = sample_factors_scores_sparse_missing_c_Eigen( Eta_tilde, prior_mean,Lambda,resid_Eta_prec,F_e_prec,Y_row_obs_sets,grainSize)
  }))
  current_state = current_state[current_state_names]

  return(current_state)
}
