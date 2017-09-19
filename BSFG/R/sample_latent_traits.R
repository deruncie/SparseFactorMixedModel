sample_latent_traits = function(BSFG_state,grainSize = 1,...) {
  data_matrices  = BSFG_state$data_matrices
  priors         = BSFG_state$priors
  run_parameters = BSFG_state$run_parameters
  run_variables  = BSFG_state$run_variables
  current_state  = BSFG_state$current_state

  current_state_names = names(current_state)
  current_state = with(c(priors,run_parameters, run_variables,data_matrices),within(current_state, {
    k = ncol(Lambda)

    XB = X %*% B
    if(inherits(XB,'Matrix')) XB = as.matrix(XB)
    if(!is.null(cis_genotypes)){
      for(j in 1:p){
        cis_X_j = cis_genotypes[[j]]
        XB[,j] = XB[,j] + cis_X_j %*% cis_effects[cis_effects_index[j]:(cis_effects_index[j+1]-1)]
      }
    }

    # -----Sample tot_Eta_prec, resid_h2, U_R ---------------- #
    #conditioning on B, F, Lambda, resid_h2, tot_Eta_prec
    Eta_tilde = Eta - XB - F %*% t(Lambda)
    scores = rep(0,p)

    # sample columns of tot_Eta_prec, resid_h2, U_R in sets with the same patterns of missing data
    for(set in seq_along(Missing_data_map)){
      cols = Missing_data_map[[set]]$Y_cols
      rows = Missing_data_map[[set]]$Y_obs
      if(length(cols) == 0 || length(rows) == 0) next

      QtEta_tilde_set = Qt_list[[set]] %**% Eta_tilde[rows,cols]
      scores[cols] = tot_prec_scores(QtEta_tilde_set,
                                     Sigma_Choleskys_list[[set]],
                                     resid_h2_index[cols],
                                     grainSize
                                    )
      tot_Eta_prec[cols] = rgamma(length(cols),shape = tot_Eta_prec_shape + length(rows)/2, rate = tot_Eta_prec_rate[cols] + 0.5*scores[cols])

      if(!length(h2_priors_resids) == ncol(h2s_matrix)) stop('wrong length of h2_priors_resids')
      if(is.null(h2_step_size)) {
        log_ps = log_p_h2s(QtEta_tilde_set,
                           tot_Eta_prec[cols],
                           Sigma_Choleskys_list[[set]],
                           h2_priors_resids,
                           grainSize)
        resid_h2_index[cols] = sample_h2s(log_ps,grainSize)
      } else{
        resid_h2_index[cols] = sample_h2s_discrete_MH_c(QtEta_tilde_set,
                                                      tot_Eta_prec[cols],
                                                      h2_priors_resids,
                                                      resid_h2_index[cols],
                                                      h2s_matrix,
                                                      Sigma_Choleskys_list[[set]],
                                                      h2_step_size,
                                                      grainSize)
      }
      resid_h2[,cols] = h2s_matrix[,resid_h2_index[cols],drop=FALSE]

      U_R[,cols] = sample_MME_ZKZts_c(Eta_tilde[rows,cols],
                                    Z[rows,],
                                    tot_Eta_prec[cols],
                                    randomEffect_C_Choleskys_list[[set]],
                                    resid_h2[,cols,drop=FALSE],
                                    resid_h2_index[cols],
                                    grainSize)
    }

    resid_Eta_prec = tot_Eta_prec / (1-colSums(resid_h2))

    # -----Sample Lambda and B_F ------------------ #
    # F, marginalizing over random effects (conditional on F_h2, tot_F_prec)
    if(b_F > 0){
      # non-QTL fixed effects
      b_F1 = ncol(X_F)
      if(length(QTL_columns_factors) > 0) {
        X_F1_cols = seq_len(b_F1)[-QTL_columns_factors]
        b_F1 = length(X_F1_cols)
        F_tilde = F - toDense(QTL_factors_Z %*% QTL_factors_X %*% B_F[-c(1:b_F1),,drop=FALSE])
      } else{
        X_F1_cols = seq_len(b_F1)
        F_tilde = F
      }
      prior_mean = matrix(0,b_F1,k)
      prior_prec = B_F_prec[X_F1_cols,,drop=FALSE] * tot_F_prec[rep(1,b_F1),,drop=FALSE]  # prior for B_F includes tot_F_prec
      B_F[1:b_F1,] = sample_MME_fixedEffects_c(Qt_list[[1]] %**% F_tilde,
                                             QtXF_list[[1]][,X_F1_cols,drop=FALSE],
                                             Sigma_Choleskys_list[[1]],
                                             F_h2_index,
                                             tot_F_prec,
                                             prior_mean,
                                             prior_prec,
                                             grainSize)
      # QTL fixed effects
      if(length(QTL_columns_factors) > 0){
        warning('QTL_columns_factors not yet implemented')
        # F_tilde = F - toDense(X_F[,X_F1_cols] %*% B_F[1:b_F1,,drop=FALSE])
        # b_F_QTL = ncol(QTL_factors_X)
        # prior_mean = matrix(0,b_F_QTL,k)
        # prior_prec = B_F_prec[QTL_columns_factors,] * tot_F_prec[rep(1,b_F_QTL),]  # prior for B_F includes tot_F_prec
        # B_F[QTL_columns_factors,] = sample_coefs_hierarchical_parallel_sparse_c_Eigen(Qt,F_tilde,QTL_factors_Z,QTL_factors_X,F_h2, tot_F_prec,s, prior_mean,prior_prec,grainSize)
      }

      XFBF = X_F %**% B_F
      F_tilde = F - XFBF # not sparse.
    } else{
      F_tilde = F
    }

    # -----Sample tot_F_prec, F_h2, U_F ---------------- #
    #conditioning on B, F, Lambda, F_h2, tot_F_prec

    Qt_F_tilde = Qt_list[[1]] %**% F_tilde

    scores = tot_prec_scores(Qt_F_tilde,
                             Sigma_Choleskys_list[[1]],
                             F_h2_index,
                             grainSize)
    if(b_F > 0) {
      scores = scores + colSums((B_F^2*B_F_prec)[!X_F_zero_variance,,drop=FALSE])   # add this if tot_F_prec part of the prior for B_F
    }
    tot_F_prec[] = rgamma(k,shape = tot_F_prec_shape + n/2+ sum(!X_F_zero_variance)/2, rate = tot_F_prec_rate + 0.5*scores)
    # tot_F_prec[] = 1

    if(!length(h2_priors_factors) == ncol(h2s_matrix)) stop('wrong length of h2_priors_factors')
    if(is.null(h2_step_size)) {
      log_ps = log_p_h2s(Qt_F_tilde,
                         tot_F_prec,
                         Sigma_Choleskys_list[[1]],
                         h2_priors_factors,
                         grainSize)
      F_h2_index = sample_h2s(log_ps,grainSize)
    } else{
      F_h2_index = sample_h2s_discrete_MH_c(Qt_F_tilde,
                                          tot_F_prec,
                                          h2_priors_factors,
                                          F_h2_index,
                                          h2s_matrix,
                                          Sigma_Choleskys_list[[1]],
                                          h2_step_size,
                                          grainSize)
    }
    F_h2[] = h2s_matrix[,F_h2_index,drop=FALSE]

    U_F[] = sample_MME_ZKZts_c(F_tilde, Z, tot_F_prec, randomEffect_C_Choleskys_list[[1]], F_h2, F_h2_index,grainSize)

    # -----Sample F----------------------- #
    #conditioning on B, U_F,U_R,Lambda, F_h2
    Eta_tilde = Eta - XB - Z %**% U_R
    F_e_prec = tot_F_prec / (1-colSums(F_h2))
    prior_mean = Z %**% U_F
    if(b_F > 0) {
      prior_mean = prior_mean + XFBF
    }

    for(set in seq_along(Missing_data_map)){
      cols = Missing_data_map[[set]]$Y_cols
      rows = Missing_data_map[[set]]$Y_obs
      if(length(cols) == 0 || length(rows) == 0) next

      F[rows,] = sample_factors_scores_c(Eta_tilde[rows,cols],
                                         prior_mean[rows,],
                                         Lambda[cols,],
                                         resid_Eta_prec[cols],
                                         F_e_prec)
    }

  }))
  current_state = current_state[current_state_names]

  return(current_state)
}

