sample_latent_traits = function(BSFG_state,grainSize = 1,...) {
  data_matrices  = BSFG_state$data_matrices
  priors         = BSFG_state$priors
  run_parameters = BSFG_state$run_parameters
  run_variables  = BSFG_state$run_variables
  current_state  = BSFG_state$current_state

  current_state_names = names(current_state)
  current_state = with(c(priors,run_parameters, run_variables,data_matrices),within(current_state, {
    k = ncol(Lambda)

    # -----Sample tot_Eta_prec, resid_h2, U_R ---------------- #
    #conditioning on B, F, Lambda, resid_h2, tot_Eta_prec
    Eta_tilde = Eta - XB - F %*% t(Lambda)
    scores = rep(0,p)
    # QtEta_tilde = Qt_list[[1]] %**% Eta_tilde
    # scores = tot_prec_scores(QtEta_tilde,
    #                                 Sigma_Choleskys_list[[1]],
    #                                 resid_h2_index,
    #                                 grainSize
    #                                )
    # tot_Eta_prec[] = rgamma(1,shape = tot_Eta_prec_shape + (n*p)/2, rate = tot_Eta_prec_rate + 0.5*sum(scores))

    # sample columns of tot_Eta_prec, resid_h2, U_R in sets with the same patterns of missing data
    for(set in seq_along(Missing_data_map)){
      cols = Missing_data_map[[set]]$Y_cols
      rows = Missing_data_map[[set]]$Y_obs
      if(length(cols) == 0 || length(rows) == 0) next

      QtEta_tilde_set = Qt_list[[set]] %**% Eta_tilde[rows,cols,drop=FALSE]
      scores[cols] = tot_prec_scores(QtEta_tilde_set,
                                     Sigma_Choleskys_list[[set]],
                                     resid_h2_index[cols],
                                     grainSize
                                    )
      shape = tot_Eta_prec_shape + length(rows)/2
      if(b_QTL > 0) {
        shape = shape + b_QTL/2
        scores[cols] = scores[cols] + colSums(B_QTL[,cols,drop=FALSE]^2 * B_QTL_prec[,cols,drop=FALSE])
      }
      if(lambda_propto_Vp) {
        # include tot_Eta_prec in prior of Lambda
        scores[cols] = scores[cols] + rowSums(Lambda[cols,,drop=FALSE]^2*Lambda_prec[cols,,drop=FALSE]*tauh[rep(1,length(cols)),])
        shape = shape + k/2
      }
      if(cauchy_sigma_tot){
        # switch to C+(0,1) on sigma = 1/sqrt(prec)
        for(i in 1:100){
          if(!'tot_Eta_prec_xi' %in% ls()) {
            print('init')
            tot_Eta_prec_xi = rep(1,p)
          }
          tot_Eta_prec[cols] = rgamma(length(cols),shape = shape - tot_Eta_prec_shape + 1/2, rate = 1/tot_Eta_prec_xi[cols] + 0.5*scores[cols])
          tot_Eta_prec_xi[cols] = 1/rgamma(length(cols),shape = 1,rate = 1 + tot_Eta_prec[cols])
        }
      } else{
        tot_Eta_prec[cols] = rgamma(length(cols),shape = shape, rate = tot_Eta_prec_rate[cols] + 0.5*scores[cols])
      }

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

      U_R[,cols] = sample_MME_ZKZts_c(Eta_tilde[rows,cols,drop=FALSE],
                                    ZL[rows,,drop=FALSE],
                                    tot_Eta_prec[cols],
                                    randomEffect_C_Choleskys_list[[set]],
                                    resid_h2[,cols,drop=FALSE],
                                    resid_h2_index[cols],
                                    grainSize)
    }


    resid_Eta_prec = tot_Eta_prec / (1-colSums(resid_h2))

    # -----Sample Lambda and B_F ------------------ #
    # F, marginalizing over random effects (conditional on F_h2, tot_F_prec)
    # only use rows in Missing_data_map[[1]]. These are all the rows with non-missing data
    rows = Missing_data_map[[1]]$Y_obs
    F_tilde = F
    if(b_F > 0){
      if(b_QTL_F > 0){
        F_tilde = F_tilde - QTL_factors_Z %**% (QTL_factors_X %*% B_QTL_F)
      }
      prior_mean = matrix(0,b_F,k)
      prior_prec = B_F_prec * tot_F_prec[rep(1,b_F),,drop=FALSE]  # prior for B_F includes tot_F_prec
      B_F[] = sample_MME_fixedEffects_c(Qt_list[[1]] %**% F_tilde[rows,,drop=FALSE],
                                        Qt1_XF,
                                        Sigma_Choleskys_list[[1]],
                                        F_h2_index,
                                        tot_F_prec,
                                        prior_mean,
                                        prior_prec,
                                        grainSize)
      XFBF = X_F %**% B_F
      F_tilde = F - XFBF
    }
    if(b_QTL_F > 0){
      prior_mean = matrix(0,b_QTL_F,k)
      prior_prec = B_QTL_F_prec# * tot_F_prec[rep(1,b_QTL_F),,drop=FALSE]  # prior for B_QTL_F includes tot_F_prec
      B_QTL_F[] = sample_MME_fixedEffects_hierarchical_c(Qt_list[[1]] %**% F_tilde[rows,,drop=FALSE],
                                                       Qt1_QTL_Factors_Z,  # only includes rows of Qt1
                                                       QTL_factors_X,
                                                       Sigma_Choleskys_list[[1]],
                                                       F_h2_index,
                                                       tot_F_prec,
                                                       prior_mean,
                                                       prior_prec,
                                                       grainSize)

      XFBF = X_F %**% B_F + QTL_factors_Z %**% (QTL_factors_X %*% B_QTL_F)
      F_tilde = F - XFBF
    }

    # -----Sample tot_F_prec, F_h2, U_F ---------------- #
    #conditioning on B, F, Lambda, F_h2, tot_F_prec

    Qt_F_tilde = Qt_list[[1]] %**% F_tilde[rows,,drop=FALSE]

    scores = tot_prec_scores(Qt_F_tilde,
                             Sigma_Choleskys_list[[1]],
                             F_h2_index,
                             grainSize)
    shape = tot_F_prec_shape + n/2
    rate = tot_F_prec_rate + scores/2
    if(b_F > 0) {
      shape = shape + sum(!X_F_zero_variance)/2
      rate = rate + colSums((B_F^2*B_F_prec)[!X_F_zero_variance,,drop=FALSE])/2   # add this if tot_F_prec part of the prior for B_F
    }
    # if(b_QTL_F > 0) {
    #   shape = shape + b_QTL_F/2
    #   rate = rate + colSums(B_QTL_F^2 * B_QTL_F_prec)/2
    # }
    tot_F_prec[] = rgamma(k,shape = shape, rate = rate)
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

    U_F[] = sample_MME_ZKZts_c(F_tilde[rows,,drop=FALSE], ZL[rows,,drop=FALSE], tot_F_prec, randomEffect_C_Choleskys_list[[1]], F_h2, F_h2_index,grainSize)

    # -----Sample F----------------------- #
    #conditioning on B, U_F,U_R,Lambda, F_h2
    Eta_tilde = Eta - XB - ZL %**% U_R
    F_e_prec = tot_F_prec / (1-colSums(F_h2))
    prior_mean = ZL %**% U_F
    if(b_F + b_QTL_F > 0) {
      prior_mean = prior_mean + XFBF
    }

    for(set in seq_along(Missing_row_data_map)){
      cols = Missing_row_data_map[[set]]$Y_cols
      rows = Missing_row_data_map[[set]]$Y_obs
      if(length(cols) > 0) {
        F[rows,] = sample_factors_scores_c(Eta_tilde[rows,cols,drop=FALSE],
                                           prior_mean[rows,,drop=FALSE],
                                           Lambda[cols,,drop=FALSE],
                                           resid_Eta_prec[cols],
                                           F_e_prec)
      } else{
        F[rows,] = prior_mean[rows,,drop=FALSE] + sweep(rstdnorm_mat(length(rows),k),2,sqrt(F_e_prec[1,]),'/')
      }
    }

  }))
  # current_state = current_state[current_state_names]

  return(current_state)
}

