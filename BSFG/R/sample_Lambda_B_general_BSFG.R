# Eta = cbind(X, F) * rbind(t(B),t(Lambda) + Z*U_F + U_E
# Qt*Eta = cbind(Qt*X,Qt*F) * rbind(t(B),t(Lambda) + Qt*Z*U_F + Qt*U_E
# note: different Qt for each column


sample_Lambda_B.general_BSFG = function(BSFG_state,grainSize = 1,...) {
  data_matrices  = BSFG_state$data_matrices
  priors         = BSFG_state$priors
  run_parameters = BSFG_state$run_parameters
  run_variables  = BSFG_state$run_variables
  current_state  = BSFG_state$current_state
  invert_aI_bZKZ = BSFG_state$run_variables$invert_aI_bZKZ
  Qt = t(invert_aI_bZKZ$Q)
  s = invert_aI_bZKZ$s

  current_state_names = names(current_state)
  current_state = with(c(priors,run_parameters, run_variables,data_matrices),within(current_state, {
    k = ncol(Lambda)

    # -----Sample Lambda and B ------------------ #
    #conditioning on W, F, marginalizing over random effects (conditional on resid_h2)
    Design = cBind(X,F)
    n_coefs = b + k
    prior_mean = matrix(0,n_coefs,p)
    if(b > 0) {
      prior_prec = rbind(B_prec,t(Plam))
    } else{ # b == 0
      prior_prec = t(Plam)
    }
    if(is.null(cis_genotypes)){
      for(set in seq_along(Missing_data_map)){
        cols = Missing_data_map[[set]]$Y_cols
        rows = Missing_data_map[[set]]$Y_obs
        coefs = sample_MME_fixedEffects_c(Qt_list[[set]] %**% Eta[rows,cols],
                                          cbind(QtX_list[[set]], Qt_list[[set]] %**% F[rows,,drop=FALSE]),
                                          Sigma_Choleskys_list[[set]],
                                          resid_h2_index[cols],
                                          tot_Eta_prec[cols],
                                          prior_mean[,cols],
                                          prior_prec[,cols],
                                          grainSize
                                          )
        if(b > 0) {
          B[,cols] = coefs[1:b,,drop=FALSE]
        }
        Lambda[cols,] = t(coefs[b + 1:k,,drop=FALSE])
      }
      # resid_prec = uncorrelated_prec_mat(resid_h2,tot_Eta_prec,s)
      # coefs = sample_coefs_parallel_sparse_c_Eigen(Qt %**% Eta,Qt %**% Design,resid_prec, prior_mean,prior_prec,grainSize)
      # if(b > 0){
      #   B[] = coefs[1:b,,drop=FALSE]
      # }
      # Lambda[] = t(coefs[b + 1:k,,drop=FALSE])
    } else{
      stop('cis effects not implemented yet')
      # result = sample_MME_fixedEffects_cis(Eta,Design,cis_genotypes,cis_effects_index,Sigma_Choleskys, resid_h2_index, tot_Eta_prec, prior_mean, prior_prec,grainSize)
      # if(b > 0){
      #   B[] = result[[1]][1:b,,drop=FALSE]
      # }
      # Lambda[] = t(result[[1]][b+1:k,,drop=FALSE])
      # cis_effects[] = result[[2]]
    }
  }))
  current_state = current_state[current_state_names]

  return(current_state)
}

