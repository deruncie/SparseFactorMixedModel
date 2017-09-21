sample_Lambda_B = function(BSFG_state,grainSize = 1,...) {
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
    n_coefs = b + k
    prior_mean = matrix(0,n_coefs,p)
    if(b > 0) {
      prior_prec = rbind(B_prec,t(Plam))
    } else{ # b == 0
      prior_prec = t(Plam)
    }
    for(set in seq_along(Missing_data_map)){
      cols = Missing_data_map[[set]]$Y_cols
      rows = Missing_data_map[[set]]$Y_obs
      if(length(cols) == 0 || length(rows) == 0) next

      if(is.null(cis_genotypes)){
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
      } else{
        result = sample_MME_fixedEffects_cis_c(Qt_list[[set]] %**% Eta[rows,cols],
                                               cbind(QtX_list[[set]], Qt_list[[set]] %**% F[rows,,drop=FALSE]),
                                               Qt_cis_genotypes_list[[set]],
                                               Sigma_Choleskys_list[[set]],
                                               resid_h2_index[cols],
                                               tot_Eta_prec[cols],
                                               prior_mean[,cols],
                                               prior_prec[,cols],
                                               cis_effects_index[cols],
                                               sum(n_cis_effects),
                                               grainSize
        )
        coefs = result[[1]]
        if(b > 0) {
          B[,cols] = coefs[1:b,,drop=FALSE]
        }
        Lambda[cols,] = t(coefs[b + 1:k,,drop=FALSE])
        cis_index = do.call(c,lapply(cols,function(x) if(n_cis_effects[x] > 0) cis_effects_index[x] + 1:n_cis_effects[x]-1))
        cis_effects[cis_index] = result[[2]][cis_index]
      }
    }
  }))
  current_state = current_state[current_state_names]

  return(current_state)
}

