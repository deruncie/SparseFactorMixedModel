sample_Lambda_B.general_BSFG = function(BSFG_state,grainSize = 1,...) {
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
      coefs = sample_MME_fixedEffects(Eta,Design,Sigma_Choleskys, resid_h2_index, tot_Eta_prec, prior_mean, prior_prec,grainSize)
      if(b > 0){
        B[] = coefs[1:b,,drop=FALSE]
      }
      Lambda[] = t(coefs[b + 1:k,,drop=FALSE])
    } else{
      for(j in 1:p){
        cis_X_j = cis_genotypes[[j]]
        Design_j = cbind(Design,cis_X_j)
        prior_mean_j = rbind(prior_mean[,j,drop=FALSE],0)
        prior_prec_j = rbind(prior_prec[,j,drop=FALSE],1e-10)
        coefs_j = sample_MME_fixedEffects(Eta[,j,drop=FALSE],Design_j,Sigma_Choleskys,  resid_h2_index[j], tot_Eta_prec[,j,drop=FALSE], prior_mean_j, prior_prec_j,grainSize)
        if(b > 0){
          B[,j] = coefs_j[1:b]
        }
        Lambda[j,] = coefs_j[b+1:k]
        cis_effects[,cis_effects_index[j]] = coefs_j[-c(1:(b+k))]
      }
    }
  }))
  current_state = current_state[current_state_names]

  return(current_state)
}

