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
    Design = cBind(X,F)
    rows = b + k
    prior_mean = matrix(0,rows,p)
    if(b > 0) {
      prior_prec = rbind(B_prec,t(Plam))
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
      result = sample_MME_fixedEffects_cis(Eta,Design,cis_genotypes,cis_effects_index,Sigma_Choleskys, resid_h2_index, tot_Eta_prec, prior_mean, prior_prec,grainSize)
      if(b > 0){
        B[] = result[[1]][1:b,,drop=FALSE]
      }
      Lambda[] = t(result[[1]][b+1:k,,drop=FALSE])
      cis_effects[] = result[[2]]
    }
  }))
  current_state = current_state[current_state_names]

  return(current_state)
}

