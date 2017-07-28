sample_latent_traits.faster_BSFG = function(BSFG_state,grainSize,...) {
  data_matrices  = BSFG_state$data_matrices
  priors         = BSFG_state$priors
  run_parameters = BSFG_state$run_parameters
  run_variables  = BSFG_state$run_variables
  current_state  = BSFG_state$current_state
  invert_aI_bZKZ = BSFG_state$run_variables$invert_aI_bZKZ
  Ut = t(invert_aI_bZKZ$U)
  s = invert_aI_bZKZ$s

  current_state_names = names(current_state)
  a =  with(c(priors,run_parameters, run_variables,data_matrices,current_state), {
                sample_latent_traits_fast(X,X_F,Z,h2s_matrix,X_F_zero_variance,
                              Eta,Lambda,F,B,B_F,U_R,U_F,resid_h2,F_h2,tot_Eta_prec,tot_F_prec,B_F_prec,
                              tot_Eta_prec_shape,tot_Eta_prec_rate,tot_F_prec_shape,tot_F_prec_rate,
                              h2_priors_resids,h2_priors_factors,
                              Ut,s,invert_aZZt_Kinv)
  })
  current_state$tot_Eta_prec
  return(current_state)
}
