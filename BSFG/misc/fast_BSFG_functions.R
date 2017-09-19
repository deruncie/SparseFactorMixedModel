# sample_h2s_discrete_MH_fast = function(UtEta,h2_divisions
#   Y,tot_Eta_prec, Sigma_Choleskys,discrete_priors,h2s_matrix,h2_index,step_size,grainSize){
#   # while this works, it is much slower than doing the full scan over all traits, at least for multiple traits
#   # testing with p=100, solving the whole set takes ~4-5x solving just 1. And this method requires doing each trait separately
#   # both methods are similarly easy to multiplex, so no advantage there either.
#   n = nrow(Y)
#   p = ncol(Y)
#   discrete_bins = length(discrete_priors)
#
#
#   # sample runif(p,0,1) before because parallel RNGs aren't consecutive.
#   r_draws = runif(p)
#
#   h2s_index = sample_h2s_discrete_MH_c(Y,tot_Eta_prec,discrete_priors,h2_index,h2s_matrix,Sigma_Choleskys,r_draws,step_size,grainSize)+1
#
#   return(h2s_index)
# }

sample_h2s_discrete_fast = function(UtEta, tot_Eta_prec, discrete_priors,s,grainSize) {
  n = nrow(UtEta)
  p = ncol(UtEta)

  log_ps = log_p_h2s_fast(UtEta, tot_Eta_prec, discrete_priors,s,grainSize)
  rs = runif(p)
  h2s_index = sample_h2s(log_ps,rs,grainSize)
  return(h2s_index)
}
