sample_B2_prec_reg_horseshoe = function(BSFG_state,...) {
  # sampling as described in Supplemental methods, except we multiply columns of Prec_lambda by delta
  # phi2 = \lambda^2 in methods
  # the delta sequence controls tau_k. We have tau_1~C+(0,tau_0), and tau_k = tau_1*prod_{l=1}^K(delta^{-1}_l)
  # delta_l controls the decrease in odds of inclusion of each element for the lth factor relative to the (l-1)th
  # Note:  Piironen and Vehtari do not include sigma^2 in prior, so I have removed
  priors         = BSFG_state$priors
  run_variables  = BSFG_state$run_variables
  run_parameters = BSFG_state$run_parameters
  current_state  = BSFG_state$current_state

  current_state = with(c(priors,run_variables,run_parameters),
                       with(B2_prior,{

                         tau_0 = prop_0/(1-prop_0) * 1/sqrt(n)

                         c2_shape = c2$nu-1
                         c2_rate = c2$V*c2$nu

                         within(current_state,{

                           # initialize variables if needed
                           if(!'B2_tau2' %in% names(current_state)){
                             if(verbose) print('initializing B_prec regularized horseshoe')
                             B2_tau2 = matrix(1,1,1)
                             B2_xi = matrix(1,1,1)
                             B2_R_phi2 = matrix(1,b2_R,p)
                             B2_R_nu = matrix(1,b2_R,p)
                             B2_F_phi2 = matrix(1,b2_F,K)
                             B2_F_nu = matrix(1,b2_F,K)
                             B2_c2 = matrix(1,1,1)
                           }

                           B2_R_2 = B2_R^2
                           B2_R_2_std = sweep(B2_R_2,2,tot_Eta_prec[1,],'*')
                           B2_F_2 = B2_F^2
                           B2_F_2_std = sweep(B2_F_2,2,tot_F_prec[1,],'*')

                           B2_R_nu[] = matrix(1/rgamma(b2_R*p,shape = 1, rate = 1 + 1/B2_R_phi2), nr = b2_R, nc = p)
                           B2_R_phi2[] = matrix(1/rgamma(b2_R*p,shape = 1, rate = 1/B2_R_nu + B2_R_2_std / (2*B2_tau2[1])),nr=b2_R,nc = p)
                           B2_F_nu[] = matrix(1/rgamma(b2_F*K,shape = 1, rate = 1 + 1/B2_F_phi2), nr = b2_F, nc = K)
                           B2_F_phi2[] = matrix(1/rgamma(b2_F*K,shape = 1, rate = 1/B2_F_nu + B2_F_2_std / (2*B2_tau2[1])),nr=b2_F,nc = K)

                           B2_xi[] = 1/rgamma(1,shape=1,rate=1/tau_0^2 + 1/B2_tau2[1])
                           B2_tau2[] = 1/rgamma(1,shape = (b2_R*p + b2_F*K + 1)/2, rate = 1/B2_xi[1] + (sum(B2_R_2_std/B2_R_phi2) + sum(B2_F_2_std/B2_F_phi2))/2)

                           B2_c2[] = 1/rgamma(1,shape = c2_shape + (b2_R*p + b2_F*K)/2,rate = c2_rate + (sum(B2_R_2) + sum(B2_F_2))/2)

                           # -----Update Plam-------------------- #
                           B2_R_prec = (B2_c2[1] + B2_tau2[1]*B2_R_phi2) / (B2_c2[1] * B2_tau2[1] * B2_R_phi2)
                           B2_F_prec = (B2_c2[1] + B2_tau2[1]*B2_F_phi2) / (B2_c2[1] * B2_tau2[1] * B2_F_phi2)

                           rm(list=c('B2_R_2','B2_R_2_std','B2_F_2','B2_F_2_std'))
                         })
                       }))
  return(current_state)
}

sample_B_prec_RE = function(BSFG_state,...){
  # treats B as a random effect - no parameter-specific shrinkage
  # error if ncol(B_F)>0
  priors         = BSFG_state$priors
  run_variables  = BSFG_state$run_variables
  current_state  = BSFG_state$current_state

  current_state = with(c(priors,run_variables),
                       with(B_prior,{
                         #   list(
                         #   # load priors
                         #   B_df   = B_prior$B_df,
                         #   B_F_df = B_prior$B_F_df,
                         #   B_QTL_df   = B_prior$B_QTL_df,
                         #   B_QTL_F_df = B_prior$B_QTL_F_df,
                         #   separate_QTL_shrinkage = B_prior$separate_QTL_shrinkage
                         # ),
                         if(b > 0) {
                           prec_shape = with(global, nu - 1)
                           prec_rate = with(global, V * nu)
                         }
                         if(b_F > 0) {
                           stop("sample_B_prec_RE doesn't work with B_F")
                         }
                         within(current_state,{

                           # initialize variables if needed
                           if(!exists('B_prec')){
                             if(b > 0) {
                               B_prec = matrix(rgamma(b,shape = prec_shape,rate=prec_rate),nrow = b, ncol = p)
                             } else{
                               B_prec = matrix(0,nrow=0,ncol=p)
                             }
                           }
                           B2 = B^2

                           B_prec[] = matrix(rgamma(b,shape = prec_shape + p/2,rate = prec_rate + rowSums(B2)/2),nrow = b, ncol = p)

                           if(length(resid_intercept) > 0){
                             B_prec[1,resid_intercept] = 1e-10
                           }
                         })
                       }))
  return(current_state)
}



sample_B_prec_ARD = function(BSFG_state,...){
  priors         = BSFG_state$priors
  run_variables  = BSFG_state$run_variables
  current_state  = BSFG_state$current_state

  current_state = with(c(priors,run_variables),
                       with(B_prior,{
                       #   list(
                       #   # load priors
                       #   B_df   = B_prior$B_df,
                       #   B_F_df = B_prior$B_F_df,
                       #   B_QTL_df   = B_prior$B_QTL_df,
                       #   B_QTL_F_df = B_prior$B_QTL_F_df,
                       #   separate_QTL_shrinkage = B_prior$separate_QTL_shrinkage
                       # ),
                              if(b > 0) {
                                tau_shape = with(global, nu - 1)
                                tau_rate = with(global, V * nu)
                              }
                              if(b_F > 0) {
                                tau_F_shape = with(global_F, nu - 1)
                                tau_F_rate = with(global_F, V * nu)
                              }
                              within(current_state,{

                                 # initialize variables if needed
                                 if(!exists('B_tau')){
                                   if(b > 0) {
                                     B_tau = matrix(c(1e-10,rgamma(b-1,shape = tau_shape, rate = tau_rate)),nrow=1)
                                   } else{
                                     B_tau = matrix(0,nrow=1,ncol=0)
                                   }
                                   if(b_F > 0) {
                                     B_F_tau = matrix(rgamma(b_F,shape = tau_F_shape, rate = tau_F_rate),nrow=1)
                                   } else{
                                     B_F_tau = matrix(0,nrow=1,ncol=0)
                                   }

                                   B_prec = matrix(B_tau,nrow = b, ncol = p)
                                   B_F_prec = matrix(B_F_tau,nrow = b_F, ncol = k)
                                 }
                                 B2 = B^2
                                 B_F2 = B_F^2 * tot_F_prec[rep(1,b_F),]  # need to account for tot_F_prec

                                 if(ncol(fixed_effects_common)>0){
                                   ib = fixed_effects_common[1,]
                                   ib_F = fixed_effects_common[2,]
                                   tau = rgamma(length(ib),
                                                shape = tau_shape + ncol(B2)/2 + ncol(B_F2)/2,
                                                rate = tau_rate +
                                                  rowSums((B2[ib,,drop=FALSE] * B_prec[ib,,drop=FALSE]/c(B_tau)[ib]))/2 +
                                                  rowSums(B_F2[ib_F,,drop=FALSE] * B_F_prec[ib_F,,drop=FALSE]/c(B_F_tau)[ib_F])/2)
                                   B_tau[1,ib] = tau
                                   B_F_tau[1,ib_F] = tau
                                 }
                                 if(length(fixed_effects_only_resid) > 0){
                                   ib = fixed_effects_only_resid
                                   tau = rgamma(length(ib),
                                                shape = tau_shape + ncol(B2)/2,
                                                rate = tau_rate +
                                                  rowSums((B2[ib,,drop=FALSE] * B_prec[ib,,drop=FALSE]/c(B_tau)[ib]))/2)
                                   B_tau[1,ib] = tau

                                 }
                                 if(length(fixed_effects_only_factors) > 0){
                                   ib = fixed_effects_only_factors
                                   tau = rgamma(length(ib),
                                                shape = tau_shape + ncol(B_F2)/2,
                                                rate = tau_rate +
                                                  rowSums(B_F2[ib,,drop=FALSE] * B_F_prec[ib,,drop=FALSE]/c(B_F_tau)[ib])/2)
                                   B_F_tau[1,ib] = tau

                                 }
                                 B_prec[] = matrix(rgamma(b*p,shape = (B_df + 1)/2,rate = (B_df + B2*c(B_tau))/2),nr = b,nc = p)
                                 B_F_prec[] = matrix(rgamma(b_F*k,shape = (B_F_df + 1)/2,rate = (B_F_df + B_F2*t(B_F_tau)[,rep(1,k)])/2),nr = b_F,nc = k)
                                 B_prec[] = B_prec*c(B_tau)
                                 if(resid_intercept){
                                   B_tau[1,1] = 1e-10
                                   B_prec[1,] = 1e-10
                                 }
                                 B_F_prec[] = B_F_prec*t(B_F_tau)[,rep(1,k)]
                                 B_F_tau[1,X_F_zero_variance] = 1e10
                                 B_F_prec[X_F_zero_variance,] = 1e10
                               })
                         }))
  return(current_state)
}

#
#
# sample_B_prec = function(BSFG_state,...){
#   priors         = BSFG_state$priors
#   run_variables  = BSFG_state$run_variables
#   current_state  = BSFG_state$current_state
#
#   current_state = with(c(priors,run_variables),within(current_state,{
#     if(b > 0) {
#       if(same_fixed_model) {  # want same tau for both
#         B2 = cbind(B,B_F)^2  # tauh
#       } else{
#         B2 = B^2
#       }
#       B_tau[1,] = rgamma(b, shape = fixed_resid_prec_shape + ncol(B2)/2, rate = fixed_resid_prec_rate + rowSums(B2)/2)
#       if(resid_intercept){
#         B_tau[1,1] = 1e-10
#       }
#       B_prec = matrix(B_tau,nrow = b, ncol = p)
#     }
#     if(b_F > 0){
#       if(same_fixed_model) {
#         B_F_tau[1,] = B_tau[1,]
#       } else{
#         B_F2 = B_F^2
#         B_F_tau[1,] = rgamma(b_F, shape = fixed_factors_prec_shape + ncol(B_F2)/2, rate = fixed_factors_prec_rate + rowSums(B_F2)/2)
#       }
#       B_F_tau[1,X_F_zero_variance] = 1e10
#       B_F_prec = matrix(B_F_tau,nrow = b_F, ncol = k)
#     }
#   }))
#   return(current_state)
# }
#
# sample_B_prec_ARD = function(BSFG_state,...){
#   priors         = BSFG_state$priors
#   run_variables  = BSFG_state$run_variables
#   current_state  = BSFG_state$current_state
#
#   current_state = with(c(priors,run_variables),
#                        with(list(
#                          # load priors
#                          B_df   = B_prior$B_df,
#                          B_F_df = B_prior$B_F_df
#                        ),within(current_state,{
#
#                          # initialize variables if needed
#                          if(!exists('B_tau')){
#                            if(b > 0) {
#                              B_tau = matrix(c(1e-10,rgamma(b-1,shape = fixed_resid_prec_shape, rate = fixed_resid_prec_rate)),nrow=1)
#                            } else{
#                              B_tau = matrix(0,nrow=1,ncol=0)
#                            }
#                            if(b_F > 0) {
#                              B_F_tau = matrix(rgamma(b_F,shape = fixed_factors_prec_shape, rate = fixed_factors_prec_rate),nrow=1)
#                            } else{
#                              B_F_tau = matrix(0,nrow=1,ncol=0)
#                            }
#
#                            B_prec = matrix(B_tau,nrow = b, ncol = p)
#                            B_F_prec = matrix(B_F_tau,nrow = b_F, ncol = k)
#                          }
#                          if(b > 0) {
#                            B2 = B^2
#                            B_tau[1,] = rgamma(b, shape = fixed_resid_prec_shape + ncol(B2)/2, rate = fixed_resid_prec_rate + rowSums((B2 * B_prec/c(B_tau)))/2)
#                            B_prec[] = matrix(rgamma(b*p,shape = (B_df + 1)/2,rate = (B_df + B2*c(B_tau))/2),nr = b,nc = p)
#                            B_prec[] = B_prec*c(B_tau)
#                            if(resid_intercept){
#                              B_tau[1,1] = 1e-10
#                              B_prec[1,] = 1e-10
#                            }
#                          }
#                          if(b_F > 0) {
#                            B_F2 = B_F^2 * tot_F_prec[rep(1,b_F),]  # need to account for tot_F_prec
#
#                            # B_F_tau[1,] = rgamma(b_F, shape = fixed_factors_prec_shape + ncol(B_F2)/2, rate = fixed_factors_prec_rate + rowSums(B_F2 * B_F_prec/c(B_F_tau))/2)
#                            # B_F_prec[] = matrix(rgamma(b_F*k,shape = (B_F_df + 1)/2,rate = (B_F_df + B_F2*c(B_F_tau))/2),nr = b_F,nc = k)
#
#                            B_F_tau[1,] = rgamma(b_F, shape = fixed_factors_prec_shape + ncol(B_F2)/2, rate = fixed_factors_prec_rate + rowSums(B_F2 * B_F_prec/t(B_F_tau)[,rep(1,k)])/2)
#                            B_F_prec[] = matrix(rgamma(b_F*k,shape = (B_F_df + 1)/2,rate = (B_F_df + B_F2*t(B_F_tau)[,rep(1,k)])/2),nr = b_F,nc = k)
#
#                            B_F_prec[] = B_F_prec*t(B_F_tau)[,rep(1,k)]
#                            B_F_tau[1,X_F_zero_variance] = 1e10
#                            B_F_prec[X_F_zero_variance,] = 1e10
#                          }
#                        })))
#   return(current_state)
# }
#
# sample_B_prec_combined_ARD = function(BSFG_state,...){
#   priors         = BSFG_state$priors
#   run_variables  = BSFG_state$run_variables
#   current_state  = BSFG_state$current_state
#
#   current_state = with(c(priors,run_variables),
#                        with(list(
#                          # load priors
#                          B_df   = B_prior$B_df,
#                          B_F_df = B_prior$B_F_df,
#                          B_QTL_df   = B_prior$B_QTL_df,
#                          B_QTL_F_df = B_prior$B_QTL_F_df,
#                          separate_QTL_shrinkage = B_prior$separate_QTL_shrinkage
#                        ),within(current_state,{
#
#                          # initialize variables if needed
#                          if(!exists('B_tau')){
#                            if(b > 0) {
#                              B_tau = matrix(c(1e-10,rgamma(b-1,shape = fixed_resid_prec_shape, rate = fixed_resid_prec_rate)),nrow=1)
#                            } else{
#                              B_tau = matrix(0,nrow=1,ncol=0)
#                            }
#                            if(b_F > 0) {
#                              B_F_tau = matrix(rgamma(b_F,shape = fixed_factors_prec_shape, rate = fixed_factors_prec_rate),nrow=1)
#                            } else{
#                              B_F_tau = matrix(0,nrow=1,ncol=0)
#                            }
#
#                            B_prec = matrix(B_tau,nrow = b, ncol = p)
#                            B_F_prec = matrix(B_F_tau,nrow = b_F, ncol = k)
#                          }
#                          B2 = B^2
#                          B_F2 = B_F^2 * tot_F_prec[rep(1,b_F),]  # need to account for tot_F_prec
#
#                          if(ncol(fixed_effects_common)>0){
#                            ib = fixed_effects_common[1,]
#                            ib_F = fixed_effects_common[2,]
#                            tau = rgamma(length(ib),
#                                        shape = fixed_resid_prec_shape[ib] + ncol(B2)/2 + ncol(B_F2)/2,
#                                        rate = fixed_resid_prec_rate[ib] +
#                                               rowSums((B2[ib,,drop=FALSE] * B_prec[ib,,drop=FALSE]/c(B_tau)[ib]))/2 +
#                                               rowSums(B_F2[ib_F,,drop=FALSE] * B_F_prec[ib_F,,drop=FALSE]/c(B_F_tau)[ib_F])/2)
#                            B_tau[1,ib] = tau
#                            B_F_tau[1,ib_F] = tau
#                          }
#                          if(length(fixed_effects_only_resid) > 0){
#                            ib = fixed_effects_only_resid
#                            tau = rgamma(length(ib),
#                                         shape = fixed_resid_prec_shape[ib] + ncol(B2)/2,
#                                         rate = fixed_resid_prec_rate[ib] +
#                                           rowSums((B2[ib,,drop=FALSE] * B_prec[ib,,drop=FALSE]/c(B_tau)[ib]))/2)
#                            B_tau[1,ib] = tau
#
#                          }
#                          if(length(fixed_effects_only_factors) > 0){
#                            ib = fixed_effects_only_factors
#                            tau = rgamma(length(ib),
#                                         shape = fixed_resid_prec_shape[ib] + ncol(B_F2)/2,
#                                         rate = fixed_resid_prec_rate[ib] +
#                                           rowSums(B_F2[ib,,drop=FALSE] * B_F_prec[ib,,drop=FALSE]/c(B_F_tau)[ib])/2)
#                            B_F_tau[1,ib] = tau
#
#                          }
#                          B_prec[] = matrix(rgamma(b*p,shape = (B_df + 1)/2,rate = (B_df + B2*c(B_tau))/2),nr = b,nc = p)
#                          B_F_prec[] = matrix(rgamma(b_F*k,shape = (B_F_df + 1)/2,rate = (B_F_df + B_F2*t(B_F_tau)[,rep(1,k)])/2),nr = b_F,nc = k)
#                          B_prec[] = B_prec*c(B_tau)
#                          if(resid_intercept){
#                            B_tau[1,1] = 1e-10
#                            B_prec[1,] = 1e-10
#                          }
#                          B_F_prec[] = B_F_prec*t(B_F_tau)[,rep(1,k)]
#                          B_F_tau[1,X_F_zero_variance] = 1e10
#                          B_F_prec[X_F_zero_variance,] = 1e10
#                        })))
#   return(current_state)
# }
#
# sample_B_prec_ARD_v2 = function(BSFG_state,...){
#   priors         = BSFG_state$priors
#   run_variables  = BSFG_state$run_variables
#   current_state  = BSFG_state$current_state
#
#   current_state = with(c(priors,run_variables),
#                        with(list(
#                          # load priors
#                          B_df   = B_prior$B_df,
#                          B_F_df = B_prior$B_F_df
#                        ),within(current_state,{
#
#                          # initialize variables if needed
#                          if(!exists('B_tau')){
#                            if(b > 0) {
#                              B_tau = matrix(c(1e-10,rgamma(b-1,shape = fixed_resid_prec_shape, rate = fixed_resid_prec_rate)),nrow=1)
#                            } else{
#                              B_tau = matrix(0,nrow=1,ncol=0)
#                            }
#                            if(b_F > 0) {
#                              B_F_beta = matrix(rgamma(b_F,shape = fixed_factors_prec_shape, rate = 1/fixed_factors_prec_rate),nrow=1)
#                            } else{
#                              B_F_beta = matrix(0,nrow=1,ncol=0)
#                            }
#
#                            B_prec = matrix(B_tau,nrow = b, ncol = p)
#                            B_F_prec = matrix(B_F_beta,nrow = b_F, ncol = k)
#                          }
#                          if(b > 0) {
#                            B2 = B^2
#                            B_tau[1,] = rgamma(b, shape = fixed_resid_prec_shape + ncol(B2)/2, rate = fixed_resid_prec_rate + rowSums((B2 * B_prec/c(B_tau)))/2)
#                            B_prec[] = matrix(rgamma(b*p,shape = (B_df + 1)/2,rate = (B_df + B2*c(B_tau))/2),nr = b,nc = p)
#                            B_prec[] = B_prec*c(B_tau)
#                            if(resid_intercept){
#                              B_tau[1,1] = 1e-10
#                              B_prec[1,] = 1e-10
#                            }
#                          }
#                          if(b_F > 0) {
#                            B_F2 = B_F^2 * tot_F_prec[rep(1,b_F),]  # need to account for tot_F_prec
#
#                            B_F_beta[1,] = rgamma(b_F,shape = fixed_factors_prec_shape + k*B_F_df, rate = 1/fixed_factors_prec_rate + rowSums(B_F_prec))
#                            B_F_prec[] = rgamma(b_F*k,shape = (B_F_df + 1)/2,rate = (t(B_F_beta)[,rep(1,k)] + B_F2)/2)
#
#                            B_F_prec[X_F_zero_variance,] = 1e10
#                          }
#                        })))
#   return(current_state)
# }
#
#
# sample_B_prec_TPB = function(BSFG_state,ncores = detectCores(),cluster=NULL,...){
#   # following https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4214276/
#   # uses GIGrvg for GIG
#   priors         = BSFG_state$priors
#   run_variables  = BSFG_state$run_variables
#   current_state  = BSFG_state$current_state
#
#   current_state = with(c(priors,run_variables),
#                        with(list(
#                          # load priors
#                          B_A = B_prior$B_A,
#                          B_B = B_prior$B_B,
#                          B_omega = B_prior$B_omega,
#
#                          B_F_A = B_prior$B_F_A,
#                          B_F_B = B_prior$B_F_B,
#                          B_F_omega = B_prior$B_F_omega,
#                          B_F_ncores = B_prior$B_F_ncores
#                        ), within(current_state,{
#                          # initialize variables if needed
#                          if(!exists('B_psi')){
#                            if(b > 0) {
#                              B_psi = matrix(rgamma(p,shape = 1/2,rate = B_omega),nrow=1)
#                              B_lambda = matrix(rgamma(b*p,shape = B_B,rate=B_psi[rep(1,b),]),b,p)
#                              B_prec = 1/matrix(rgamma(b*p,shape = B_A,rate = B_lambda),b,p)
#                            }
#                          }
#                          if(!exists('B_F_psi')){
#                            if(b_F > 0) {
#                              B_F_psi = matrix(rgamma(k,shape = 1/2,rate = B_F_omega),nrow=1)
#                              B_F_lambda = matrix(rgamma(b_F*k,shape = B_F_B,rate=B_F_psi[rep(1,b_F),]),b_F,k)
#                              B_F_prec = 1/matrix(rgamma(b_F*k,shape = B_F_A,rate = B_F_lambda),b_F,k)
#                            }
#                          }
#                          if(b > 0) { # B[i,j] ~ N(0,B_prec[i,j]^{-1})
#                            # recover()
#                            B_psi[1,] = rgamma(p,shape=b*B_B + 1/2, rate = colSums(B_lambda) + B_omega)
#                            B_lambda[] = rgamma(b*p,shape = B_A + B_B, rate = 1/B_prec + B_psi[rep(1,b),])
#                            B_prec[] = 1/rgig_multiple(b*p,lambda = rep(B_A-1/2,b*p),chi = c(B^2),psi = 2*c(B_lambda))
#                            B_prec[B_prec==0] = 1e-10
#                            if(resid_intercept){
#                              B_prec[1,] = 1e-10
#                            }
#                          }
#                          if(b_F > 0) { # B_F[i,j] ~ N(0,B_F_prec[i,j]^{-1} * tot_F_prec[j]^{-1})
#                            B_F_psi[1,] = 1e-7#rgamma(k,shape=sum(!X_F_zero_variance)*B_F_B + 1/2, rate = colSums(B_F_lambda[!X_F_zero_variance,]) + B_F_omega)
#                            B_F_lambda[] = rgamma(b_F*k,shape = B_F_A + B_F_B, rate = 1/B_F_prec + B_F_psi[rep(1,b_F),])
#                            B_F_prec[] = 1/rgig_multiple(b_F*k,lambda = rep(B_F_A-1/2,b_F*k),chi = c(B_F^2*tot_F_prec[rep(1,b_F),]),psi = 2*c(B_F_lambda))
#                            B_F_prec[B_F_prec < 1e-10] = 1e-10
#                            B_F_prec[B_F_prec > 1e10] = 1e10
#                            B_F_prec[X_F_zero_variance,] = 1e10
#                          }
#                        })))
#   return(current_state)
# }
#
#
# sample_B_prec_QTLBEN = function(BSFG_state,...){
#   priors         = BSFG_state$priors
#   run_variables  = BSFG_state$run_variables
#   current_state  = BSFG_state$current_state
#
#   QTL_columns_resid = BSFG_state$data_matrices$QTL_columns_resid
#   QTL_columns_factors = BSFG_state$data_matrices$QTL_columns_factors
#
#   current_state = with(c(priors,run_variables),within(current_state,{
#     # load priors
#     B_df   = B_prior$B_df
#     B_F_df = B_prior$B_F_df
#
#     # initialize variables if needed
#     if(!exists('B_tau')){
#       if(b > 0) {
#         B_tau = matrix(c(1e-10,rgamma(b-1,shape = priors$fixed_resid_prec_shape, rate = priors$fixed_resid_prec_rate)),nrow=1)
#       } else{
#         B_tau = matrix(0,ncol=0,nrow=1)
#       }
#       if(b_F > 0) {
#         B_F_tau = matrix(rgamma(b_F,shape = priors$fixed_factors_prec_shape, rate = priors$fixed_factors_prec_rate),nrow=1)
#       } else{
#         B_F_tau = matrix(0,ncol=0,nrow=1)
#       }
#
#       B_prec = matrix(B_tau,nrow = b, ncol = p)
#       B_F_prec = matrix(B_F_tau,nrow = b_F, ncol = k)
#     }
#
#     if(b > 0) {
#       non_QTL_resid = 1:nrow(B)
#       if(!is.null(QTL_columns_resid)){
#         non_QTL_resid = non_QTL_resid[-QTL_columns_resid]
#       }
#       B2 = B^2
#       B_tau[1,] = rgamma(b, shape = fixed_resid_prec_shape + ncol(B2)/2, rate = fixed_resid_prec_rate + rowSums((B2 * B_prec/c(B_tau)))/2)
#       # B_tau[1,non_QTL_resid] = rgamma(length(non_QTL_resid), shape = fixed_resid_prec_shape + ncol(B2)/2,
#       # rate = fixed_resid_prec_rate + rowSums((B2[non_QTL_resid,,drop=FALSE] * B_prec[non_QTL_resid,,drop=FALSE]/c(B_tau[non_QTL_resid])))/2)
#       B_prec[non_QTL_resid,] = matrix(rgamma(length(non_QTL_resid)*p,shape = (B_df + 1)/2,rate = (B_df + B2[non_QTL_resid,,drop=FALSE]*c(B_tau[non_QTL_resid]))/2),nr = length(non_QTL_resid),nc = p)
#       B_prec[non_QTL_resid,] = B_prec[non_QTL_resid,,drop=FALSE]*c(B_tau[non_QTL_resid])
#
#       # inverse Gaussian Bayesian Elastic Net from Li and Lin 2010
#       B_F_std = B_F[QTL_columns_resid,,drop=FALSE]
#       QTL_lambda1 = B_df
#       QTL_lambda2 = B_tau[QTL_columns_resid]
#       tau = matrix(1/rinvgauss(n=length(QTL_columns_resid)*p,
#                                mean = sqrt(QTL_lambda1)/(2*QTL_lambda2*abs(B_F_std)),
#                                shape = QTL_lambda1/(4*QTL_lambda2))+1,
#                    nr = length(QTL_columns_resid),nc = p)
#       B_prec[QTL_columns_resid,] = tau/(tau-1) * QTL_lambda2
#
#       if(resid_intercept){
#         B_tau[1,1] = 1e-10
#         B_prec[1,] = 1e-10
#       }
#     }
#     if(b_F > 0) {
#       non_QTL_factor = 1:nrow(B_F)
#       if(!is.null(QTL_columns_factors)){
#         non_QTL_factor = non_QTL_factor[-QTL_columns_factors]
#       }
#       B_F2 = B_F^2
#       B_F_tau[1,] = rgamma(b_F, shape = fixed_factors_prec_shape + ncol(B_F2)/2, rate = fixed_factors_prec_rate + rowSums(B_F2 * B_F_prec/c(B_F_tau))/2)
#       # B_F_tau[1,non_QTL_factor] = rgamma(length(non_QTL_factor), shape = fixed_factors_prec_shape + ncol(B_F2)/2,
#       #                                 rate = fixed_factors_prec_rate + rowSums((B_F2[non_QTL_factors,,drop=FALSE] * B_F_prec[non_QTL_factors,,drop=FALSE]/c(B_F_tau[non_QTL_factors])))/2)
#       B_F_prec[non_QTL_factor,] = matrix(rgamma(length(non_QTL_factor)*k,shape = (B_df + 1)/2,rate = (B_df + B_F2[non_QTL_factor,,drop=FALSE]*c(B_F_tau[non_QTL_factor]))/2),nr = length(non_QTL_factor),nc = k)
#       B_F_prec[non_QTL_factor,] = B_F_prec[non_QTL_factor,,drop=FALSE]*c(B_F_tau[non_QTL_factor])
#
#       # inverse Gaussian Bayesian Elastic Net from Li and Lin 2010
#       B_F_tau[QTL_columns_factors] = 1
#       QTL_lambda1 = B_df
#       QTL_lambda2 = B_df
#       tau = matrix(1/rinvgauss(n=length(QTL_columns_factors)*k,
#                                mean = sqrt(QTL_lambda1)/(2*QTL_lambda2*abs(B_F[QTL_columns_factors,,drop=FALSE])),
#                                shape = rep(QTL_lambda1*B_F_tau[QTL_columns_factors]/(4*QTL_lambda2),k)
#       )+1,
#       nr = length(QTL_columns_factors),nc = k)
#       B_F_prec[QTL_columns_factors,] = tau/(tau-1) * B_F_tau[QTL_columns_factors] * QTL_lambda2
#       B_F_prec[B_F_prec<0] = 1e-10
#       if(min(B_F_prec) < 0) recover()
#       rm(tau)
#
#       B_F_tau[1,X_F_zero_variance] = 1e10
#       B_F_prec[X_F_zero_variance,] = 1e10
#     }
#   }))
#   return(current_state)
# }
