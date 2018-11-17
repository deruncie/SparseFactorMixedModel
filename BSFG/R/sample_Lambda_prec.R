sample_Lambda_prec_horseshoe = function(BSFG_state,...) {
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
                       with(Lambda_prior,{

                         if(!exists('delta_iterations_factor')) delta_iterations_factor = 100

                         # delta_1_shape = delta_1$shape  delta_1 = 1
                         # delta_1_rate  = delta_1$rate
                         delta_l_shape = delta_l$shape
                         delta_l_rate  = delta_l$rate

                         tau_0 = prop_0/(1-prop_0) * 1/sqrt(n)

                         within(current_state,{

                           # initialize variables if needed
                           if(!'Lambda_tau2' %in% names(current_state)){
                             if(verbose) print('initializing Lambda_prec regularized horseshoe')
                             Lambda_tau2 = matrix(1,1,1)
                             Lambda_xi = matrix(1,1,1)
                             Lambda_phi2 = matrix(1,p,K)
                             Lambda_nu = matrix(1,p,K)
                             delta = with(priors,matrix(c(1,rgamma(K-1,shape = delta_l_shape,rate = delta_l_rate)),nrow=1))
                             Lambda_prec = matrix(1,p,K)
                             trunc_point_delta = 1
                             Lambda_m_eff = matrix(1,1,K)
                           }

                           Lambda2 = Lambda^2
                           Lambda2_std = Lambda2 * (tot_Eta_prec[1,]) / 2

                           Lambda_nu[] = matrix(1/rgamma(p*K,shape = 1, rate = 1 + 1/Lambda_phi2), nr = p, nc = K)
                           Lambda2_std_delta = sweep(Lambda2_std,2, cumprod(delta),'*')
                           Lambda_phi2[] = matrix(1/rgamma(p*K,shape = 1, rate = 1/Lambda_nu + Lambda2_std_delta / Lambda_tau2[1]),nr=p,nc = K)

                           scores = colSums(Lambda2_std / Lambda_phi2)
                           # for(i in 1:delta_iterations_factor) {
                           #   cumprod_delta = cumprod(delta[1,])
                           #   Lambda_tau2[] = 1/rgamma(1,shape = (p*K+1)/2, rate = 1/Lambda_xi[1] + sum(cumprod_delta*scores))
                           #   Lambda_xi[] = 1/rgamma(1,shape = 1,rate = 1/tau_0^2 + 1/Lambda_tau2[1])
                           #   for(h in 2:K) {
                           #     delta[h] = rgamma(1,shape = delta_l_shape + p/2*(K-h+1),rate = delta_l_rate + sum(cumprod_delta[h:K]*scores[h:K])/(Lambda_tau2[1]*delta[h]))
                           #     cumprod_delta = cumprod(delta[1,])
                           #   }
                           # }
                           new_samples = sample_tau2_delta_c_Eigen_v2(Lambda_tau2[1],Lambda_xi[1],delta,scores,
                                                                      tau_0,delta_l_shape,delta_l_rate,
                                                                      p,delta_iterations_factor)
                           Lambda_tau2[] = new_samples$tau2
                           Lambda_xi[] = new_samples$xi
                           delta[] = new_samples$delta

                           # -----Update Plam-------------------- #
                           Lambda_prec[] = 1/(Lambda_tau2[1] * sweep(Lambda_phi2,2,cumprod(delta),'/'))

                           # ----- Calcualte m_eff -------------- #
                           kappa = 1/(1+n/Lambda_prec)
                           Lambda_m_eff[] = colSums(1-kappa)

                           rm(list = c('Lambda2','Lambda2_std','Lambda2_std_delta','scores','new_samples','kappa'))
                         })
                       }))
  return(current_state)
}


sample_Lambda_prec_ARD = function(BSFG_state,...) {
  priors         = BSFG_state$priors
  run_variables  = BSFG_state$run_variables
  run_parameters = BSFG_state$run_parameters
  current_state  = BSFG_state$current_state

  current_state = with(c(priors,run_variables,run_parameters),
                       with(Lambda_prior,{

                         if(!exists('delta_iteractions_factor')) delta_iteractions_factor = 100

                         delta_1_shape = delta_1$shape
                         delta_1_rate  = delta_1$rate
                         delta_2_shape = delta_2$shape
                         delta_2_rate  = delta_2$rate

                         within(current_state,{

                         # initialize variables if needed
                         if(!exists('delta')){
                           delta = with(priors,matrix(c(rgamma(1,shape = delta_1_shape,rate = delta_1_rate),rgamma(K-1,shape = delta_2_shape,rate = delta_2_rate)),nrow=1))
                           # tauh  = matrix(cumprod(delta),nrow=1)
                           Lambda_phi = Lambda_prec = matrix(1,p,K)
                           # Plam = sweep(Lambda_prec,2,tauh,'*')
                           # Lambda[] = Lambda / sqrt(Plam)
                           trunc_point_delta = 1
                         }

                         Lambda2 = Lambda^2
                         Lambda2_std = Lambda2 * (tot_Eta_prec[1,]) #/ 2
                         tauh = cumprod(delta)

                         # Lambda_phi[] = matrix(rgamma(p*K,shape = (Lambda_df + 1)/2,rate = (Lambda_df + sweep(Lambda2,2,tauh,'*'))/2),nr = p,nc = K)
                         Lambda_phi[] = matrix(rgamma(p*K,shape = (Lambda_df + 1)/2,rate = (Lambda_df + sweep(Lambda2_std,2,tauh,'*'))/2),nr = p,nc = K)
                         # if(lambda_propto_Vp){
                         #   Lambda2 = Lambda2 * tot_Eta_prec[1,]
                         # }
                         # Lambda_prec[] = matrix(rgamma(p*K,shape = (Lambda_df + 1)/2,rate = (Lambda_df + sweep(Lambda2,2,tauh,'*'))/2),nr = p,nc = K)

                         # # trait one is special?
                         # Lambda_prec[1,] = 1e-10

                         # # -----Sample delta, update tauh------ #
                         scores = 0.5*colSums(Lambda2_std*Lambda_phi)
                         shapes = c(delta_1_shape + 0.5*p*K,
                                    delta_2_shape + 0.5*p*((K-1):1))
                         times = delta_iteractions_factor
                         # randg_draws = matrix(rgamma(times*K,shape = shapes,rate = 1),nr=times,byrow=T)
                         # delta[] = sample_delta_c_Eigen( delta,tauh,scores,delta_1_rate,delta_2_rate,randg_draws)
                         randu_draws = matrix(runif(times*K),nr=times)
                         delta[] = sample_trunc_delta_c_Eigen( delta,tauh,scores,shapes,delta_1_rate,delta_2_rate,randu_draws,trunc_point_delta)
                         tauh[]  = matrix(cumprod(delta),nrow=1)

                         Lambda_prec[] = sweep(Lambda_phi,2,tauh,'*')
                         # # # -----Update Plam-------------------- #
                         # Plam[] = sweep(Lambda_prec,2,tauh,'*')
                         # if(lambda_propto_Vp){
                         #  Plam[] = Plam * tot_Eta_prec[1,]
                         # }
                       })
                }))
  return(current_state)
}

# sample_Lambda_prec_ARD_v2 = function(BSFG_state,...) {
#   # hierarchical delta/tauh
#   priors         = BSFG_state$priors
#   run_variables  = BSFG_state$run_variables
#   run_parameters = BSFG_state$run_parameters
#   current_state  = BSFG_state$current_state
#
#   current_state = with(c(priors,run_variables,run_parameters),
#                        with(list(
#                          # load priors
#                          Lambda_df     = Lambda_prior$Lambda_df,
#                          delta_1_rate  = Lambda_prior$delta_1$rate,
#                          delta_1_shape = Lambda_prior$delta_1$shape,
#                          delta_2_rate  = Lambda_prior$delta_2$rate,
#                          delta_2_shape = Lambda_prior$delta_2$shape
#                        ),within(current_state,{
#
#                          # initialize variables if needed
#                          if(!exists('delta')){
#                            delta = with(priors,matrix(c(rgamma(1,shape = delta_1_shape,rate = delta_1_rate),rgamma(K-1,shape = delta_2_shape,rate = delta_2_rate)),nrow=1))
#                            tauh  = matrix(cumprod(delta),nrow=1)
#                            Lambda_prec = Plam = matrix(1,p,K)
#                          }
#
#                          Lambda_prec[] = matrix(rgamma(p*K,shape = (Lambda_df + 1)/2,rate = tauh[rep(1,p),] + Lambda^2/2),nr = p,nc = K)
#
#                          # # trait one is special?
#                          # Lambda_prec[1,] = 1e-10
#
#                          # # -----Sample delta, update tauh------ #
#                          shapes = c(delta_1_shape + Lambda_df*0.5*p*K,
#                                     delta_2_shape + Lambda_df*0.5*p*((K-1):1))
#                          times = delta_iteractions_factor
#                          randg_draws = matrix(rgamma(times*K,shape = shapes,rate = 1),nr=times,byrow=T)
#                          scores = colSums(Lambda_prec)
#                          delta[] = sample_delta_c_Eigen( delta,tauh,scores,delta_1_rate,delta_2_rate,randg_draws)
#                          tauh[]  = matrix(cumprod(delta),nrow=1)
#
#                          # # -----Update Plam-------------------- #
#                          Plam[] = Lambda_prec
#                        })))
#   return(current_state)
# }
#
# sample_Lambda_prec_ARD_group = function(BSFG_state,...) {
#   priors         = BSFG_state$priors
#   run_variables  = BSFG_state$run_variables
#   run_parameters = BSFG_state$run_parameters
#   current_state  = BSFG_state$current_state
#
#   current_state = with(c(priors,run_variables,run_parameters),
#                        with(list(
#                          # load priors
#                          Lambda_df     = Lambda_prior$Lambda_df,
#                          delta_1_rate  = Lambda_prior$delta_1$rate,
#                          delta_1_shape = Lambda_prior$delta_1$shape,
#                          delta_2_rate  = Lambda_prior$delta_2$rate,
#                          delta_2_shape = Lambda_prior$delta_2$shape,
#                          prior_Lambda_group_shape = Lambda_prior$prior_Lambda_group_shape,
#                          prior_Lambda_group_rate  = Lambda_prior$prior_Lambda_group_rate,
#                          group_factor = Lambda_prior$group_factor
#                        ),within(current_state,{
#
#                          # initialize variables if needed
#                          if(!exists('delta')){
#                            delta = with(priors,matrix(c(rgamma(1,shape = delta_1_shape,rate = delta_1_rate),rgamma(K-1,shape = delta_2_shape,rate = delta_2_rate)),nrow=1))
#                            tauh  = matrix(cumprod(delta),nrow=1)
#                            Lambda_prec = Plam = matrix(1,p,K)
#                          }
#
#                          if(!exists(group_shrinkage)){
#                            Lambda_groups = match(levels())
#                          }
#
#                          Lambda2 = Lambda^2
#                          shrinkage = tauh[1:nrow(Lambda),] * group_shrinkage
#                          Lambda_prec[] = matrix(rgamma(p*K,shape = (Lambda_df + 1)/2,rate = (Lambda_df + sweep(Lambda2,2,tauh,'*'))/2),nr = p,nc = K)
#
#                          # # trait one is special?
#                          # Lambda_prec[1,] = 1e-10
#
#                          # # -----Sample delta, update tauh------ #
#                          shapes = c(delta_1_shape + 0.5*p*K,
#                                     delta_2_shape + 0.5*p*((K-1):1))
#                          times = delta_iteractions_factor
#                          randg_draws = matrix(rgamma(times*K,shape = shapes,rate = 1),nr=times,byrow=T)
#                          scores = 0.5*colSums(Lambda2*Lambda_prec)
#                          delta[] = sample_delta_c_Eigen( delta,tauh,scores,delta_1_rate,delta_2_rate,randg_draws)
#                          tauh[]  = matrix(cumprod(delta),nrow=1)
#
#                          # # -----Update Plam-------------------- #
#                          Plam[] = sweep(Lambda_prec,2,tauh,'*')
#                          if('Plam_filter' %in% ls()) Plam[] = Plam * Plam_filter
#                        })))
#   return(current_state)
# }
#
#
# sample_Lambda_prec_TPB = function(BSFG_state,ncores = detectCores(),cluster = NULL,...) {
#   priors         = BSFG_state$priors
#   run_variables  = BSFG_state$run_variables
#   run_parameters = BSFG_state$run_parameters
#   current_state  = BSFG_state$current_state
#
#   current_state = with(c(priors,run_variables,run_parameters),
#                        with(list(
#                          # load priors
#                          Lambda_A      = Lambda_prior$Lambda_A,
#                          Lambda_B      = Lambda_prior$Lambda_B,
#                          delta_1_rate  = Lambda_prior$delta_1$rate,
#                          delta_1_shape = Lambda_prior$delta_1$shape,
#                          delta_2_rate  = Lambda_prior$delta_2$rate,
#                          delta_2_shape = Lambda_prior$delta_2$shape
#                        ),within(current_state,{
#
#                          # initialize variables if needed
#                          if(!exists('delta')){
#                            delta = with(priors,matrix(c(rgamma(1,shape = delta_1_shape,rate = delta_1_rate),rgamma(K-1,shape = delta_2_shape,rate = delta_2_rate)),nrow=1))
#                            tauh  = matrix(cumprod(delta),nrow=1)
#                            Lambda_psi = matrix(rgamma(K,shape = 1/2,rate = 1),nrow=1)
#                            Lambda_lambda = matrix(rgamma(p*K,shape = Lambda_B,rate=Lambda_psi[rep(1,p),]),p,K)
#                            Lambda_prec = 1/matrix(rgamma(p*K,shape = Lambda_A,rate = Lambda_lambda),p,K)
#                            Plam = matrix(1,p,K)
#                          }
#
#                          Lambda2 = Lambda^2
#
#                          Lambda_psi[1,] = 1#rgamma(K,shape=p*Lambda_B + 1/2, rate = colSums(Lambda_lambda) + 1)
#                          Lambda_lambda[] = rgamma(p*K,shape = Lambda_A + Lambda_B, rate = 1/Lambda_prec + Lambda_psi[rep(1,p),])
#                          Lambda_prec[] = 1/rgig_multiple(p*K,lambda = rep( Lambda_A-1/2,p*K),chi = c(Lambda2*tauh[rep(1,p),]),psi = 2*c(Lambda_lambda))
#                          Lambda_prec[Lambda_prec < 1e-10] = 1e-10
#                          Lambda_prec[Lambda_prec > 1e10] = 1e10
#
#                          # # trait one is special?
#                          # Lambda_prec[1,] = 1e-10
#
#                          # # -----Sample delta, update tauh------ #
#                          shapes = c(delta_1_shape + 0.5*p*K,
#                                     delta_2_shape + 0.5*p*((K-1):1))
#                          times = delta_iteractions_factor
#                          randg_draws = matrix(rgamma(times*K,shape = shapes,rate = 1),nr=times,byrow=T)
#                          scores = 0.5*colSums(Lambda2*Lambda_prec)
#                          delta[] = sample_delta_c_Eigen( delta,tauh,scores,delta_1_rate,delta_2_rate,randg_draws)
#                          tauh[]  = matrix(cumprod(delta),nrow=1)
#
#                          # # -----Update Plam-------------------- #
#                          Plam[] = Lambda_prec * tauh[rep(1,p),]
#                        })))
#   return(current_state)
# }
#
# sample_Lambda_prec_horseshoe = function(BSFG_state,...) {
#   # using algorithm from Makalic and Schmidt (2015)
#   # adapting following guidance from Piironen and Vehtari for prior on tau2 (omega2)
#   # for now, will remove omega2 and just use delta/tauh as before
#   # idea is to calibrate prior on tau2 based on expectation of meff (proportion of non-zero entries)
#   # Lambda_omega2/tauh | tauh ~ C+(0,1/tauh) controls the sparsity of the factors
#   # delta_1 controls the proportion in the 1st factor. This number should be approximately fixed.
#   # delta_i controls the decrease in odds of inclusion of each element for the ith factor relative to the (i-1)th
#   # Note:  Piironen and Vehtari do not include sigma^2 in prior, so I have removed
#   priors         = BSFG_state$priors
#   run_variables  = BSFG_state$run_variables
#   run_parameters = BSFG_state$run_parameters
#   current_state  = BSFG_state$current_state
#
#   current_state = with(c(priors,run_variables,run_parameters),
#                        with(Lambda_prior,{
#
#                          if(!exists('delta_iteractions_factor')) delta_iteractions_factor = 100
#
#                          delta_1_shape = delta_1$shape
#                          delta_1_rate  = delta_1$rate
#                          delta_2_shape = delta_2$shape
#                          delta_2_rate  = delta_2$rate
#
#                          within(current_state,{
#
#                            # initialize variables if needed
#                            if(!exists('delta')){
#                              Lambda_omega2 = matrix(1,1,1)
#                              Lambda_xi = matrix(1,1,1)
#                              Lambda_phi2 = matrix(1,p,K)
#                              Lambda_nu = matrix(1,p,K)
#                              delta = with(priors,matrix(c(rgamma(1,shape = delta_1_shape,rate = delta_1_rate),rgamma(K-1,shape = delta_2_shape,rate = delta_2_rate)),nrow=1))
#                              tauh  = matrix(cumprod(delta),nrow=1)
#                              Plam = matrix(1,p,K)
#                              Lambda_prec = matrix(1,p,K)
#                              trunc_point_delta = 1
#                            }
#
#
#                            if(lambda_propto_Vp){
#                              Lambda2_std = Lambda^2 * (tot_Eta_prec[1,])
#                            } else{
#                              Lambda2_std = Lambda^2
#                            }
#
#                            Lambda_phi2[] = matrix(1/rgamma(p*K,shape = 1, rate = 1/Lambda_nu + Lambda2_std * (tauh/Lambda_omega2[1]/2)[rep(1,p),]),nr=p,nc = K)
#                            Lambda_nu[] = matrix(1/rgamma(p*K,shape = 1, rate = 1 + 1/Lambda_phi2), nr = p, nc = K)
#
#                            # -----Sample delta, update tauh------ #
#                            scores = colSums(Lambda2_std / Lambda_phi2)/2
#                            shapes = c((p*K + 1) / 2,1,
#                                       delta_1_shape + 0.5*p*K,
#                                       delta_2_shape + 0.5*p*((K-1):1))
#                            times = delta_iteractions_factor
#                            # randg_draws = matrix(rgamma(times*(2+K),shape = shapes,rate = 1),nr=times,byrow=T)
#                            # deltas = sample_delta_omega_c_Eigen(delta,tauh,Lambda_omega2,Lambda_xi,scores,delta_1_rate,delta_2_rate,randg_draws)
#                            randu_draws = matrix(runif(times*(2+K)),nr=times)
#                            deltas = sample_trunc_delta_omega_c_Eigen( delta,tauh,Lambda_omega2,Lambda_xi,scores,shapes,delta_1_rate,delta_2_rate,randu_draws,trunc_point_delta)
#                            Lambda_omega2[] = deltas[[1]]
#                            Lambda_xi[] = deltas[[2]]
#                            delta[] = deltas[[3]]
#                            tauh[]  = matrix(cumprod(delta),nrow=1)
#
#
#                            # -----Update Plam-------------------- #
#                            Lambda_prec[] = (1/Lambda_omega2[1])/Lambda_phi2
#                            Plam[] = Lambda_prec * tauh[rep(1,p),]
#                            if(lambda_propto_Vp){
#                              Plam[] = Plam * tot_Eta_prec[1,]
#                            }
#                          })
#                        }))
#   return(current_state)
# }
#
