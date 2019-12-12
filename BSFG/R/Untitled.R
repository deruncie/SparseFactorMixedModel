sample_B2_prec_horseshoe2 = function(BSFG_state,...) {
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

                         within(current_state,{

                           # initialize variables if needed
                           if(!any(c('B2_R_tau2','B2_F_tau2') %in% names(current_state))){
                             if(verbose) print('initializing B_prec regularized horseshoe')
                             B2_R_xi = matrix(1/rgamma(K,shape=1/2,rate=1/tau_0^2),nr=1)
                             B2_R_tau2 = matrix(1/rgamma(K,shape = 1/2, rate = 1/B2_R_xi[1,]),nr=1)
                             B2_R_nu = matrix(1/rgamma(b2_R*K,shape = 1/2, rate = 1), nr = b2_R, nc = K)
                             B2_R_phi2 = matrix(1/rgamma(b2_R*K,shape = 1/2, rate = 1/B2_R_nu),nr=b2_R,nc = K)
                             B2_R_prec = 1 / sweep(B2_R_phi2,2,B2_R_tau2[1,],'*')
                             B2_R = rstdnorm_mat(b2_R,K)/sqrt(B2_R_prec)


                             B2_F_xi = matrix(1/rgamma(K,shape=1/2,rate=1/tau_0^2),nr=1)
                             B2_F_tau2 = matrix(1/rgamma(K,shape = 1/2, rate = 1/B2_F_xi[1,]),nr=1)
                             B2_F_nu = matrix(1/rgamma(b2_F*K,shape = 1/2, rate = 1), nr = b2_F, nc = K)
                             B2_F_phi2 = matrix(1/rgamma(b2_F*K,shape = 1/2, rate = 1/B2_F_nu),nr=b2_F,nc = K)
                             B2_F_prec = 1 / sweep(B2_F_phi2,2,B2_F_tau2[1,],'*')
                             B2_F = rstdnorm_mat(b2_F,K)/sqrt(B2_F_prec)
                           } else {

                             B2_R_2 = B2_R^2
                             B2_R_2_std = sweep(B2_R_2,2,tot_Eta_prec[1,],'*')
                             B2_F_2 = B2_F^2
                             B2_F_2_std = sweep(B2_F_2,2,tot_F_prec[1,],'*')

                             B2_R_nu[] = matrix(1/rgamma(b2_R*p,shape = 1, rate = 1 + 1/B2_R_phi2), nr = b2_R, nc = p)
                             B2_R_phi2[] = matrix(1/rgamma(b2_R*p,shape = 1, rate = 1/B2_R_nu + sweep(B2_R_2_std,2,(2*B2_R_tau2),'/')),nr=b2_R,nc = p)
                             B2_F_nu[] = matrix(1/rgamma(b2_F*K,shape = 1, rate = 1 + 1/B2_F_phi2), nr = b2_F, nc = K)
                             B2_F_phi2[] = matrix(1/rgamma(b2_F*K,shape = 1, rate = 1/B2_F_nu + sweep(B2_F_2_std,2,(2*B2_F_tau2),'/')),nr=b2_F,nc = K)

                             B2_R_xi[] = 1/rgamma(p,shape=1,rate=1/tau_0^2 + 1/B2_R_tau2[1,])
                             B2_R_tau2[] = 1/rgamma(p,shape = (b2_R + 1)/2, rate = 1/B2_R_xi[1,] + (colSums(B2_R_2_std/B2_R_phi2))/2)

                             B2_F_xi[] = 1/rgamma(K,shape=1,rate=1/tau_0^2 + 1/B2_F_tau2[1,])
                             B2_F_tau2[] = 1/rgamma(K,shape = (b2_F + 1)/2, rate = 1/B2_F_xi[1,] + (colSums(B2_F_2_std/B2_F_phi2))/2)

                             rm(list=c('B2_R_2','B2_R_2_std','B2_F_2','B2_F_2_std'))
                           }


                           # -----Update Plam-------------------- #
                           B2_R_prec = 1 / sweep(B2_R_phi2,2,B2_R_tau2[1,],'*')
                           B2_F_prec = 1 / sweep(B2_F_phi2,2,B2_F_tau2[1,],'*')

                         })
                       }))
  return(current_state)
}
