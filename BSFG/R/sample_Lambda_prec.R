sample_Lambda_prec_ARD = function(BSFG_state,...) {
  priors         = BSFG_state$priors
  run_variables  = BSFG_state$run_variables
  run_parameters = BSFG_state$run_parameters
  current_state  = BSFG_state$current_state

  current_state = with(c(priors,run_variables,run_parameters),
                       with(list(
                         # load priors
                         Lambda_df     = Lambda_prior$Lambda_df,
                         delta_1_rate  = Lambda_prior$delta_1$rate,
                         delta_1_shape = Lambda_prior$delta_1$shape,
                         delta_2_rate  = Lambda_prior$delta_2$rate,
                         delta_2_shape = Lambda_prior$delta_2$shape
                       ),within(current_state,{

                         # initialize variables if needed
                         if(!exists('delta')){
                           delta = with(priors,matrix(c(rgamma(1,shape = delta_1_shape,rate = delta_1_rate),rgamma(k-1,shape = delta_2_shape,rate = delta_2_rate)),nrow=1))
                           tauh  = matrix(cumprod(delta),nrow=1)
                           Lambda_prec = matrix(1,p,k)
                           Plam = sweep(Lambda_prec,2,tauh,'*')
                           # Lambda[] = Lambda / sqrt(Plam)
                         }

                         Lambda2 = Lambda^2
                         Lambda_prec[] = matrix(rgamma(p*k,shape = (Lambda_df + 1)/2,rate = (Lambda_df + sweep(Lambda2,2,tauh,'*'))/2),nr = p,nc = k)

                         # # trait one is special?
                         # Lambda_prec[1,] = 1e-10

                         # # -----Sample delta, update tauh------ #
                         scores = 0.5*colSums(Lambda2*Lambda_prec)
                         shapes = c(delta_1_shape + 0.5*p*k,
                                    delta_2_shape + 0.5*p*((k-1):1))
                         times = delta_iteractions_factor
                         randg_draws = matrix(rgamma(times*k,shape = shapes,rate = 1),nr=times,byrow=T)
                         delta[] = sample_delta_c_Eigen( delta,tauh,scores,delta_1_rate,delta_2_rate,randg_draws)
                         tauh[]  = matrix(cumprod(delta),nrow=1)

                         # # -----Update Plam-------------------- #
                         Plam[] = sweep(Lambda_prec,2,tauh,'*')
                         if('Plam_filter' %in% ls()) Plam[] = Plam * Plam_filter
                       })))
  return(current_state)
}

sample_Lambda_prec_ARD_v2 = function(BSFG_state,...) {
  # hierarchical delta/tauh
  priors         = BSFG_state$priors
  run_variables  = BSFG_state$run_variables
  run_parameters = BSFG_state$run_parameters
  current_state  = BSFG_state$current_state

  current_state = with(c(priors,run_variables,run_parameters),
                       with(list(
                         # load priors
                         Lambda_df     = Lambda_prior$Lambda_df,
                         delta_1_rate  = Lambda_prior$delta_1$rate,
                         delta_1_shape = Lambda_prior$delta_1$shape,
                         delta_2_rate  = Lambda_prior$delta_2$rate,
                         delta_2_shape = Lambda_prior$delta_2$shape
                       ),within(current_state,{

                         # initialize variables if needed
                         if(!exists('delta')){
                           delta = with(priors,matrix(c(rgamma(1,shape = delta_1_shape,rate = delta_1_rate),rgamma(k-1,shape = delta_2_shape,rate = delta_2_rate)),nrow=1))
                           tauh  = matrix(cumprod(delta),nrow=1)
                           Lambda_prec = Plam = matrix(1,p,k)
                         }

                         Lambda_prec[] = matrix(rgamma(p*k,shape = (Lambda_df + 1)/2,rate = tauh[rep(1,p),] + Lambda^2/2),nr = p,nc = k)

                         # # trait one is special?
                         # Lambda_prec[1,] = 1e-10

                         # # -----Sample delta, update tauh------ #
                         shapes = c(delta_1_shape + Lambda_df*0.5*p*k,
                                    delta_2_shape + Lambda_df*0.5*p*((k-1):1))
                         times = delta_iteractions_factor
                         randg_draws = matrix(rgamma(times*k,shape = shapes,rate = 1),nr=times,byrow=T)
                         scores = colSums(Lambda_prec)
                         delta[] = sample_delta_c_Eigen( delta,tauh,scores,delta_1_rate,delta_2_rate,randg_draws)
                         tauh[]  = matrix(cumprod(delta),nrow=1)

                         # # -----Update Plam-------------------- #
                         Plam[] = Lambda_prec
                       })))
  return(current_state)
}

sample_Lambda_prec_ARD_group = function(BSFG_state,...) {
  priors         = BSFG_state$priors
  run_variables  = BSFG_state$run_variables
  run_parameters = BSFG_state$run_parameters
  current_state  = BSFG_state$current_state

  current_state = with(c(priors,run_variables,run_parameters),
                       with(list(
                         # load priors
                         Lambda_df     = Lambda_prior$Lambda_df,
                         delta_1_rate  = Lambda_prior$delta_1$rate,
                         delta_1_shape = Lambda_prior$delta_1$shape,
                         delta_2_rate  = Lambda_prior$delta_2$rate,
                         delta_2_shape = Lambda_prior$delta_2$shape,
                         prior_Lambda_group_shape = Lambda_prior$prior_Lambda_group_shape,
                         prior_Lambda_group_rate  = Lambda_prior$prior_Lambda_group_rate,
                         group_factor = Lambda_prior$group_factor
                       ),within(current_state,{

                         # initialize variables if needed
                         if(!exists('delta')){
                           delta = with(priors,matrix(c(rgamma(1,shape = delta_1_shape,rate = delta_1_rate),rgamma(k-1,shape = delta_2_shape,rate = delta_2_rate)),nrow=1))
                           tauh  = matrix(cumprod(delta),nrow=1)
                           Lambda_prec = Plam = matrix(1,p,k)
                         }

                         if(!exists(group_shrinkage)){
                           Lambda_groups = match(levels())
                         }

                         Lambda2 = Lambda^2
                         shrinkage = tauh[1:nrow(Lambda),] * group_shrinkage
                         Lambda_prec[] = matrix(rgamma(p*k,shape = (Lambda_df + 1)/2,rate = (Lambda_df + sweep(Lambda2,2,tauh,'*'))/2),nr = p,nc = k)

                         # # trait one is special?
                         # Lambda_prec[1,] = 1e-10

                         # # -----Sample delta, update tauh------ #
                         shapes = c(delta_1_shape + 0.5*p*k,
                                    delta_2_shape + 0.5*p*((k-1):1))
                         times = delta_iteractions_factor
                         randg_draws = matrix(rgamma(times*k,shape = shapes,rate = 1),nr=times,byrow=T)
                         scores = 0.5*colSums(Lambda2*Lambda_prec)
                         delta[] = sample_delta_c_Eigen( delta,tauh,scores,delta_1_rate,delta_2_rate,randg_draws)
                         tauh[]  = matrix(cumprod(delta),nrow=1)

                         # # -----Update Plam-------------------- #
                         Plam[] = sweep(Lambda_prec,2,tauh,'*')
                         if('Plam_filter' %in% ls()) Plam[] = Plam * Plam_filter
                       })))
  return(current_state)
}


sample_Lambda_prec_TPB = function(BSFG_state,ncores = detectCores(),cluster = NULL,...) {
  priors         = BSFG_state$priors
  run_variables  = BSFG_state$run_variables
  run_parameters = BSFG_state$run_parameters
  current_state  = BSFG_state$current_state

  current_state = with(c(priors,run_variables,run_parameters),
                       with(list(
                         # load priors
                         Lambda_A      = Lambda_prior$Lambda_A,
                         Lambda_B      = Lambda_prior$Lambda_B,
                         delta_1_rate  = Lambda_prior$delta_1$rate,
                         delta_1_shape = Lambda_prior$delta_1$shape,
                         delta_2_rate  = Lambda_prior$delta_2$rate,
                         delta_2_shape = Lambda_prior$delta_2$shape
                       ),within(current_state,{

                         # initialize variables if needed
                         if(!exists('delta')){
                           delta = with(priors,matrix(c(rgamma(1,shape = delta_1_shape,rate = delta_1_rate),rgamma(k-1,shape = delta_2_shape,rate = delta_2_rate)),nrow=1))
                           tauh  = matrix(cumprod(delta),nrow=1)
                           Lambda_psi = matrix(rgamma(k,shape = 1/2,rate = 1),nrow=1)
                           Lambda_lambda = matrix(rgamma(p*k,shape = Lambda_B,rate=Lambda_psi[rep(1,p),]),p,k)
                           Lambda_prec = 1/matrix(rgamma(p*k,shape = Lambda_A,rate = Lambda_lambda),p,k)
                           Plam = matrix(1,p,k)
                         }

                         Lambda2 = Lambda^2

                         Lambda_psi[1,] = 1#rgamma(k,shape=p*Lambda_B + 1/2, rate = colSums(Lambda_lambda) + 1)
                         Lambda_lambda[] = rgamma(p*k,shape = Lambda_A + Lambda_B, rate = 1/Lambda_prec + Lambda_psi[rep(1,p),])
                         Lambda_prec[] = 1/rgig_multiple(p*k,lambda = rep( Lambda_A-1/2,p*k),chi = c(Lambda2*tauh[rep(1,p),]),psi = 2*c(Lambda_lambda))
                         Lambda_prec[Lambda_prec < 1e-10] = 1e-10
                         Lambda_prec[Lambda_prec > 1e10] = 1e10

                         # # trait one is special?
                         # Lambda_prec[1,] = 1e-10

                         # # -----Sample delta, update tauh------ #
                         shapes = c(delta_1_shape + 0.5*p*k,
                                    delta_2_shape + 0.5*p*((k-1):1))
                         times = delta_iteractions_factor
                         randg_draws = matrix(rgamma(times*k,shape = shapes,rate = 1),nr=times,byrow=T)
                         scores = 0.5*colSums(Lambda2*Lambda_prec)
                         delta[] = sample_delta_c_Eigen( delta,tauh,scores,delta_1_rate,delta_2_rate,randg_draws)
                         tauh[]  = matrix(cumprod(delta),nrow=1)

                         # # -----Update Plam-------------------- #
                         Plam[] = Lambda_prec * tauh[rep(1,p),]
                       })))
  return(current_state)
}
