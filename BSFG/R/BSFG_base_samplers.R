#' Run BSFG Gibbs sampler
#'
#' Run MCMC chain for a specified number of iterations
#'
#' @param BSFG_state BSFG_state object of current chain
#' @param n_samples Number of iterations to add to the chain (not number of posterior samples to draw.
#'     This is determined by n_samples / thin)
#' @param grainSize Minimum size of sub-problems for dividing among processes. Sent to RcppParallel
sample_BSFG = function(BSFG_state,n_samples,grainSize = 1,...) {
  data_matrices  = BSFG_state$data_matrices
  priors         = BSFG_state$priors
  run_parameters = BSFG_state$run_parameters
  run_variables  = BSFG_state$run_variables

  # ----------------------------------------------- #
  # -----------Reset Global Random Number Stream--- #
  # ----------------------------------------------- #
  do.call("RNGkind",as.list(BSFG_state$RNG$RNGkind))  ## must be first!
  assign(".Random.seed", BSFG_state$RNG$Random.seed, .GlobalEnv)

  # ----------------------------------------------- #
  # ----------------Set up run--------------------- #
  # ----------------------------------------------- #
  save_freq    = run_parameters$save_freq
  burn         = run_parameters$burn
  thin         = run_parameters$thin
  start_i      = BSFG_state$current_state$nrun

  # ----------------------------------------------- #
  # ---Extend posterior matrices for new samples--- #
  # ----------------------------------------------- #

  sp = (start_i + n_samples - burn)/thin - BSFG_state$Posterior$total_samples
  BSFG_state$Posterior = expand_Posterior(BSFG_state$Posterior,max(0,sp))

  # ----------------------------------------------- #
  # --------------start gibbs sampling------------- #
  # ----------------------------------------------- #

  start_time = Sys.time()
  for(i in start_i+(1:n_samples)){
    BSFG_state$current_state$nrun = i
    BSFG_state$current_state = BSFG_state$current_state[!sapply(BSFG_state$current_state,is.null)]

    # ----- Sample Lambda ---------------- #
    BSFG_state$current_state = sample_Lambda_B(BSFG_state,grainSize = grainSize,...)

    # BSFG_state$current_state$Lambda[1,2:min(5,BSFG_state$current_state$k)] = 1  # TEMPORARY!!!

    # ----- Sample other factor model parameters  ---------------- #
    BSFG_state$current_state = sample_latent_traits(BSFG_state,grainSize = grainSize,...)

    # -----Sample Lambda_prec ------------- #
    BSFG_state$current_state = sample_Lambda_prec(BSFG_state)

    # -----Sample prec_B ------------- #
    BSFG_state$current_state = sample_prec_B_ARD(BSFG_state)
    # BSFG_state$current_state = sample_prec_B_QTLBEN(BSFG_state)

    # ----- sample Eta ----- #
    observation_model_state = run_parameters$observation_model(run_parameters$observation_model_parameters,BSFG_state)$state
    BSFG_state$current_state[names(observation_model_state)] = observation_model_state

    # -- adapt number of factors to samples ---#
    # if(i > 200 && i < burn && runif(1) < with(BSFG_state$run_parameters,1/exp(b0 + b1*i))){  # adapt with decreasing probability per iteration
    #   BSFG_state$current_state = update_k(BSFG_state)
    # }

    # -- save sampled values (after thinning) -- #
    if( (i-burn) %% thin == 0 && i > burn) {
      BSFG_state$Posterior = save_posterior_sample(BSFG_state)
    }
  }
  end_time = Sys.time()
  print(end_time - start_time)
  BSFG_state$current_state$total_time = BSFG_state$current_state$total_time + end_time - start_time

  # ----------------------------------------------- #
  # ------------Save state for restart------------- #
  # ----------------------------------------------- #

  current_state = BSFG_state$current_state
  save(current_state,file='current_state.RData')

  BSFG_state$RNG = list(
    Random.seed = .Random.seed,
    RNGkind = RNGkind()
  )
  return(BSFG_state)
}

sample_Lambda_B = function(BSFG_state,...){
  UseMethod("sample_Lambda_B",BSFG_state)
}

sample_latent_traits = function(BSFG_state,...){
  UseMethod("sample_latent_traits",BSFG_state)
}

sample_Lambda_prec = function(BSFG_state) {
  priors         = BSFG_state$priors
  run_variables  = BSFG_state$run_variables
  run_parameters = BSFG_state$run_parameters
  current_state  = BSFG_state$current_state

  current_state = with(c(priors,run_variables,run_parameters),within(current_state,{
    Lambda2 = Lambda^2
    Lambda_prec = matrix(rgamma(p*k,shape = (Lambda_df + 1)/2,rate = (Lambda_df + sweep(Lambda2,2,tauh,'*'))/2),nr = p,nc = k)

    # # trait one is special?
    # Lambda_prec[1,] = 1e-10

    # # -----Sample delta, update tauh------ #
    shapes = c(delta_1_shape + 0.5*p*k,
               delta_2_shape + 0.5*p*((k-1):1))
    times = delta_iteractions_factor
    randg_draws = matrix(rgamma(times*k,shape = shapes,rate = 1),nr=times,byrow=T)
    delta[] = sample_delta_c_Eigen( delta,tauh,Lambda_prec,delta_1_shape,delta_1_rate,delta_2_shape,delta_2_rate,randg_draws,Lambda2)
    tauh[]  = matrix(cumprod(delta),nrow=1)

    # # -----Update Plam-------------------- #
    Plam[] = sweep(Lambda_prec,2,tauh,'*')
  }))
  return(current_state)
}

sample_prec_B = function(BSFG_state){
  priors         = BSFG_state$priors
  run_variables  = BSFG_state$run_variables
  current_state  = BSFG_state$current_state

  current_state = with(c(priors,run_variables),within(current_state,{
    if(b > 0) {
      if(same_fixed_model) {  # want same tau for both
        B2 = cbind(B,B_F)^2  # tauh
      } else{
        B2 = B^2
      }
      tau_B[1,] = rgamma(b, shape = fixed_resid_prec_shape + ncol(B2)/2, rate = fixed_resid_prec_rate + rowSums(B2)/2)
      if(resid_intercept){
        tau_B[1,1] = 1e-10
      }
      prec_B = matrix(tau_B,nrow = b, ncol = p)
    }
    if(b_F > 0){
      if(same_fixed_model) {
        tau_B_F[1,] = tau_B[1,]
      } else{
        B_F2 = B_F^2
        tau_B_F[1,] = rgamma(b_F, shape = fixed_factors_prec_shape + ncol(B_F2)/2, rate = fixed_factors_prec_rate + rowSums(B_F2)/2)
      }
      tau_B_F[1,X_F_zero_variance] = 1e10
      prec_B_F = matrix(tau_B_F,nrow = b_F, ncol = k)
    }
  }))
  return(current_state)
}

sample_prec_B_ARD = function(BSFG_state){
  priors         = BSFG_state$priors
  run_variables  = BSFG_state$run_variables
  current_state  = BSFG_state$current_state

  current_state = with(c(priors,run_variables),within(current_state,{
    if(b > 0) {
      B2 = B^2
      tau_B[1,] = rgamma(b, shape = fixed_resid_prec_shape + ncol(B2)/2, rate = fixed_resid_prec_rate + rowSums((B2 * prec_B/c(tau_B)))/2)
      prec_B[] = matrix(rgamma(b*p,shape = (B_df + 1)/2,rate = (B_df + B2*c(tau_B))/2),nr = b,nc = p)
      prec_B[] = prec_B*c(tau_B)
      if(resid_intercept){
        tau_B[1,1] = 1e-10
        prec_B[1,] = 1e-10
      }
    }
    if(b_F > 0) {
      B_F2 = B_F^2
      tau_B_F[1,] = rgamma(b_F, shape = fixed_factors_prec_shape + ncol(B_F2)/2, rate = fixed_factors_prec_rate + rowSums(B_F2 * prec_B_F/c(tau_B_F))/2)
      prec_B_F[] = matrix(rgamma(b_F*k,shape = (B_F_df + 1)/2,rate = (B_F_df + B_F2*c(tau_B_F))/2),nr = b_F,nc = k)
      prec_B_F[] = prec_B_F*c(tau_B_F)
      tau_B_F[1,X_F_zero_variance] = 1e10
      prec_B_F[X_F_zero_variance,] = 1e10
    }
  }))
  return(current_state)
}

sample_prec_B_QTLBEN = function(BSFG_state){
  priors         = BSFG_state$priors
  run_variables  = BSFG_state$run_variables
  current_state  = BSFG_state$current_state

  QTL_columns_resid = BSFG_state$data_matrices$QTL_columns_resid
  QTL_columns_factors = BSFG_state$data_matrices$QTL_columns_factors

  current_state = with(c(priors,run_variables),within(current_state,{
    if(b > 0) {
      non_QTL_resid = 1:nrow(B)
      if(!is.null(QTL_columns_resid)){
        non_QTL_resid = non_QTL_resid[-QTL_columns_resid]
      }
      B2 = B^2
      tau_B[1,] = rgamma(b, shape = fixed_resid_prec_shape + ncol(B2)/2, rate = fixed_resid_prec_rate + rowSums((B2 * prec_B/c(tau_B)))/2)
      # tau_B[1,non_QTL_resid] = rgamma(length(non_QTL_resid), shape = fixed_resid_prec_shape + ncol(B2)/2,
      # rate = fixed_resid_prec_rate + rowSums((B2[non_QTL_resid,,drop=FALSE] * prec_B[non_QTL_resid,,drop=FALSE]/c(tau_B[non_QTL_resid])))/2)
      prec_B[non_QTL_resid,] = matrix(rgamma(length(non_QTL_resid)*p,shape = (B_df + 1)/2,rate = (B_df + B2[non_QTL_resid,,drop=FALSE]*c(tau_B[non_QTL_resid]))/2),nr = length(non_QTL_resid),nc = p)
      prec_B[non_QTL_resid,] = prec_B[non_QTL_resid,,drop=FALSE]*c(tau_B[non_QTL_resid])

      # inverse Gaussian Bayesian Elastic Net from Li and Lin 2010
      B_F_std = B_F[QTL_columns_resid,,drop=FALSE]
      QTL_lambda1 = B_df
      QTL_lambda2 = tau_B[QTL_columns_resid]
      tau = matrix(1/rinvgauss(n=length(QTL_columns_resid)*p,
                               mean = sqrt(QTL_lambda1)/(2*QTL_lambda2*abs(B_F_std)),
                               shape = QTL_lambda1/(4*QTL_lambda2))+1,
                   nr = length(QTL_columns_resid),nc = p)
      prec_B[QTL_columns_resid,] = tau/(tau-1) * QTL_lambda2

      if(resid_intercept){
        tau_B[1,1] = 1e-10
        prec_B[1,] = 1e-10
      }
    }
    if(b_F > 0) {
      non_QTL_factor = 1:nrow(B_F)
      if(!is.null(QTL_columns_factors)){
        non_QTL_factor = non_QTL_factor[-QTL_columns_factors]
      }
      B_F2 = B_F^2
      tau_B_F[1,] = rgamma(b_F, shape = fixed_factors_prec_shape + ncol(B_F2)/2, rate = fixed_factors_prec_rate + rowSums(B_F2 * prec_B_F/c(tau_B_F))/2)
      # tau_B_F[1,non_QTL_factor] = rgamma(length(non_QTL_factor), shape = fixed_factors_prec_shape + ncol(B_F2)/2,
      #                                 rate = fixed_factors_prec_rate + rowSums((B_F2[non_QTL_factors,,drop=FALSE] * prec_B_F[non_QTL_factors,,drop=FALSE]/c(tau_B_F[non_QTL_factors])))/2)
      prec_B_F[non_QTL_factor,] = matrix(rgamma(length(non_QTL_factor)*k,shape = (B_df + 1)/2,rate = (B_df + B_F2[non_QTL_factor,,drop=FALSE]*c(tau_B_F[non_QTL_factor]))/2),nr = length(non_QTL_factor),nc = k)
      prec_B_F[non_QTL_factor,] = prec_B_F[non_QTL_factor,,drop=FALSE]*c(tau_B_F[non_QTL_factor])

      # inverse Gaussian Bayesian Elastic Net from Li and Lin 2010
      tau_B_F[QTL_columns_factors] = 1
      QTL_lambda1 = B_df
      QTL_lambda2 = B_df
      tau = matrix(1/rinvgauss(n=length(QTL_columns_factors)*k,
                               mean = sqrt(QTL_lambda1)/(2*QTL_lambda2*abs(B_F[QTL_columns_factors,,drop=FALSE])),
                               shape = rep(QTL_lambda1*tau_B_F[QTL_columns_factors]/(4*QTL_lambda2),k)
      )+1,
      nr = length(QTL_columns_factors),nc = k)
      prec_B_F[QTL_columns_factors,] = tau/(tau-1) * tau_B_F[QTL_columns_factors] * QTL_lambda2
      prec_B_F[prec_B_F<0] = 1e-10
      if(min(prec_B_F) < 0) recover()
      rm(tau)

      tau_B_F[1,X_F_zero_variance] = 1e10
      prec_B_F[X_F_zero_variance,] = 1e10
    }
  }))
  return(current_state)
}
