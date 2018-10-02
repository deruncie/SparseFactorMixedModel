#' Run BSFG Gibbs sampler
#'
#' Run MCMC chain for a specified number of iterations
#'
#' @param BSFG_state BSFG_state object of current chain
#' @param n_samples Number of iterations to add to the chain (not number of posterior samples to draw.
#'     This is determined by n_samples / thin)
#' @param grainSize Minimum size of sub-problems for dividing among processes. Sent to RcppParallel
sample_BSFG = function(BSFG_state,n_samples,grainSize = 1,verbose=TRUE,...) {
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

  if(verbose) pb = txtProgressBar(min=start_i,max = start_i+n_samples,style=3)
  start_time = Sys.time()
  for(i in start_i+(1:n_samples)){
    BSFG_state$current_state$nrun = i
    BSFG_state$current_state = BSFG_state$current_state[!sapply(BSFG_state$current_state,is.null)]

    # ----- Sample model parameters  except precisions ---------------- #
    BSFG_state$current_state = sample_latent_traits(BSFG_state,grainSize = grainSize,...)

    # -----Sample Lambda_prec ------------- #
    BSFG_state$current_state = BSFG_state$priors$Lambda_prior$sampler(BSFG_state,...)

    # -----Sample B2_prec ------------- #
    BSFG_state$current_state = BSFG_state$priors$B2_prior$sampler(BSFG_state,...)

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
    if(verbose) setTxtProgressBar(pb, i)
  }
  end_time = Sys.time()
  if(verbose) close(pb)
  print(end_time - start_time)
  BSFG_state$current_state$total_time = BSFG_state$current_state$total_time + end_time - start_time

  # ----------------------------------------------- #
  # ------------Save state for restart------------- #
  # ----------------------------------------------- #

  current_state = BSFG_state$current_state
  saveRDS(current_state,file=sprintf('%s/current_state.rds',BSFG_state$run_ID))

  BSFG_state$RNG = list(
    Random.seed = .Random.seed,
    RNGkind = RNGkind()
  )
  return(BSFG_state)
}
