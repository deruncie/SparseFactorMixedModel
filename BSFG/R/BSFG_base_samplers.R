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

    # ----- Sample Lambda ---------------- #
    BSFG_state$current_state = sample_Lambda_B(BSFG_state,grainSize = grainSize,...)

    # BSFG_state$current_state$Lambda[1,2:min(5,BSFG_state$current_state$k)] = 1  # TEMPORARY!!!

    # ----- Sample other factor model parameters  ---------------- #
    BSFG_state$current_state = sample_latent_traits(BSFG_state,grainSize = grainSize,...)

    # -----Sample Lambda_prec ------------- #
    BSFG_state$current_state = BSFG_state$priors$Lambda_prior$sampler(BSFG_state,...)

    # -----Sample B_prec ------------- #
    BSFG_state$current_state = BSFG_state$priors$B_prior$sampler(BSFG_state,...)
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
      Lambda_prec = Plam = matrix(1,p,k)
    }

    Lambda2 = Lambda^2
    Lambda_prec[] = matrix(rgamma(p*k,shape = (Lambda_df + 1)/2,rate = (Lambda_df + sweep(Lambda2,2,tauh,'*'))/2),nr = p,nc = k)

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

                         Lambda_psi[1,] = rgamma(k,shape=p*Lambda_B + 1/2, rate = colSums(Lambda_lambda) + 1)
                         Lambda_lambda[] = rgamma(p*k,shape = Lambda_A + Lambda_B, rate = 1/Lambda_prec + Lambda_psi[rep(1,p),])

                         if(!exists('Lambda_ncores')) Lambda_ncores = ncores
                         if(is.null(cluster[1])) {
                           Lambda_prec[] = 1/do.call(cbind,mclapply(1:k,function(j) {
                             sapply(1:p,function(i) GIGrvg::rgig(n=1,lambda = Lambda_A-1/2, chi = Lambda2[i,j]*tauh[j], psi = 2*Lambda_lambda[i,j]))
                           },mc.cores = Lambda_ncores))  # this replaces Lambda_tau - BSFG requires precision, not variance.
                         } else {
                           rm_cluster=FALSE
                           if(is.na(cluster[1])){
                             rm_cluster=TRUE
                             cluster = makeCluster(1)
                             print(cluster)
                           }
                           vars = c('k','p','Lambda_A','Lambda2','tauh','Lambda_lambda','Lambda_ncores')
                           clusterEnv = lapply(vars,function(x) eval(parse(text=x)))
                           names(clusterEnv) = vars
                           clusterEnv = as.environment(clusterEnv)
                           clusterExport(cluster,vars,envir = clusterEnv)
                           clusterEvalQ(cluster,library(parallel))
                           Lambda_prec[] = 1/do.call(cbind,parLapply(cluster,1,function(x) {
                             mclapply(1:k,function(j) {
                               sapply(1:p,function(i) GIGrvg::rgig(n=1,lambda = Lambda_A-1/2, chi = Lambda2[i,j]*tauh[j], psi = 2*Lambda_lambda[i,j]))
                             },mc.cores = Lambda_ncores)
                           })[[1]])
                           clusterEvalQ(cluster,gc)
                           if(rm_cluster){
                             stopCluster(cluster)
                             cluster=NA
                           }
                           rm(clusterEnv)
                           rm(vars)
                           gc()
                         }

                         Lambda_prec[Lambda_prec < 1e-10] = 1e-10
                         Lambda_prec[Lambda_prec > 1e10] = 1e10

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
                       })))
  return(current_state)
}

sample_B_prec = function(BSFG_state,...){
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
      B_tau[1,] = rgamma(b, shape = fixed_resid_prec_shape + ncol(B2)/2, rate = fixed_resid_prec_rate + rowSums(B2)/2)
      if(resid_intercept){
        B_tau[1,1] = 1e-10
      }
      B_prec = matrix(B_tau,nrow = b, ncol = p)
    }
    if(b_F > 0){
      if(same_fixed_model) {
        B_F_tau[1,] = B_tau[1,]
      } else{
        B_F2 = B_F^2
        B_F_tau[1,] = rgamma(b_F, shape = fixed_factors_prec_shape + ncol(B_F2)/2, rate = fixed_factors_prec_rate + rowSums(B_F2)/2)
      }
      B_F_tau[1,X_F_zero_variance] = 1e10
      B_F_prec = matrix(B_F_tau,nrow = b_F, ncol = k)
    }
  }))
  return(current_state)
}

sample_B_prec_ARD = function(BSFG_state,...){
  priors         = BSFG_state$priors
  run_variables  = BSFG_state$run_variables
  current_state  = BSFG_state$current_state

  current_state = with(c(priors,run_variables),
                       with(list(
                                  # load priors
                                  B_df   = B_prior$B_df,
                                  B_F_df = B_prior$B_F_df
                       ),within(current_state,{

    # initialize variables if needed
    if(!exists('B_tau')){
      if(b > 0) {
        B_tau = matrix(c(1e-10,rgamma(b-1,shape = priors$fixed_resid_prec_shape, rate = priors$fixed_resid_prec_rate)),nrow=1)
      } else{
        B_tau = matrix(0,nrow=1,ncol=0)
      }
      if(b_F > 0) {
        B_F_tau = matrix(rgamma(b_F,shape = priors$fixed_factors_prec_shape, rate = priors$fixed_factors_prec_rate),nrow=1)
      } else{
        B_F_tau = matrix(0,nrow=1,ncol=0)
      }

      B_prec = matrix(B_tau,nrow = b, ncol = p)
      B_F_prec = matrix(B_F_tau,nrow = b_F, ncol = k)
    }
    if(b > 0) {
      B2 = B^2
      B_tau[1,] = rgamma(b, shape = fixed_resid_prec_shape + ncol(B2)/2, rate = fixed_resid_prec_rate + rowSums((B2 * B_prec/c(B_tau)))/2)
      B_prec[] = matrix(rgamma(b*p,shape = (B_df + 1)/2,rate = (B_df + B2*c(B_tau))/2),nr = b,nc = p)
      B_prec[] = B_prec*c(B_tau)
      if(resid_intercept){
        B_tau[1,1] = 1e-10
        B_prec[1,] = 1e-10
      }
    }
    if(b_F > 0) {
      B_F2 = B_F^2 * tot_F_prec[rep(1,b_F),]  # need to account for tot_F_prec
      B_F_tau[1,] = rgamma(b_F, shape = fixed_factors_prec_shape + ncol(B_F2)/2, rate = fixed_factors_prec_rate + rowSums(B_F2 * B_F_prec/c(B_F_tau))/2)
      B_F_prec[] = matrix(rgamma(b_F*k,shape = (B_F_df + 1)/2,rate = (B_F_df + B_F2*c(B_F_tau))/2),nr = b_F,nc = k)
      B_F_prec[] = B_F_prec*c(B_F_tau)
      B_F_tau[1,X_F_zero_variance] = 1e10
      B_F_prec[X_F_zero_variance,] = 1e10
    }
  })))
  return(current_state)
}


sample_B_prec_TPB = function(BSFG_state,ncores = detectCores(),cluster=NULL,...){
  # following https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4214276/
  # uses GIGrvg for GIG
  priors         = BSFG_state$priors
  run_variables  = BSFG_state$run_variables
  current_state  = BSFG_state$current_state

  current_state = with(c(priors,run_variables),
                       with(list(
                         # load priors
                         B_A = B_prior$B_A,
                         B_B = B_prior$B_B,
                         B_omega = B_prior$B_omega,

                         B_F_A = B_prior$B_F_A,
                         B_F_B = B_prior$B_F_B,
                         B_F_omega = B_prior$B_F_omega
                       ), within(current_state,{
    # initialize variables if needed
    if(!exists('B_psi')){
      if(b > 0) {
        B_psi = matrix(rgamma(p,shape = 1/2,rate = B_omega),nrow=1)
        B_lambda = matrix(rgamma(b*p,shape = B_B,rate=B_psi[rep(1,b),]),b,p)
        B_prec = 1/matrix(rgamma(b*p,shape = B_A,rate = B_lambda),b,p)
      }
    }
    if(!exists('B_F_psi')){
      if(b_F > 0) {
        B_F_psi = matrix(rgamma(k,shape = 1/2,rate = B_F_omega),nrow=1)
        B_F_lambda = matrix(rgamma(b_F*k,shape = B_F_B,rate=B_F_psi[rep(1,b_F),]),b_F,k)
        B_F_prec = 1/matrix(rgamma(b_F*k,shape = B_F_A,rate = B_F_lambda),b_F,k)
      }
    }
    if(b > 0) { # B[i,j] ~ N(0,B_prec[i,j]^{-1})
      # recover()
      B_psi[1,] = rgamma(p,shape=b*B_B + 1/2, rate = colSums(B_lambda) + B_omega)
      B_lambda[] = rgamma(b*p,shape = B_A + B_B, rate = sweep(1/B_prec,2,B_psi,'+'))

      if(!exists('B_ncores')) B_ncores = ncores
      if(is.null(cluster)) {
        B_prec[] = 1/do.call(cbind,mclapply(1:p,function(j) {
          sapply(1:b,function(i) GIGrvg::rgig(n=1,lambda = B_A-1/2, chi = B[i,j]^2, psi = 2*B_lambda[i,j]))
        },mc.cores = B_ncores))  # this replaces B_tau - BSFG requires precision, not variance.
      } else{
        rm_cluster=FALSE
        if(is.na(cluster[1])){
          rm_cluster=TRUE
          cluster = makeCluster(1)
          print(cluster)
        }
        vars = c('b','p','B_A','B','B_lambda','B_ncores')
        clusterEnv = lapply(vars,function(x) eval(parse(text=x)))
        names(clusterEnv) = vars
        clusterEnv = as.environment(clusterEnv)
        clusterExport(cluster,vars,envir = clusterEnv)
        clusterEvalQ(cluster,library(parallel))
        B_prec[] = 1/do.call(cbind,parLapply(cluster,1,function(x) {
          mclapply(1:p,function(j) {
            sapply(1:b,function(i) GIGrvg::rgig(n=1,lambda = B_A-1/2, chi = B[i,j]^2, psi = 2*B_lambda[i,j]))
          },mc.cores = B_ncores)
        })[[1]])
        clusterEvalQ(cluster,gc)
        if(rm_cluster){
          stopCluster(cluster)
          cluster=NA
        }
        rm(clusterEnv)
        rm(vars)
        gc()
      }

      B_prec[B_prec==0] = 1e-10
      if(resid_intercept){
        B_prec[1,] = 1e-10
      }
    }
    if(b_F > 0) { # B_F[i,j] ~ N(0,B_F_prec[i,j]^{-1} * tot_F_prec[j]^{-1})
      B_F_psi[1,] = rgamma(k,shape=b_F*B_F_B + 1/2, rate = colSums(B_F_lambda) + B_F_omega)
      B_F_lambda[] = rgamma(b_F*k,shape = B_F_A + B_F_B, rate = sweep(1/B_F_prec,2,B_F_psi,'+'))

      if(!exists('B_F_ncores')) B_F_ncores = ncores
      if(is.null(cluster)) {
        B_F_prec[] = 1/do.call(cbind,mclapply(1:k,function(j) {
          sapply(1:b_F,function(i) GIGrvg::rgig(n=1,lambda = B_F_A-1/2, chi = B_F[i,j]^2*tot_F_prec[j], psi = 2*B_F_lambda[i,j]))
        },mc.cores = B_F_ncores))  # this replaces B_F_tau - BSFG requires precision, not variance.
      } else {
        rm_cluster=FALSE
        if(is.na(cluster[1])){
          rm_cluster=TRUE
          cluster = makeCluster(1)
          print(cluster)
        }
        vars = c('k','b_F','B_F_A','B_F','tot_F_prec','B_F_lambda','B_F_ncores')
        clusterEnv = lapply(vars,function(x) eval(parse(text=x)))
        names(clusterEnv) = vars
        clusterEnv = as.environment(clusterEnv)
        clusterExport(cluster,vars,envir = clusterEnv)
        clusterEvalQ(cluster,library(parallel))
        B_F_prec[] = 1/do.call(cbind,parLapply(cluster,1,function(x) {
          mclapply(1:k,function(j) {
            sapply(1:b_F,function(i) GIGrvg::rgig(n=1,lambda = B_F_A-1/2, chi = B_F[i,j]^2*tot_F_prec[j], psi = 2*B_F_lambda[i,j]))
          },mc.cores = B_F_ncores)
        })[[1]])
        clusterEvalQ(cluster,gc)
        if(rm_cluster){
          stopCluster(cluster)
          cluster=NA
        }
        rm(clusterEnv)
        rm(vars)
        gc()
      }

      B_F_prec[B_F_prec==0] = 1e-10
      B_F_prec[X_F_zero_variance,] = 1e10
    }
  })))
  return(current_state)
}


sample_B_prec_QTLBEN = function(BSFG_state,...){
  priors         = BSFG_state$priors
  run_variables  = BSFG_state$run_variables
  current_state  = BSFG_state$current_state

  QTL_columns_resid = BSFG_state$data_matrices$QTL_columns_resid
  QTL_columns_factors = BSFG_state$data_matrices$QTL_columns_factors

  current_state = with(c(priors,run_variables),within(current_state,{
    # load priors
    B_df   = B_prior$B_df
    B_F_df = B_prior$B_F_df

    # initialize variables if needed
    if(!exists('B_tau')){
      if(b > 0) {
        B_tau = matrix(c(1e-10,rgamma(b-1,shape = priors$fixed_resid_prec_shape, rate = priors$fixed_resid_prec_rate)),nrow=1)
      } else{
        B_tau = matrix(0,ncol=0,nrow=1)
      }
      if(b_F > 0) {
        B_F_tau = matrix(rgamma(b_F,shape = priors$fixed_factors_prec_shape, rate = priors$fixed_factors_prec_rate),nrow=1)
      } else{
        B_F_tau = matrix(0,ncol=0,nrow=1)
      }

      B_prec = matrix(B_tau,nrow = b, ncol = p)
      B_F_prec = matrix(B_F_tau,nrow = b_F, ncol = k)
    }

    if(b > 0) {
      non_QTL_resid = 1:nrow(B)
      if(!is.null(QTL_columns_resid)){
        non_QTL_resid = non_QTL_resid[-QTL_columns_resid]
      }
      B2 = B^2
      B_tau[1,] = rgamma(b, shape = fixed_resid_prec_shape + ncol(B2)/2, rate = fixed_resid_prec_rate + rowSums((B2 * B_prec/c(B_tau)))/2)
      # B_tau[1,non_QTL_resid] = rgamma(length(non_QTL_resid), shape = fixed_resid_prec_shape + ncol(B2)/2,
      # rate = fixed_resid_prec_rate + rowSums((B2[non_QTL_resid,,drop=FALSE] * B_prec[non_QTL_resid,,drop=FALSE]/c(B_tau[non_QTL_resid])))/2)
      B_prec[non_QTL_resid,] = matrix(rgamma(length(non_QTL_resid)*p,shape = (B_df + 1)/2,rate = (B_df + B2[non_QTL_resid,,drop=FALSE]*c(B_tau[non_QTL_resid]))/2),nr = length(non_QTL_resid),nc = p)
      B_prec[non_QTL_resid,] = B_prec[non_QTL_resid,,drop=FALSE]*c(B_tau[non_QTL_resid])

      # inverse Gaussian Bayesian Elastic Net from Li and Lin 2010
      B_F_std = B_F[QTL_columns_resid,,drop=FALSE]
      QTL_lambda1 = B_df
      QTL_lambda2 = B_tau[QTL_columns_resid]
      tau = matrix(1/rinvgauss(n=length(QTL_columns_resid)*p,
                               mean = sqrt(QTL_lambda1)/(2*QTL_lambda2*abs(B_F_std)),
                               shape = QTL_lambda1/(4*QTL_lambda2))+1,
                   nr = length(QTL_columns_resid),nc = p)
      B_prec[QTL_columns_resid,] = tau/(tau-1) * QTL_lambda2

      if(resid_intercept){
        B_tau[1,1] = 1e-10
        B_prec[1,] = 1e-10
      }
    }
    if(b_F > 0) {
      non_QTL_factor = 1:nrow(B_F)
      if(!is.null(QTL_columns_factors)){
        non_QTL_factor = non_QTL_factor[-QTL_columns_factors]
      }
      B_F2 = B_F^2
      B_F_tau[1,] = rgamma(b_F, shape = fixed_factors_prec_shape + ncol(B_F2)/2, rate = fixed_factors_prec_rate + rowSums(B_F2 * B_F_prec/c(B_F_tau))/2)
      # B_F_tau[1,non_QTL_factor] = rgamma(length(non_QTL_factor), shape = fixed_factors_prec_shape + ncol(B_F2)/2,
      #                                 rate = fixed_factors_prec_rate + rowSums((B_F2[non_QTL_factors,,drop=FALSE] * B_F_prec[non_QTL_factors,,drop=FALSE]/c(B_F_tau[non_QTL_factors])))/2)
      B_F_prec[non_QTL_factor,] = matrix(rgamma(length(non_QTL_factor)*k,shape = (B_df + 1)/2,rate = (B_df + B_F2[non_QTL_factor,,drop=FALSE]*c(B_F_tau[non_QTL_factor]))/2),nr = length(non_QTL_factor),nc = k)
      B_F_prec[non_QTL_factor,] = B_F_prec[non_QTL_factor,,drop=FALSE]*c(B_F_tau[non_QTL_factor])

      # inverse Gaussian Bayesian Elastic Net from Li and Lin 2010
      B_F_tau[QTL_columns_factors] = 1
      QTL_lambda1 = B_df
      QTL_lambda2 = B_df
      tau = matrix(1/rinvgauss(n=length(QTL_columns_factors)*k,
                               mean = sqrt(QTL_lambda1)/(2*QTL_lambda2*abs(B_F[QTL_columns_factors,,drop=FALSE])),
                               shape = rep(QTL_lambda1*B_F_tau[QTL_columns_factors]/(4*QTL_lambda2),k)
      )+1,
      nr = length(QTL_columns_factors),nc = k)
      B_F_prec[QTL_columns_factors,] = tau/(tau-1) * B_F_tau[QTL_columns_factors] * QTL_lambda2
      B_F_prec[B_F_prec<0] = 1e-10
      if(min(B_F_prec) < 0) recover()
      rm(tau)

      B_F_tau[1,X_F_zero_variance] = 1e10
      B_F_prec[X_F_zero_variance,] = 1e10
    }
  }))
  return(current_state)
}
