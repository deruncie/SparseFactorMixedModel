fast_BSFG_sampler_comb = function(BSFG_state,n_samples) {
  # -- Daniel Runcie -- #
  #
  # Gibbs sampler for genetic covariance estimation based on mixed effects
  # model, with missing data
  #
  # Code for:
  # 
  # Runcie and Mukherjee (2013) Dissecting high-dimensional traits
  # with Bayesian sparse factor analysis of genetic covariance matrices.
  # GENETICS.
  #
  # (c) July 30, 2013
  #
  # code based on original provided by Anirban Bhattacharya
  #
  #         
  # This function implements the BSF-G partially collapsed Gibbs sampler.
  # Variable notation follows the Runcie and Mukherjee manuscript as closely
  # as possible.
  #
  # All input and initialization functions are carried out by the function
  # fast_BSFG_sampler_init See the documentation of that function for details.
  # 
  # The sampler is designed to do a short-medium length run and then return
  # the state of the chain. After assessing the progress of the chain,
  # (optionally), the Posterior matrices can be reset, and the chain
  # re-started where it left off. Right now, the random seed is not saved. Probably should
  # be if I can figure out how.
  # 
  # This function takes the following inputs:
  #     data_matrices: struct holding the (should be imutable) data, design and incidence matrices:
  #          Y:  		Full phenotype data. n x p
  # 		   X: 		Fixed effect design matrix. n x b
  #          Z_1:     Random effect 1 incidence matrix. n x r1
  #          Z_2:     Random effect 2 incidence matrix. n x r2
  #          Y_missing: incidence matrix of missing data in Y
  #     start_i: iteration number of the end of the last run.
  #     draw_iter: frequency of updating diagnostic plots
  #     burn: number of burnin samples
  #     sp: total number of samples to collect
  #     thin: thinning rate of chain
  #     simulation: boolean. Is this a simulation?
  #     params: struct with chain parameters.
  #     priors: struct with all relevant prior hyperparameters
  #     Posterior: struct with posterior matrices, or running posterior means. 
  #            Note: not sure how well posterior means work after re-starting chain. Probably not well.
  #     current_state: current (initial) conditions of all model parameters
  # 
  # Several diagnostic plots are produced during the run. 
  #     Their interpretation is described within the source codes:
  #         draw_simulation_diagnostics.m: For simulated data with known true values
  #         draw_results_diagnostics.m: Otherwise
  #         
  
  data_matrices  = BSFG_state$data_matrices
  params         = BSFG_state$params
  priors         = BSFG_state$priors
  Posterior      = BSFG_state$Posterior
  current_state  = BSFG_state$current_state
  run_parameters = BSFG_state$run_parameters
  run_variables  = BSFG_state$run_variables
  traitnames     = BSFG_state$traitnames
  
  # ----------------------------------------------- #
  # ----------------Load data matrices------------- #
  # ----------------------------------------------- #
  
  Y  			  = data_matrices$Y
  Z_1     	= data_matrices$Z_1
  Z_2     	= data_matrices$Z_2
  X       	= data_matrices$X
  Y_missing = data_matrices$Y_missing
  
  p   = run_variables$p
  n   = run_variables$n
  r   = run_variables$r
  r2  = run_variables$r2
  b   = run_variables$b
  
  
  # ----------------------------------------------- #
  # ----------------Load priors-------------------- #
  # ----------------------------------------------- #
  
  resid_Y_prec_shape =   priors$resid_Y_prec_shape
  resid_Y_prec_rate  =   priors$resid_Y_prec_rate
  E_a_prec_shape     =   priors$E_a_prec_shape
  E_a_prec_rate      =   priors$E_a_prec_rate
  W_prec_shape       =   priors$W_prec_shape
  W_prec_rate        =   priors$W_prec_rate
  Lambda_df          =   priors$Lambda_df
  delta_1_shape      =   priors$delta_1_shape
  delta_1_rate       =   priors$delta_1_rate
  delta_2_shape      =   priors$delta_2_shape
  delta_2_rate       =   priors$delta_2_rate
  h2_priors_factors  =   priors$h2_priors_factors
  
  
  # ----------------------------------------------- #
  # ----------------Load current state------------- #
  # ----------------------------------------------- #
  resid_Y_prec =   current_state$resid_Y_prec
  F            =   current_state$F
  Lambda       =   current_state$Lambda
  E_a          =   current_state$E_a
  W            =   current_state$W
  Lambda_prec  =   current_state$Lambda_prec
  delta        =   current_state$delta
  tauh         =   current_state$tauh
  E_a_prec     =   current_state$E_a_prec
  W_prec       =   current_state$W_prec
  Plam         =   current_state$Plam
  F_h2         =   current_state$F_h2
  F_a          =   current_state$F_a
  B            =   current_state$B
  start_i      =   current_state$nrun
  k = ncol(F)
  
  # ----------------------------------------------- #
  # -----------Reset Global Random Number Stream--- #
  # ----------------------------------------------- #
  do.call("RNGkind",as.list(BSFG_state$RNG$RNGkind))  ## must be first!
  assign(".Random.seed", BSFG_state$RNG$Random.seed, .GlobalEnv)
  
  # ----------------------------------------------- #
  # -----------Load pre-calcualted matrices-------- #
  # ----------------------------------------------- #
  Ainv                             = run_variables$Ainv
  A_2_inv                          = run_variables$A_2_inv
  invert_aI_bZAZ                   = run_variables$invert_aI_bZAZ  
  invert_aPXA_bDesignDesignT       = run_variables$invert_aPXA_bDesignDesignT 
  invert_aZZt_Ainv                 = run_variables$invert_aZZt_Ainv
  invert_aPXA_bDesignDesignT_rand2 = run_variables$invert_aPXA_bDesignDesignT_rand2  
  
  
  # ----------------------------------------------- #
  # ----------------Set up run--------------------- #
  # ----------------------------------------------- #
  #     b0,b1: parameters controlling rate of adaptation of factor model size
  #     h2_divisions: number of discrete steps for each factor heritability parameter
  #     epsilon: truncation point for factor loadings during adaptation
  b0           = run_parameters$b0
  b1           = run_parameters$b1
  h2_divisions = run_parameters$h2_divisions
  epsilon      = run_parameters$epsilon
  prop         = run_parameters$prop
  save_freq    = run_parameters$save_freq
  burn         = run_parameters$burn
  thin         = run_parameters$thin
  draw_iter    = run_parameters$draw_iter
  simulation   = run_parameters$simulation
  
  
  # ----------------------------------------------- #
  # ---Extend posterior matrices for new samples--- #
  # ----------------------------------------------- #
  # change posterior into a list (can the elements in dataframe be lists?)
  sp = (start_i + n_samples - burn)/thin - ncol(Posterior$Lambda)
  if(sp > 0){
    Posterior$Lambda        = cbind(Posterior$Lambda,matrix(0,nr = nrow(Posterior$Lambda),nc = sp))
    Posterior$F             = cbind(Posterior$F,matrix(0,nr = nrow(Posterior$F),nc = sp))
    Posterior$F_a           = cbind(Posterior$F_a,matrix(0,nr = nrow(Posterior$F_a),nc = sp))
    Posterior$delta         = cbind(Posterior$delta,matrix(0,nr = nrow(Posterior$delta),nc = sp))
    Posterior$F_h2          = cbind(Posterior$F_h2,matrix(0,nr = nrow(Posterior$F_h2),nc = sp))
    Posterior$resid_Y_prec  = cbind(Posterior$resid_Y_prec,matrix(0,nr = nrow(Posterior$resid_Y_prec),nc = sp))
    Posterior$E_a_prec      = cbind(Posterior$E_a_prec,matrix(0,nr = nrow(Posterior$E_a_prec),nc = sp))
    Posterior$W_prec        = cbind(Posterior$W_prec,matrix(0,nr = nrow(Posterior$W_prec),nc = sp))
  }
  
  # ----------------------------------------------- #
  # --------------start gibbs sampling------------- #
  # ----------------------------------------------- #
  
  start_time = Sys.time()
  for(i in start_i+(1:n_samples)){
    
    # -----fill in missing phenotypes----- #
    #conditioning on everything else
    if(sum(Y_missing)>0) {
      meanTraits = X %*% B + F %*% t(Lambda) + Z_1 %*% E_a + Z_2 %*% W
      resids = matrix(rnorm(p*n,0,sqrt(1/resid_Y_prec)),nr = n,nc = p,byrow=T)
      Y[Y_missing] = meanTraits[Y_missing] + resids[Y_missing]
    }
    # recover()
    
    # -----Sample Lambda------------------ #
    #conditioning on W, B, F, marginalizing over E _a
    Y_tilde = Y - X %*% B - Z_2 %*% W
    Y_tilde = as.matrix(Y_tilde)
    # Lambda = sample_Lambda( Y_tilde,F,resid_Y_prec, E_a_prec,Plam,invert_aI_bZAZ )
    
    Lambda = sample_Lambda_comb( Y_tilde,F,resid_Y_prec, E_a_prec,Plam,invert_aI_bZAZ )

    # -----Sample Lambda_prec------------- #
    Lambda2 = Lambda^2
    Lambda_prec = matrix(rgamma(p*k,shape = (Lambda_df + 1)/2,rate = (Lambda_df + sweep(Lambda2,2,tauh,'*'))/2),nr = p,nc = k)
    
    # -----Sample delta, update tauh------ #
    delta = sample_delta( delta,tauh,Lambda_prec,delta_1_shape,delta_1_rate,delta_2_shape,delta_2_rate,Lambda2 )
    tauh  = cumprod(delta)
    
    # -----Update Plam-------------------- #
    Plam = sweep(Lambda_prec,2,tauh,'*')
    
    # -----Sample resid_Y_prec ---------------- #
    # should be able to sample this here marginalizing over E_a in the same way that E_a_prec is done below.
    
    # -----Sample E_a_prec---------------- #
    #conditioning on W, B, F, marginalizing over E_a 
    # E_a_prec = sample_prec_discrete_conditional(Y - X %*% B-F %*% t(Lambda) - Z_2 %*% W,h2_divisions,h2_priors_resids,invert_aI_bZAZ,resid_Y_prec)
    # E_a_prec = sample_prec_discrete_conditional_c(Y - X %*% B-F %*% t(Lambda) - Z_2 %*% W,h2_divisions,h2_priors_resids,invert_aI_bZAZ,resid_Y_prec)
    
    # pops should be a output from setup datasets and for pop in pops, the row number for 
    # each pop should also included in the setup dataset
    # read pops (number of populations)
    # read row number for each pop
    # in the loop, replace Y as Y[pops,], etc (instead of using lists because we need to use all Y together to sample lambda)
   
    # -----Sample B and E_a--------------- #
    #conditioning on W, F, Lambda
    for(pop in pops){
      Y_tilde = Y[[pop]] - F[[pop]] %*% t(Lambda) - Z_2[[pop]] %*% W[[pop]]
      Y_tilde = as.matrix(Y_tilde)
      # location_sample = sample_means( Y_tilde, resid_Y_prec, E_a_prec, invert_aPXA_bDesignDesignT )
      location_sample = sample_means( Y_tilde, resid_Y_prec[[pop]], E_a_prec[[pop]], invert_aPXA_bDesignDesignT[[pop]] )
      B[[pop]]   = location_sample[1:b,]
      if(length(B[[pop]])==p){
        B[[pop]] = matrix(0,nr=0,nc=p)
      }
      E_a[[pop]] = location_sample[b+(1:r),]
      
      
      # -----Sample W ---------------------- #
      #conditioning on B, E_a, F, Lambda
      if(ncol(Z_2) > 0) {
        Y_tilde = Y[[pop]] - X[[pop]] %*% B[[pop]] - Z_1[[pop]] %*% E_a[[pop]] - F[[pop]] %*% t(Lambda)
        Y_tilde = as.matrix(Y_tilde)
        # location_sample = sample_means( Y_tilde, resid_Y_prec, W_prec, invert_aPXA_bDesignDesignT_rand2 )
        location_sample = sample_means( Y_tilde, resid_Y_prec[[pop]], W_prec[[pop]], invert_aPXA_bDesignDesignT_rand2[[pop]] )
        W[[pop]] = location_sample
      }
      
      # -----Sample F_h2-------------------- #
      #conditioning on F, marginalizing over F_a
      # F_h2 = sample_h2s_discrete(F,h2_divisions,h2_priors_factors,invert_aI_bZAZ)
      F_h2[[pop]] = sample_h2s_discrete(F[[pop]],h2_divisions,h2_priors_factors,invert_aI_bZAZ[[pop]])
      
      # -----Sample F_a--------------------- #
      #conditioning on F, F_h2
      # F_a = sample_F_a_c(F,Z_1,F_h2,invert_aZZt_Ainv)
      F_a[[pop]] = sample_F_a(F[[pop]],Z_1[[pop]],F_h2[[pop]],invert_aZZt_Ainv[[pop]]);
      
      # -----Sample F----------------------- #
      #conditioning on B,F_a,E_a,W,Lambda, F_h2
      Y_tilde = Y[[pop]] - X[[pop]] %*% B[[pop]] - Z_1[[pop]] %*% E_a[[pop]] - Z_2[[pop]] %*% W[[pop]]
      Y_tilde = as.matrix(Y_tilde)
      # F1 = sample_factors_scores( Y_tilde, Z_1,Lambda,resid_Y_prec,F_a,F_h2 )
      F[[pop]] = sample_factors_scores(Y_tilde, Z_1[[pop]],Lambda,resid_Y_prec[[pop]],F_a[[pop]],F_h2[[pop]] )
      
      
      # -----Sample resid_Y_prec------------ #
      Y_tilde = Y[[pop]] - X[[pop]] %*% B[[pop]] - F[[pop]] %*% t(Lambda) - Z_1[[pop]] %*% E_a[[pop]] - Z_2[[pop]] %*% W[[pop]]
      resid_Y_prec[[pop]] = rgamma(p,shape = resid_Y_prec_shape + n/2, rate = resid_Y_prec_rate+1/2*colSums(Y_tilde^2)) #model residual precision
      # n is the sum of all populations
      
      # -----Sample E_a_prec---------------- #
      #     #random effect 1 (D) residual precision
      E_a_prec[[pop]] = rgamma(p,shape = E_a_prec_shape + r/2, rate = E_a_prec_rate + 1/2*diag(t(E_a[[pop]]) %*% Ainv[[pop]] %*% E_a[[pop]]))
      
      # -----Sample W_prec------------------ #
      if(ncol(Z_2) > 0) {
        W_prec[[pop]] =  rgamma(p, W_prec_shape + r2/2,rate = W_prec_rate + 1/2*diag(t(W[[pop]]) %*% A_2_inv[[pop]] %*% W[[pop]]))
      }
    }
      # -- adapt number of factors to samples ---#
      # adapt_result = update_k( F,Lambda,F_a,F_h2,Lambda_prec,Plam,delta,tauh,Z_1,Lambda_df,delta_2_shape,delta_2_rate,b0,b1,i,epsilon,prop )
      # F           = adapt_result$F
      # Lambda      = adapt_result$Lambda
      # F_a         = adapt_result$F_a
      # F_h2        = adapt_result$F_h2
      # Lambda_prec = adapt_result$Lambda_prec
      # Plam        = adapt_result$Plam
      # delta       = adapt_result$delta
      # tauh        = adapt_result$tauh
      # k = ncol(F)
      
      # -- save sampled values (after thinning) -- #
      if( (i-burn) %% thin == 0 && i > burn) {
        
        sp_num = (i-burn)/thin
        
        Posterior = save_posterior_samples_comb( sp_num,Posterior,Lambda,F,F_a,B,W,E_a,delta,F_h2,resid_Y_prec,E_a_prec,W_prec)
        
        # reorder elements in Lambda
        Lambda_pos  = Posterior$Lambda
        Lam_rowmean = rowMeans(Lambda_pos)
        Lambda_hat  = matrix(Lam_rowmean,nr=p)
        
        
        if(sp_num %% save_freq == 0) save(Posterior,file='Posterior.RData')
      }
        
      
      # -- provide run diagnostics and plots -- #
      # if(i %% draw_iter  == 0) {
      # directory=strread(pwd,'#s','delimiter','/')
      # disp(directory(end))
      # disp(i)
      # elapsed=toc
      # #output some running statistics on the current factors and their genetic variances
      # # each row represents a factor:
      # #    delta[i], [i], F_h2[i], V_P(F_i), V_A(F_i)
      # [delta,[1:size(F,2)]',F_h2,sum(F.^2)'./(size(F,1)-1),sum(F_a.^2)'./(size(F_a,1)-1)]
      # disp(strcat('Time remaining:',num2str((nrun-i) * (elapsed/i) * 1/60)))
      
      #make some plots of some running statistics
      
      #   draw_results_diagnostics(sp_num,run_parameters,run_variables,Lambda, F_h2, Posterior,E_a_prec,resid_Y_prec,traitnames)
      # }
    }

     end_time = Sys.time()
     print(end_time - start_time)
     # save(Posterior,file = 'Posterior.RData')


  # ----------------------------------------------- #
  # ------------Save state for restart------------- #
  # ----------------------------------------------- #
  
  current_state[["resid_Y_prec"]]  = resid_Y_prec
  current_state[["Lambda_prec"]]   = Lambda_prec
  current_state[["delta"]]         = delta
  current_state[["tauh"]]          = tauh
  current_state[["Plam"]]          = Plam
  current_state[["Lambda"]]        = Lambda
  current_state[["F_h2"]]          = F_h2
  current_state[["E_a_prec"]]      = E_a_prec
  current_state[["W_prec"]]        = W_prec
  current_state[["F_a"]]           = F_a
  current_state[["F"]]             = F
  current_state[["E_a"]]           = E_a
  current_state[["B"]]             = B
  current_state[["W"]]             = W
  current_state[["nrun"]] 			= i
  
  # save(current_state,file='current_state.RData')
  

  BSFG_state$current_state = current_state
  BSFG_state$Posterior = Posterior
  BSFG_state$RNG = list(
    Random.seed = .Random.seed,
    RNGkind = RNGkind()
  )
  
  save(BSFG_state,file = 'BSFG_state.RData')
  return(BSFG_state)
} 