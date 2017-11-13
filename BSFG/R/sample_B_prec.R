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
                             B_tau = matrix(c(1e-10,rgamma(b-1,shape = fixed_resid_prec_shape, rate = fixed_resid_prec_rate)),nrow=1)
                           } else{
                             B_tau = matrix(0,nrow=1,ncol=0)
                           }
                           if(b_F > 0) {
                             B_F_tau = matrix(rgamma(b_F,shape = fixed_factors_prec_shape, rate = fixed_factors_prec_rate),nrow=1)
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

                           # B_F_tau[1,] = rgamma(b_F, shape = fixed_factors_prec_shape + ncol(B_F2)/2, rate = fixed_factors_prec_rate + rowSums(B_F2 * B_F_prec/c(B_F_tau))/2)
                           # B_F_prec[] = matrix(rgamma(b_F*k,shape = (B_F_df + 1)/2,rate = (B_F_df + B_F2*c(B_F_tau))/2),nr = b_F,nc = k)

                           B_F_tau[1,] = rgamma(b_F, shape = fixed_factors_prec_shape + ncol(B_F2)/2, rate = fixed_factors_prec_rate + rowSums(B_F2 * B_F_prec/t(B_F_tau)[,rep(1,k)])/2)
                           B_F_prec[] = matrix(rgamma(b_F*k,shape = (B_F_df + 1)/2,rate = (B_F_df + B_F2*t(B_F_tau)[,rep(1,k)])/2),nr = b_F,nc = k)

                           B_F_prec[] = B_F_prec*t(B_F_tau)[,rep(1,k)]
                           B_F_tau[1,X_F_zero_variance] = 1e10
                           B_F_prec[X_F_zero_variance,] = 1e10
                         }
                       })))
  return(current_state)
}

sample_B_prec_combined_ARD = function(BSFG_state,...){
  priors         = BSFG_state$priors
  run_variables  = BSFG_state$run_variables
  current_state  = BSFG_state$current_state

  current_state = with(c(priors,run_variables),
                       with(list(
                         # load priors
                         B_df   = B_prior$B_df,
                         B_F_df = B_prior$B_F_df,
                         B_QTL_df   = B_prior$B_QTL_df,
                         B_QTL_F_df = B_prior$B_QTL_F_df
                       ),within(current_state,{

                         # initialize variables if needed
                         if(!exists('B_tau')){
                           if(b > 0) {
                             B_tau = matrix(c(1e-10,rgamma(b-1,shape = fixed_resid_prec_shape, rate = fixed_resid_prec_rate)),nrow=1)
                           } else{
                             B_tau = matrix(0,nrow=1,ncol=0)
                           }
                           if(b_F > 0) {
                             B_F_tau = matrix(rgamma(b_F,shape = fixed_factors_prec_shape, rate = fixed_factors_prec_rate),nrow=1)
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
                                       shape = fixed_resid_prec_shape[ib] + ncol(B2)/2 + ncol(B_F2)/2,
                                       rate = fixed_resid_prec_rate[ib] +
                                              rowSums((B2[ib,,drop=F] * B_prec[ib,,drop=F]/c(B_tau)[ib]))/2 +
                                              rowSums(B_F2[ib_F,,drop=F] * B_F_prec[ib_F,,drop=F]/c(B_F_tau)[ib_F])/2)
                           B_tau[1,ib] = tau
                           B_F_tau[1,ib_F] = tau
                         }
                         if(length(fixed_effects_only_resid) > 0){
                           ib = fixed_effects_only_resid
                           tau = rgamma(length(ib),
                                        shape = fixed_resid_prec_shape[ib] + ncol(B2)/2,
                                        rate = fixed_resid_prec_rate[ib] +
                                          rowSums((B2[ib,,drop=F] * B_prec[ib,,drop=F]/c(B_tau)[ib]))/2)
                           B_tau[1,ib] = tau

                         }
                         if(length(fixed_effects_only_factors) > 0){
                           ib = fixed_effects_only_factors
                           tau = rgamma(length(ib),
                                        shape = fixed_resid_prec_shape[ib] + ncol(B_F2)/2,
                                        rate = fixed_resid_prec_rate[ib] +
                                          rowSums(B_F2[ib,,drop=F] * B_F_prec[ib,,drop=F]/c(B_F_tau)[ib])/2)
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


                         if(b_QTL + b_QTL_F > 0){
                           if(!exists('B_QTL_tau')){
                             B_QTL_tau = matrix(1,1,b_QTL)
                             B_QTL_F_tau = matrix(1,1,b_QTL_F)
                             B_QTL_prec = matrix(1,b_QTL,p)
                             B_QTL_F_prec = matrix(1,b_QTL_F,k)
                           }
                           B_QTL2 = B_QTL^2
                           B_QTL_F2 = B_QTL_F^2 * tot_F_prec[rep(1,b_QTL_F),]

                           tau = rgamma(1,
                                        shape = fixed_factors_QTL_prec_shape[1] + length(B_QTL2)/2 + length(B_QTL2)/2,
                                        rate = fixed_factors_QTL_prec_rate[1] +
                                                  sum(B_QTL2 * B_QTL_prec/c(B_QTL_tau))/2 +
                                                  sum(B_QTL_F2 * B_QTL_F_prec/c(B_QTL_F_tau))/2)
                           B_QTL_tau[] = tau
                           B_QTL_F_tau[] = tau

                           B_QTL_prec[] = matrix(rgamma(b_QTL*p,shape = (B_QTL_df + 1)/2,rate = (B_QTL_df + B_QTL2*c(B_QTL_tau))/2),nr = b_QTL,nc = p)
                           B_QTL_F_prec[] = matrix(rgamma(b_QTL_F*k,shape = (B_QTL_F_df + 1)/2,rate = (B_QTL_F_df + B_QTL_F2*t(B_QTL_F_tau)[,rep(1,k)])/2),nr = b_QTL_F,nc = k)
                           B_QTL_F_prec[] = B_QTL_F_prec * t(B_QTL_F_tau)[,rep(1,k)]
                         }

                       })))
  return(current_state)
}

sample_B_prec_ARD_v2 = function(BSFG_state,...){
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
                             B_tau = matrix(c(1e-10,rgamma(b-1,shape = fixed_resid_prec_shape, rate = fixed_resid_prec_rate)),nrow=1)
                           } else{
                             B_tau = matrix(0,nrow=1,ncol=0)
                           }
                           if(b_F > 0) {
                             B_F_beta = matrix(rgamma(b_F,shape = fixed_factors_prec_shape, rate = 1/fixed_factors_prec_rate),nrow=1)
                           } else{
                             B_F_beta = matrix(0,nrow=1,ncol=0)
                           }

                           B_prec = matrix(B_tau,nrow = b, ncol = p)
                           B_F_prec = matrix(B_F_beta,nrow = b_F, ncol = k)
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

                           B_F_beta[1,] = rgamma(b_F,shape = fixed_factors_prec_shape + k*B_F_df, rate = 1/fixed_factors_prec_rate + rowSums(B_F_prec))
                           B_F_prec[] = rgamma(b_F*k,shape = (B_F_df + 1)/2,rate = (t(B_F_beta)[,rep(1,k)] + B_F2)/2)

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
                         B_F_omega = B_prior$B_F_omega,
                         B_F_ncores = B_prior$B_F_ncores
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
                           B_lambda[] = rgamma(b*p,shape = B_A + B_B, rate = 1/B_prec + B_psi[rep(1,b),])
                           B_prec[] = 1/rgig_multiple(b*p,lambda = rep(B_A-1/2,b*p),chi = c(B^2),psi = 2*c(B_lambda))
                           B_prec[B_prec==0] = 1e-10
                           if(resid_intercept){
                             B_prec[1,] = 1e-10
                           }
                         }
                         if(b_F > 0) { # B_F[i,j] ~ N(0,B_F_prec[i,j]^{-1} * tot_F_prec[j]^{-1})
                           B_F_psi[1,] = 1e-7#rgamma(k,shape=sum(!X_F_zero_variance)*B_F_B + 1/2, rate = colSums(B_F_lambda[!X_F_zero_variance,]) + B_F_omega)
                           B_F_lambda[] = rgamma(b_F*k,shape = B_F_A + B_F_B, rate = 1/B_F_prec + B_F_psi[rep(1,b_F),])
                           B_F_prec[] = 1/rgig_multiple(b_F*k,lambda = rep(B_F_A-1/2,b_F*k),chi = c(B_F^2*tot_F_prec[rep(1,b_F),]),psi = 2*c(B_F_lambda))
                           B_F_prec[B_F_prec < 1e-10] = 1e-10
                           B_F_prec[B_F_prec > 1e10] = 1e10
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
