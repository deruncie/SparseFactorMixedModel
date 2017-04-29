sample_Lambda_B.fast_BSFG = function(BSFG_state,grainSize,...) {
  data_matrices  = BSFG_state$data_matrices
  priors         = BSFG_state$priors
  run_parameters = BSFG_state$run_parameters
  run_variables  = BSFG_state$run_variables
  current_state  = BSFG_state$current_state
  invert_aI_bZKZ = BSFG_state$run_variables$invert_aI_bZKZ
  Ut = t(invert_aI_bZKZ$U)
  s = invert_aI_bZKZ$s

  current_state_names = names(current_state)
  current_state = with(c(priors,run_parameters, run_variables,data_matrices),within(current_state, {
    k = ncol(Lambda)
    # UtEta = as.matrix(Ut %*% Eta)

    # -----Sample Lambda and B ------------------ #
    #conditioning on F, marginalizing over U_R
    Design = as.matrix(cbind(X,F))
    # UtDesign = as.matrix(Ut %*% Design)
    rows = b + k
    prior_mean = matrix(0,rows,p)
    if(b > 0) {
      prior_prec = rbind(prec_B,t(Plam))
    } else{ # b == 0
      prior_prec = t(Plam)
    }
    if(is.null(cis_genotypes)){
      # coefs = sample_coefs_parallel_sparse_c( UtEta,UtDesign,resid_h2, tot_Eta_prec,s, prior_mean,prior_prec,grainSize)
      randn_theta = matrix(rnorm(rows*p),rows)
      randn_e = matrix(rnorm(n*p),n)
      coefs = sample_coefs_parallel_sparse_c_Eigen( Ut,Eta,Design,resid_h2, tot_Eta_prec,s, prior_mean,prior_prec,randn_theta,randn_e,grainSize)
      if(b > 0){
        B[] = coefs[1:b,,drop=FALSE]
      }
      Lambda[] = t(coefs[b + 1:k,,drop=FALSE])
    } else{
      for(j in 1:p){
        cis_X_j = cis_genotypes[[j]]
        if(var(cis_X_j) > 0) {   # Temporary fix
          # UtDesign_j = Ut %*% cis_X_j
          # if(is(UtDesign_j,'Matrix')) UtDesign_j = UtDesign_j@x
          Design_j = cbind(Design,cis_X_j)
          prior_mean_j = rbind(prior_mean[,j,drop=FALSE],matrix(0,ncol(cis_X_j)))
          prior_prec_j = rbind(prior_prec[,j,drop=FALSE],matrix(1e-10,ncol(cis_X_j)))
          randn_theta = matrix(rnorm(ncol(Design_j)),ncol(Design_j))
          randn_e = matrix(rnorm(n),n)
          # coefs_j = sample_coefs_parallel_sparse_c(UtEta[,j,drop=FALSE],cbind(UtDesign,UtDesign_j),resid_h2[,j,drop=FALSE], tot_Eta_prec[,j,drop=FALSE],s,prior_mean_j,prior_prec_j,grainSize)
          coefs_j = sample_coefs_parallel_sparse_c_Eigen(Ut,Eta[,j,drop=FALSE],Design_j,
                                                         resid_h2[,j,drop=FALSE], tot_Eta_prec[,j,drop=FALSE],
                                                         s,prior_mean_j,prior_prec_j,
                                                         randn_theta,randn_e,
                                                         grainSize)
          if(b > 0){
            B[,j] = coefs_j[1:b]
          }
          Lambda[j,] = coefs_j[b+1:k]
          cis_effects[,cis_effects_index[j]] = coefs_j[-c(1:(b+k))]
        } else{
          prior_mean_j = prior_mean[,j,drop=FALSE]
          prior_prec_j = prior_prec[,j,drop=FALSE]
          randn_theta = matrix(rnorm(rows),rows)
          randn_e = matrix(rnorm(n),n)
          # coefs_j = sample_coefs_parallel_sparse_c(UtEta[,j,drop=FALSE],UtDesign,resid_h2[,j,drop=FALSE], tot_Eta_prec[,j,drop=FALSE],s,prior_mean_j,prior_prec_j,grainSize)
          coefs_j = sample_coefs_parallel_sparse_c_Eigen(Ut,Eta[,j,drop=FALSE],Design,
                                                         resid_h2[,j,drop=FALSE], tot_Eta_prec[,j,drop=FALSE],
                                                         s,prior_mean_j,prior_prec_j,
                                                         randn_theta,randn_e,
                                                         grainSize)
          if(b > 0){
            B[,j] = coefs_j[1:b]
          }
          Lambda[j,] = coefs_j[b+1:k]
          cis_effects[,cis_effects_index[j]] = 0
        }
      }
    }
  }))
  current_state = current_state[current_state_names]

  return(current_state)
}
