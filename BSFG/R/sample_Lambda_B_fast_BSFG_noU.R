sample_Lambda_B.fast_BSFG_noU = function(BSFG_state,grainSize,...) {
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

    # -----Sample Lambda and B ------------------ #
    #conditioning on F, marginalizing over U_R
    Design = cBind(X,F)
    rows = b + k
    prior_mean = matrix(0,rows,p)
    if(b > 0) {
      prior_prec = rbind(B_prec,t(Plam))
    } else{ # b == 0
      prior_prec = t(Plam)
    }
    if(is.null(cis_genotypes)){
      # recover()
      randn_theta = matrix(rnorm(rows*p),rows)
      randn_e = matrix(rnorm(n*p),n)
      # UtE = toDense(Ut %*% Eta)
      # UtD = toDense(Ut %*% Design)
      # microbenchmark(
        # {
          coefs = sample_coefs_parallel_sparse_c_Eigen( Ut,Eta,Design,resid_h2, tot_Eta_prec,s, prior_mean,prior_prec,randn_theta,randn_e,grainSize)
      #     coefs2 = sample_coefs_parallel_sparse_c_Eigen2( Ut,Eta,Design,resid_h2, tot_Eta_prec,s, prior_mean,prior_prec,randn_theta,randn_e,grainSize)
      # )
        # },
      # {
      # },{
        # h2_mat = resid_h2[rep(1,n),]
        # R_var = (h2_mat*s + 1-h2_mat) / tot_Eta_prec[rep(1,n),]
      # coefs2 = sample_coefs_parallel_sparse_c_Eigen2( Ut,Eta,Design,R_var, prior_mean,prior_prec,randn_theta,randn_e,grainSize)
      #   },{
      #     UtE = toDense(Ut %*% Eta)
      #     UtD = toDense(Ut %*% Design)
      #     coefs2 = sample_coefs_parallel_sparse_c_Eigen2(UtE,UtD,resid_h2, tot_Eta_prec,s, prior_mean,prior_prec,randn_theta,randn_e,grainSize)
      #   },{
      # R_prec = uncorrelated_prec_mat(resid_h2,tot_Eta_prec,s)
      # # coefs3 = sample_coefs_parallel_sparse_c_Eigen3( (Ut %*% Eta),(Ut %*% Design),R_prec, prior_mean,prior_prec,randn_theta,randn_e,grainSize)
      # coefs3 = sample_coefMat_uncorrelated_parallel_Eigen(UtE,UtD,R_prec, prior_mean,prior_prec,randn_theta,randn_e,grainSize)
      # })
    # )
      # coefs = rbind(B,t(Lambda))
      # groups = seq(1,rows,by=5)
      # randn_e = matrix(rnorm(n*p*length(groups)),n*length(groups))
      # coefs=sample_coefs_parallel_sparse_c_Eigen_group( Ut,Eta,Design,coefs,groups,resid_h2, tot_Eta_prec,s, prior_mean,prior_prec,randn_theta,randn_e,grainSize)
      if(b > 0){
        B[] = coefs[1:b,,drop=FALSE]
      }
      Lambda[] = t(coefs[b + 1:k,,drop=FALSE])
    } else{
      randn_theta = matrix(rnorm(rows*p),rows)
      randn_e = matrix(rnorm(n*p),n)
      randn_cis = rnorm(cis_effects_index[length(cis_effects_index)]-1)
      result = sample_cis_coefs_parallel_sparse_c_Eigen(Ut,Eta,toDense(Design),cis_genotypes,resid_h2,tot_Eta_prec,s,prior_mean,prior_prec,
                                                        randn_theta,randn_e,randn_cis,cis_effects_index-1,1)
      if(b > 0){
        B[] = result[[1]][1:b,,drop=FALSE]
      }
      Lambda[] = t(result[[1]][b+1:k,,drop=FALSE])
      cis_effects[] = result[[2]]
    }
  }))
  current_state = current_state[current_state_names]

  return(current_state)
}
