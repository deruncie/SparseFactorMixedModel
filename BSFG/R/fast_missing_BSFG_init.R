initialize_BSFG.fast_missing_BSFG = function(BSFG_state, K_mats = NULL, chol_Ki_mats = NULL,verbose=T,
                                     invert_aI_bZKZ   = NULL,
                                     invert_aZZt_Kinv = NULL,...){

    Y_missing  = BSFG_state$run_parameters$observation_model_parameters$Y_missing
    Z          = BSFG_state$data_matrices$Z
    h2s_matrix = BSFG_state$data_matrices$h2s_matrix

    RE_names   = rownames(h2s_matrix)
    n_RE       = length(RE_names)
    stopifnot(n_RE == 1)  # must be only 1 random effect

    run_parameters = BSFG_state$run_parameters
    run_variables  = BSFG_state$run_variables

# ------------------------------------ #
# ----Precalculate some matrices------ #
# ------------------------------------ #

    #invert the random effect covariance matrices
    K = K_mats[[1]]
    chol_Kinv = chol_Ki_mats[[1]]

    #pre-calculate transformation parameters to diagonalize aI + bZKZ for fast
    #inversion: inv(aI + bZKZ) = U*diag(1/(b*s+a))*U'
    #uses singular value decomposition of ZKZ for stability when ZKZ is low
    #rank

    if(is.null(invert_aI_bZKZ)) {
      ZKZt = Z %*% K %*% t(Z)
      # identify observations
      if(verbose) print("svd's for missing data")
      Y_obs = apply(!Y_missing,2,which)
      unique_Y_obs = unique(c(list(1:nrow(Y_missing)),Y_obs))
      unique_Y_obs_str = lapply(unique_Y_obs,paste,collapse='')
      Y_obs_index = sapply(Y_obs,function(x) which(unique_Y_obs_str == paste(x,collapse='')))
      invert_aI_bZKZ = lapply(unique_Y_obs,function(x) {
        if(length(x) == 0){
          return(list(
            Ut = as(Matrix(matrix(0,0,0),sparse=T),'dgCMatrix'),
            s = c()
          ))
        }
        result = svd(ZKZt[x,x])
        return(list(
          Ut = t(as(Matrix(result$u,sparse=T),'dgCMatrix')),
          s = result$d
        ))
      })
    }

    #genetic effect variances of factor traits
    # diagonalizing a*Z  '*Z   + b*Kinv for fast inversion
    #diagonalize mixed model equations for fast inversion:
    # inv(a*Z  '*Z   + b*Kinv) = U*diag(1./(a.*s1+b.*s2))*U'
    #similar to fixed effects + random effects 1 above, but no fixed effects.

    if(is.null(invert_aZZt_Kinv)) {
        ZtZ = crossprod(Z)
        svd_ZZt = svd(ZtZ)
        ZZt_sqrt = t(sweep(svd_ZZt$u,2,sqrt(svd_ZZt$d),'*'))

        GSVD_R = function(K,B){
          K_invB = t(solve(t(B),t(K)))
          svd_K_invB = svd(K_invB)
          d = svd_K_invB$d
          U = svd_K_invB$u
          V = svd_K_invB$v
          norm_factor = sqrt(1+d^2)
          c = d/norm_factor
          s = 1/norm_factor
          X = sweep(t(B) %*% V,2,norm_factor,'*')

          return(list(U=svd_K_invB$u, V = svd_K_invB$v,
                          X = X,c=c,s=s))
        }
        result = GSVD_R(ZZt_sqrt,as.matrix(chol_Kinv))

        invert_aZZt_Kinv = list(
            U = as(drop0(Matrix(t(solve(result$X)),sparse=T),tol = run_parameters$drop0_tol),'dgCMatrix'),
            # U = t(solve(result$X)),
    			s1 = result$c^2,
    			s2 = result$s^2
    		)
    }

# ----------------------------- #
# ----Save run parameters------ #
# ----------------------------- #

    run_variables = c(run_variables,list(
      chol_Kinv        = chol_Kinv,
      Y_obs                  = Y_obs,
      Y_obs_index            = Y_obs_index,
			invert_aI_bZKZ   = invert_aI_bZKZ,
			invert_aZZt_Kinv = invert_aZZt_Kinv
    ))

    RNG = list(
    	Random.seed = .Random.seed,
    	RNGkind = RNGkind()
    )

    BSFG_state$run_variables = run_variables
    BSFG_state$RNG = RNG

    return(BSFG_state)
}
