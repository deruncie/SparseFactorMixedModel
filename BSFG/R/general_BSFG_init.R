initialize_BSFG.general_BSFG = function(BSFG_state, K_mats = NULL, chol_Ki_mats = NULL,
                                        Sigma_Choleskys = NULL, randomEffect_C_Choleskys = NULL,
									ncores = detectCores(),verbose=T,...){

    Z_matrices = BSFG_state$data_matrices$Z_matrices
    Z          = BSFG_state$data_matrices$Z
    h2s_matrix = BSFG_state$data_matrices$h2s_matrix

    RE_names   = rownames(h2s_matrix)
    n_RE       = length(RE_names)

    run_parameters = BSFG_state$run_parameters
    run_variables  = BSFG_state$run_variables

# ------------------------------------ #
# ----Precalculate ZKZts, chol_Ks ---- #
# ------------------------------------ #

  # setup of symbolic Cholesky of C
  if(!is.null(randomEffect_C_Choleskys)){
    if(length(randomEffect_C_Choleskys) != ncol(h2s_matrix)) stop('wrong number of randomEffect_C_Choleskys provided')
    if(verbose) print('using provided randomEffect_C_Choleskys')
  } else{
    if(verbose) {
      print('Pre-calculating matrices')
      print('creating randomEffects_C')
    }

    ZtZ = as(forceSymmetric(drop0(crossprod(Z),tol = run_parameters$drop0_tol)),'dgCMatrix')
    # randomEffect_C_Choleskys_c = new(randomEffect_C_Cholesky_database,lapply(chol_Ki_mats,function(x) as(x,'dgCMatrix')),h2s_matrix,as(ZtZ,'dgCMatrix'),run_parameters$drop0_tol,1)
    # randomEffect_C_Choleskys = lapply(1:ncol(h2s_matrix),function(i) {
    #   list(chol_C = randomEffect_C_Choleskys_c$get_chol_Ci(i),
    #        chol_K_inv = randomEffect_C_Choleskys_c$get_chol_K_inv_i(i)
    #   )
    # })

    # do calculations in several sets
    group_size = 2*detectCores()
    n_groups = ceiling(ncol(h2s_matrix)/group_size)
    col_groups = tapply(1:ncol(h2s_matrix),gl(n_groups,group_size,ncol(h2s_matrix)),function(x) x)
    randomEffect_C_Choleskys_c_list = list()
    if(verbose) pb = txtProgressBar(min=0,max = n_groups,style=3)
    for(i in 1:length(col_groups)){
      randomEffect_C_Choleskys_c_list[[i]] = new(randomEffect_C_Cholesky_database,lapply(chol_Ki_mats,function(x) as(x,'dgCMatrix')),h2s_matrix[,col_groups[[i]],drop=FALSE],ZtZ,run_parameters$drop0_tol,1)
      if(verbose) setTxtProgressBar(pb, i)
    }
    if(verbose) close(pb)
    randomEffect_C_Choleskys = do.call(c,lapply(1:length(col_groups),function(j) {
      randomEffect_C_Choleskys_c = randomEffect_C_Choleskys_c_list[[j]]
      lapply(1:length(col_groups[[j]]),function(i) {
        list(chol_C = randomEffect_C_Choleskys_c$get_chol_Ci(i),
             chol_K_inv = randomEffect_C_Choleskys_c$get_chol_K_inv_i(i)
        )
      })
    }))
  }

	# Sigma
	if(!is.null(Sigma_Choleskys)) {
	  if(length(Sigma_Choleskys) != ncol(h2s_matrix)) stop('wrong number of Sigma_Choleskys provided')
	  if(verbose) print('using provided Sigma_Choleskys')
	} else{
  	if(verbose) print('creating Sigma_Choleskys')
  	ZKZts = list()
  	for(i in 1:n_RE){
  		ZKZts[[i]] = forceSymmetric(drop0(Z_matrices[[i]] %*% K_mats[[i]] %*% t(Z_matrices[[i]]),tol = run_parameters$drop0_tol))
  	}
  	Sigma_Choleskys_c = new(Sigma_Cholesky_database,lapply(ZKZts,function(x) as(x,'dgCMatrix')),h2s_matrix,run_parameters$drop0_tol,1)
  	Sigma_Choleskys = lapply(1:ncol(h2s_matrix),function(i) {
  	  list(log_det = Sigma_Choleskys_c$get_log_det(i),
  	       chol_Sigma = Sigma_Choleskys_c$get_chol_Sigma(i))
  	})

  	# do calculations in several sets
  	# group_size = 2*detectCores()
  	# n_groups = ceiling(ncol(h2s_matrix)/group_size)
  	# col_groups = tapply(1:ncol(h2s_matrix),gl(n_groups,group_size,ncol(h2s_matrix)),function(x) x)
  	# Sigma_Choleskys_c_list = list()
  	# if(verbose) pb = txtProgressBar(min=0,max = n_groups,style=3)
  	# for(i in 1:length(col_groups)){
  	#   Sigma_Choleskys_c_list[[i]] = new(Sigma_Cholesky_database,lapply(ZKZts,function(x) as(x,'dgCMatrix')),h2s_matrix[,col_groups[[i]],drop=FALSE],run_parameters$drop0_tol,1)
  	#   if(verbose) setTxtProgressBar(pb, i)
  	# }
  	# if(verbose) close(pb)
  	# Sigma_Choleskys = do.call(c,lapply(1:length(col_groups),function(j) {
  	#   Sigma_Choleskys_c = Sigma_Choleskys_c_list[[j]]
  	#   lapply(1:length(col_groups[[j]]),function(i) {
  	#     list(log_det = Sigma_Choleskys_c$get_log_det(i),
  	#          chol_Sigma = Sigma_Choleskys_c$get_chol_Sigma(i))
  	#   })
  	# }))
	}

	if(verbose) print('done')


# ----------------------------- #
# ----Save run parameters------ #
# ----------------------------- #

	run_variables = c(run_variables,list(
	  Sigma_Choleskys          = Sigma_Choleskys,
	  randomEffect_C_Choleskys = randomEffect_C_Choleskys
    ))

    RNG = list(
    	Random.seed = .Random.seed,
    	RNGkind = RNGkind()
    )

    BSFG_state$run_variables = run_variables
    BSFG_state$RNG = RNG

    return(BSFG_state)
}

