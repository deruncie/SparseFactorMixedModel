load_simulation_data = function(file = NULL){
    if(is.null(file)){
      if(file.exists('../setup.RData')) {
        load('../setup.RData')
        for(i in 1:10) names(setup) = sub('.','_',names(setup),fixed=T)
      } else{
        if (!requireNamespace("R.matlab", quietly = TRUE)) {
          stop("R.matlab needed to load setup.mat. Please install it.",
               call. = FALSE)
        }
        setup = readMat('../setup.mat')
        for(i in 1:10) names(setup) = sub('.','_',names(setup),fixed=T)
      }
    } else{
      load(file)
    }
    Y = setup$Y
    A = setup$A
    r = dim(A)[1]
    n = nrow(Y)
    data = data.frame(Group = gl(r,n/r))
    rownames(A) = data$Group
    if(is.null(colnames(Y))) colnames(Y) = paste('Trait',1:ncol(Y),sep='_')
    return(list(Y = Y, data = data, A_mats = list(Group = A),setup = setup))
}


update_k = function( current_state, priors,run_parameters,data_matrices) {
# adapt the number of factors by dropping factors with only small loadings
# if they exist, or adding new factors sampled from the prior if all factors
# appear important. The rate of adaptation decreases through the chain,
# controlled by b0 and b1. Should work correctly over continuations of
# previously stopped chains.

	current_state_members = names(current_state)
	current_state = with(c(priors,run_parameters,data_matrices),within(current_state, {
		i = nrun
		n = nrow(F)
		k = ncol(F)
		r = nrow(F_a)
		p = nrow(Lambda)
		gene_rows = 1:p

		prob = 1/exp(b0 + b1*i)                # probability of adapting
		uu = runif(1)
		lind = colMeans(abs(Lambda) < epsilon)    # proportion of elements in each column less than eps in magnitude
		vec = lind >= prop
		num = sum(vec)       # number of redundant columns

		if(uu < prob && i>200){
			if(i > 20 && num == 0 && all(lind < 0.995) && k < 2*p) { #add a column
				k             = k+1
				Lambda_prec   = cbind(Lambda_prec,rgamma(p,shape = Lambda_df/2, rate = Lambda_df/2))
				delta         = cbind(delta,rgamma(1,shape = delta_2_shape,rate = delta_2_rate))
				tauh          = matrix(cumprod(delta),nrow = 1)
				Plam          = sweep(Lambda_prec,2,tauh,'*')
				Lambda        = cbind(Lambda,rnorm(p,0,sqrt(1/Plam[,k])))
				F_h2_index    = c(F_h2_index,sample(1:ncol(h2s_matrix),1))
				F_h2          = h2s_matrix[,F_h2_index,drop=FALSE]
				tot_F_prec    = cbind(tot_F_prec,1)
				F_a           = cbind(F_a,rnorm(r,0,sqrt(sum(F_h2[,k]))))
				F             = cbind(F,rnorm(n,as.matrix(Z %*% F_a[,k]),sqrt(1-sum(F_h2[,k]))))
			} else if(num > 0) { # drop redundant columns
				nonred = which(vec == 0) # non-redundant loadings columns
				while(length(nonred) < 2) {
					nonred = c(nonred,which(vec != 0)[1])
					vec[nonred[length(nonred)]] = 0
				}
				k = length(nonred)
				Lambda = Lambda[,nonred,drop=FALSE]
				Lambda_prec = Lambda_prec[,nonred,drop=FALSE]
				F = F[,nonred,drop=FALSE]
				for(red in which(vec == 1)){
					if(red == length(vec)) next
					# combine deltas so that the shrinkage of kept columns doesnt
					# decrease after dropping redundant columns
					delta[red+1] = delta[red+1]*delta[red]
				}
				if(is.null(dim(delta))) recover()
				delta = delta[,nonred,drop=FALSE]
				tauh = matrix(cumprod(delta),nrow=1)
				Plam = sweep(Lambda_prec,2,tauh,'*')
				F_h2 = F_h2[,nonred,drop=FALSE]
				F_h2_index = F_h2_index[nonred]
				tot_F_prec = tot_F_prec[,nonred,drop=FALSE]
				F_a = F_a[,nonred,drop=FALSE]
			}
		}
	}))
	current_state = current_state[current_state_members]

	return(current_state)
}

reorder_factors = function(BSFG_state){
	# re-orders factors in decreasing size of Lambda %*% F
	# based on current state
	# also re-orders Posterior

	current_state = BSFG_state$current_state

	Lambda = current_state$Lambda
	F = current_state$F

	# size is sum lambda_ij^2 * var(F_i)
	sizes = colSums(Lambda^2) * colMeans(F^2)
	factor_order = order(sizes,decreasing=T)

	reorder_params = c('Lambda','Lambda_prec','Plam',
						'delta','tauh',
						'F','F_a','F_h2','F_a_prec','F_e_prec','tot_F_prec'
						)

	# reorder currrent state
	for(param in reorder_params){
		if(! param %in% names(current_state)) next
		current_state[[param]] = current_state[[param]][,factor_order,drop=F]
	}
	current_state$delta = matrix(c(current_state$tauh[1],current_state$tauh[-1]/current_state$tauh[-length(current_state$tauh)]),nrow=1)
	BSFG_state$current_state = current_state

	# reorder Posterior
	Posterior = BSFG_state$Posterior

	for(param in reorder_params){
		if(! param %in% names(Posterior)) next
		if(dim(Posterior[[param]])[3] == 0) next
		Posterior[[param]] = Posterior[[param]][,factor_order,,drop=FALSE]
	}

	BSFG_state$Posterior = Posterior

	return(BSFG_state)
}

save_posterior_samples = function( sp_num, current_state, Posterior) {
	require(abind)

	# All parameters in current are matrices.
	# Posterior is a list of arrays.
	# All factor parameters are matrices with ncol == k
	# Posterior arrays are expanded / contracted as the number of factors changes (with update_k)

  Posterior$sp_num = sp_num

	current_state = within(current_state,{
		# transform variables so that the variance of each column of F is 1.
		F_var = 1/tot_F_prec
		F_a = sweep(F_a,2,sqrt(F_var),'/')
		F = sweep(F,2,sqrt(F_var),'/')
		Lambda = sweep(Lambda,2,sqrt(F_var),'*')
	})

	sp = dim(Posterior$Lambda)[3]

	for(param in Posterior$sample_params){
		ncol_current = ncol(current_state[[param]])
		ncol_Posterior = dim(Posterior[[param]])[2]
		if(ncol_current > ncol_Posterior){
			Posterior[[param]] = abind(Posterior[[param]],array(0,dim = c(nrow(current_state[[param]]),ncol_current-ncol_Posterior,sp)),along=2)
		}
		if(ncol_current < ncol_Posterior) {
			Posterior[[param]] = Posterior[[param]][,1:ncol_current,,drop=F]
		}
		Posterior[[param]][,,sp_num] = current_state[[param]]
	}

	for(param in Posterior$posteriorMean_params){
		Posterior[[param]] = Posterior[[param]]*(sp_num - 1) + current_state[[param]]/sp_num
	}

	return(Posterior)
}


initialize_Posterior = function(Posterior,current_state){
	  for(param in Posterior$sample_params){
    	Posterior[[param]] = array(0,dim = c(dim(current_state[[param]]),0))
    }

    for(param in Posterior$posteriorMean_params) {
    	Posterior[[param]] = array(0,dim = dim(current_state[[param]]))
    }
    for(param in Posterior$per_trait_params){
    	dimnames(Posterior[[param]])[[2]] = current_state$traitnames
    }
    dimnames(Posterior$Lambda)[[1]] = current_state$traitnames
    Posterior$sp_num = 0
    Posterior
}

expand_Posterior = function(Posterior,size){
	require(abind)

	for(param in Posterior$sample_params){
		Posterior[[param]] = abind(Posterior[[param]],array(NA,dim = c(dim(Posterior[[param]])[1:2],size)),along = 3)
	}
	Posterior
}

clear_Posterior = function(BSFG_state) {
	# resets Posterior samples if burnin was not sufficient
	Posterior = BSFG_state$Posterior
	run_parameters = BSFG_state$run_parameters

	current_samples = dim(Posterior[[Posterior$sample_params[1]]])[3]

	if(current_samples > 0) {
    	run_parameters$burn = run_parameters$burn + run_parameters$thin*current_samples
    }

    Posterior = initialize_Posterior(Posterior,BSFG_state$current_state)

    BSFG_state$Posterior = Posterior
    BSFG_state$run_parameters = run_parameters
    return(BSFG_state)
}
