save_posterior_samples_ipx = function( sp_num, current_state, Posterior) {
	# Posterior,Lambda,F,F_a,B,W,E_a,delta,F_h2,resid_Y_prec,E_a_prec,W_prec) {
	# save posteriors. Full traces are kept of the more interesting parameters.
	# Only the posterior means are kept of less interesting parameters. These
	# should be correctly calculated over several re-starts of the sampler.

	Posterior = with(current_state, {
		# transform variables so that the variance of each column of F is 1.
		F_var = 1/tot_F_prec
		F_a = sweep(F_a,2,sqrt(F_var),'/')
		F = sweep(F,2,sqrt(F_var),'/')
		Lambda = sweep(Lambda,2,sqrt(F_var),'*')
		sp = ncol(Posterior$Lambda)
		#save factor samples
		if(length(Lambda) > nrow(Posterior$Lambda)){
			# expand factor sample matrix if necessary
			Posterior$Lambda = rbind(Posterior$Lambda, 	matrix(0,nr = length(Lambda)-nrow(Posterior$Lambda),nc = sp))
			Posterior$F      = rbind(Posterior$F, 	   	matrix(0,nr = length(F)     -nrow(Posterior$F),		nc = sp))
			Posterior$F_a    = rbind(Posterior$F_a, 	matrix(0,nr = length(F_a) 	-nrow(Posterior$F_a),	nc = sp))
			Posterior$delta  = rbind(Posterior$delta, 	matrix(0,nr = length(delta) -nrow(Posterior$delta),	nc = sp))
			Posterior$F_h2   = rbind(Posterior$F_h2, 	matrix(0,nr = length(F_h2) 	-nrow(Posterior$F_h2),	nc = sp))
			Posterior$tot_F_prec   = rbind(Posterior$tot_F_prec, 	matrix(0,nr = length(tot_F_prec) 	-nrow(Posterior$tot_F_prec),	nc = sp))
		}
		Posterior$Lambda[1:length(Lambda),sp_num] = c(Lambda)
		Posterior$F[1:length(F),sp_num]     = c(F)
		Posterior$F_a[1:length(F_a),sp_num] = c(F_a)
		Posterior$delta[1:length(delta),sp_num] = delta
		Posterior$F_h2[1:length(F_h2),sp_num] = F_h2
		Posterior$tot_F_prec[1:length(tot_F_prec),sp_num] = tot_F_prec

		Posterior$tot_Y_prec[,sp_num] = tot_Y_prec
		Posterior$resid_h2[,sp_num] = resid_h2
		Posterior$resid_Y_prec[,sp_num] = tot_Y_prec/(1-resid_h2)
		Posterior$E_a_prec[,sp_num]     = tot_Y_prec/resid_h2
		Posterior$W_prec[,sp_num]       = W_prec

		# save B,U,W
		Posterior$B   = (Posterior$B*(sp_num-1) + B)/sp_num
		Posterior$E_a = (Posterior$E_a*(sp_num-1) + E_a)/sp_num
		Posterior$W   = (Posterior$W*(sp_num-1) + W)/sp_num
		Posterior
	})
	return(Posterior)
}
             


sample_delta_ipx = function( delta,tauh,Lambda_prec,delta_1_shape,delta_1_rate,delta_2_shape,delta_2_rate,Lambda2,times = 1 ) {
	#sample delta and tauh parameters that control the magnitudes of higher
	#index factor loadings.

	# matlab = readMat('sample_delta_data.mat')
	# for(i in 1:10) names(matlab) = sub('.','_',names(matlab),fixed=T)
	# delta = matlab$delta
	# tauh = matlab$tauh
	# Lambda_prec = matlab$Lambda_prec
	# delta_1_shape = matlab$delta_1_shape
	# delta_1_rate = matlab$delta_1_rate
	# delta_2_shape = matlab$delta_2_shape
	# delta_2_rate = matlab$delta_2_rate
	# Lambda2	 = matlab$Lambda2	

	k = length(tauh)
	mat = Lambda_prec * Lambda2
	n_genes = nrow(mat)

	for(i in 1:times) {

		shape = delta_1_shape + 0.5*n_genes*k
		rate = delta_1_rate + 0.5*(1/delta[1])*sum(tauh*colSums(mat))
		delta[1] = rgamma(1,shape = shape,rate = rate)
		tauh = cumprod(delta)

		for(h in 2:(k-1)) {
			shape = delta_2_shape + 0.5*n_genes*(k-h+1)
			if(h<k){
				rate = delta_2_rate + 0.5*(1/delta[h])*sum(tauh[h:k]*colSums(mat[,h:k]))
			} else{
				rate = delta_2_rate + 0.5*(1/delta[h])*sum(tauh[h:k]*sum(mat[,h:k]))    	
			}
			delta[h] = rgamma(1,shape = shape, rate = rate)
			tauh = cumprod(delta)
		}
	}
	
	return(delta)
}


clear_Posterior = function(BSFG_state) {
	# resets Posterior samples if burnin was not sufficient
	Posterior = BSFG_state$Posterior
	run_parameters = BSFG_state$run_parameters

	if(!is.null(ncol(Posterior$Lambda))) {    
    	run_parameters$burn = run_parameters$burn + run_parameters$thin*ncol(Posterior$Lambda)
    }

    p = nrow(Posterior$resid_Y_prec)
    b = nrow(Posterior$B)
    n = nrow(Posterior$W)
    r = nrow(Posterior$E_a)
    r2 = nrow(Posterior$W)
    
    Posterior = list(
		    Lambda        = matrix(0,nr=0,nc=0),
		    F_a           = matrix(0,nr=0,nc=0),
		    F             = matrix(0,nr=0,nc=0),
		    delta         = matrix(0,nr=0,nc=0),
            tot_F_prec    = matrix(0,nr=0,nc=0),
		    F_h2          = matrix(0,nr=0,nc=0),
		    tot_Y_prec    = matrix(0,nr = p,nc = 0),
		    resid_h2      = matrix(0,nr = p,nc = 0),
		    resid_Y_prec  = matrix(0,nr = p,nc = 0),
		    E_a_prec      = matrix(0,nr = p,nc = 0),
		    W_prec        = matrix(0,nr = p,nc = 0),
		    B             = matrix(0,nr = b,nc = p),
		    W             = matrix(0,nr = r2,nc = p),
		    E_a           = matrix(0,nr = r,nc = p)
    	)

    BSFG_state$Posterior = Posterior
    BSFG_state$run_parameters = run_parameters
    return(BSFG_state)

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
				k=k+1
				Lambda_prec = cbind(Lambda_prec,rgamma(p,shape = Lambda_df/2, rate = Lambda_df/2))
				delta[k] = rgamma(1,shape = delta_2_shape,rate = delta_2_rate)
				tauh = cumprod(delta)
				Plam = sweep(Lambda_prec,2,tauh,'*')
				Lambda = cbind(Lambda,rnorm(p,0,sqrt(1/Plam[,k])))
				F_h2[k] = runif(1)
				tot_F_prec[k] = 1
				F_a = cbind(F_a,rnorm(r,0,sqrt(F_h2[k])))
				F = cbind(F,rnorm(n,as.matrix(Z_1 %*% F_a[,k]),sqrt(1-F_h2[k])))
			} else if(num > 0) { # drop redundant columns
				nonred = which(vec == 0) # non-redundant loadings columns
				while(length(nonred) < 2) {
					nonred = c(nonred,which(vec != 0)[1])
					vec[nonred[length(nonred)]] = 0
				} 
				k = length(nonred)
				Lambda = Lambda[,nonred]
				Lambda_prec = Lambda_prec[,nonred]
				F = F[,nonred]
				for(red in which(vec == 1)){
					if(red == k) next
					# combine deltas so that the shrinkage of kept columns doesnt
					# decrease after dropping redundant columns
					delta[red+1] = delta[red+1]*delta[red]
				}
				delta = delta[nonred]
				tauh = cumprod(delta)
				Plam = sweep(Lambda_prec,2,tauh,'*')
				F_h2 = F_h2[nonred]
				tot_F_prec = tot_F_prec[nonred]
				F_a = F_a[,nonred]
			}
		}
	}))
	current_state = current_state[current_state_members]

	return(current_state)
}
