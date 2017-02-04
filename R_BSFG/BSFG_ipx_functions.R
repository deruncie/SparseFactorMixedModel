save_posterior_samples_ipx = function( sp_num, current_state, Posterior) {
	# Posterior,Lambda,F,F_a,B,W,E_a,delta,F_h2,resid_Y_prec,E_a_prec,W_prec) {
	# save posteriors. Full traces are kept of the more interesting parameters.
	# Only the posterior means are kept of less interesting parameters. These
	# should be correctly calculated over several re-starts of the sampler.

	current_state = within(current_state,{
		# transform variables so that the variance of each column of F is 1.
		F_var = 1/tot_F_prec
		F_a = sweep(F_a,2,sqrt(F_var),'/')
		F = sweep(F,2,sqrt(F_var),'/')
		Lambda = sweep(Lambda,2,sqrt(F_var),'*')		
	})

	Posterior = with(current_state, {
		sp = ncol(Posterior$Lambda)

		for(param in c('Lambda','F','F_a','delta','F_h2','tot_F_prec')){
			if(length(current_state[[param]]) > nrow(Posterior[[param]])){
				Posterior[[param]] = rbind(Posterior[[param]],matrix(0,nr = length(current_state[[param]]) - nrow(Posterior[[param]]),nc = sp))
			}
			if(length(current_state[[param]]) < nrow(Posterior[[param]])){
				Posterior[[param]] = Posterior[[param]][1:length(current_state[[param]]),]
			}
			Posterior[[param]][,sp_num] = c(current_state[[param]])
		}

		Posterior$tot_Y_prec[,sp_num]   = tot_Y_prec
		Posterior$resid_h2[,sp_num]     = c(resid_h2)
		Posterior$resid_Y_prec[,sp_num] = tot_Y_prec/(1-resid_h2)
		Posterior$E_a_prec[,sp_num]     = tot_Y_prec/resid_h2

		# save B,U,W
		Posterior$B   = (Posterior$B*(sp_num-1) + B)/sp_num
		Posterior$E_a = (Posterior$E_a*(sp_num-1) + E_a)/sp_num
		Posterior
	})
	return(Posterior)
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
		    B             = matrix(0,nr = b,nc = p),
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
				F = cbind(F,rnorm(n,as.matrix(Z %*% F_a[,k]),sqrt(1-F_h2[k])))
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
