save_posterior_samples_ipx = function( sp_num, current_state, Posterior) {
	# Posterior,Lambda,F,F_a,B,W,E_a,delta,F_h2,resid_Y_prec,E_a_prec,W_prec) {
	# save posteriors. Full traces are kept of the more interesting parameters.
	# Only the posterior means are kept of less interesting parameters. These
	# should be correctly calculated over several re-starts of the sampler.

	Posterior = with(current_state, {
		# transform variables so that the variance of each column of F is 1.
		F_var = 1/F_a_prec + 1/F_e_prec
		# F_h2 = F_e_prec / (F_e_prec + F_a_prec)
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
		Posterior$resid_Y_prec[,sp_num] = resid_Y_prec
		Posterior$E_a_prec[,sp_num]     = E_a_prec
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
		    F_h2          = matrix(0,nr=0,nc=0),
            tot_Y_prec    = matrix(0,nr=0,nc=0),
            tot_F_prec    = matrix(0,nr=0,nc=0),
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