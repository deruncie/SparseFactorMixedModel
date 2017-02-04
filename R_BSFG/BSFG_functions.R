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
		if(is.null(dim(current_state[[param]]))){
			current_state[[param]] = current_state[[param]][factor_order]
		} else if(dim(current_state[[param]])[2] == 1) {
			current_state[[param]] = current_state[[param]][factor_order,]			
		} else {
			current_state[[param]] = current_state[[param]][,factor_order]			
		}
	}
	current_state$delta = c(current_state$tauh[1],current_state$tauh[-1]/current_state$tauh[-length(current_state$tauh)])
	BSFG_state$current_state = current_state

	# reorder Posterior
	Posterior = BSFG_state$Posterior
	if(ncol(Posterior$Lambda) == 0) return(BSFG_state)

	p = nrow(Lambda)
	k = dim(Posterior$Lambda)[1]/p
	if(length(factor_order) < k) factor_order = c(factor_order,seq(length(factor_order)+1,k))
	
	for(param in reorder_params){
		if(!param %in% names(Posterior)) next
		n = dim(Posterior[[param]])[1]/k
		index = matrix(1:(n*k),nrow = n)
		Posterior[[param]] = Posterior[[param]][c(index[,factor_order]),]
	}
	BSFG_state$Posterior = Posterior

	return(BSFG_state)
}

save_posterior_samples = function( sp_num, current_state, Posterior) {
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
		Posterior$resid_Y_prec[,sp_num] = tot_Y_prec/(1-colSums(matrix(resid_h2,ncol = nrow(Lambda))))
		Posterior$E_a_prec[,sp_num]     = tot_Y_prec/colSums(matrix(resid_h2,ncol = nrow(Lambda)))

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
