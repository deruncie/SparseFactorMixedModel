sample_MME_single_diagA = function(y, W, C, RinvSqW, prior_mean,prior_prec,Cholesky_R,chol_R,R_Perm,tot_Y_prec) {
	# R is aZAZ + bI
	# 	 then form chol_R
	# 	 check whether solve(t(chol_Rinv),theta) or chol_R %*% theta is better.
	# G is diagonal (fixed effect) - prior prec
	# prior_prec must be > 0
	n = length(y)
	n_theta = length(prior_prec)
	theta_star = prior_mean + rnorm(n_theta)/sqrt(prior_prec)
	e_star = (chol_R) %*% rnorm(n) / sqrt(tot_Y_prec)
	W_theta_star = W %*% theta_star
	y_resid = y - W_theta_star - e_star@x
	if(!is.null(R_Perm)) {
		y_resid_p = R_Perm %*% y_resid
	} else{
		y_resid_p = y_resid
	}

	WtRinvy = crossprod(RinvSqW, solve(Cholesky_R,y_resid_p,'L')) * tot_Y_prec

	theta_tilda = solve(C,WtRinvy)

	theta = theta_tilda@x + theta_star
	theta
}
sample_MME_single_diagA = compiler::cmpfun(sample_MME_single_diagA)

sample_MME_multiple_diagR = function(Y,W,Cholesky_C,pe,chol_A_inv,tot_Y_prec){
	n = nrow(Y)
	p = ncol(Y)
	n_theta = nrow(chol_A_inv)
	theta_star = solve(chol_A_inv,matrix(rnorm(n_theta*p),ncol = p))
	# theta_star = matrix(theta_star@x,ncol = p) + prior_mean
	e_star = matrix(rnorm(n*p)/sqrt(pe),ncol = p,byrow=T)
	W_theta_star = W %*% theta_star
	Y_resid = Y - W_theta_star@x - e_star
	WtRiy = crossprod(W,sweep(Y_resid,2,pe,'*'))

	theta_tilda = solve(Cholesky_C,WtRiy)

	theta = sweep(matrix(theta_tilda@x,nc=p),2,tot_Y_prec,'/') + theta_star@x

	return(theta)
}
sample_MME_multiple_diagR = compiler::cmpfun(sample_MME_multiple_diagR)

sample_MME_fixedEffects = function(Y,W,Sigma_Choleskys, Sigma_Perm, h2s_index, tot_Y_prec, prior_mean, prior_prec,ncores){
	require(parallel)
	# using method described in MCMC Course notes
	p = ncol(Y)
	n = nrow(Y)
	b = ncol(W)
	Wp = W
	Yp = Y
	if(!is.null(Sigma_Perm)) {
		Wp = Sigma_Perm %*% W
		Yp = Sigma_Perm %*% Y
		Yp = matrix(Yp,nrow(Y))
	}
	res = mclapply(1:p,function(j) {
		Cholesky_R = Sigma_Choleskys[[h2s_index[j]]]$Cholesky_Sigma
		chol_R = Sigma_Choleskys[[h2s_index[j]]]$chol_Sigma
		RinvSqW = solve(Cholesky_R,Wp,'L')
		C = crossprod(RinvSqW) * tot_Y_prec[j]
		diag(C) = diag(C) + prior_prec[,j]
		theta_j = sample_MME_single_diagA(Y[,j], W, C, RinvSqW, prior_mean[,j],prior_prec[,j],Cholesky_R,chol_R,Sigma_Perm,tot_Y_prec[j])
		theta_j
	},mc.cores = ncores)
	res = do.call(cbind,res)
	res
}
sample_MME_fixedEffects = compiler::cmpfun(sample_MME_fixedEffects)

sample_MME_ZAZts = function(Y, W, tot_Y_prec, randomEffect_C_Choleskys, h2s, h2s_index, chol_Ai_mats,ncores){
	# using method described in MCMC Course notes
	require(parallel)
	Y = as.matrix(Y)
	p = ncol(Y)
	n = nrow(Y)
	b = ncol(W)
	pes = tot_Y_prec / (1-colSums(h2s))

	unique_h2s = unique(h2s_index)
	unique_h2s_index = sapply(unique_h2s,function(x) which(h2s_index == x))	
	thetas = mclapply(seq_along(unique_h2s),function(j){
		Cholesky_C = randomEffect_C_Choleskys[[unique_h2s[j]]]$Cholesky_C
		chol_A_inv = randomEffect_C_Choleskys[[unique_h2s[j]]]$chol_A_inv * sqrt(tot_Y_prec[j])
		traits_j = unique_h2s_index[[j]]
		sample_MME_multiple_diagR(as.matrix(Y[,traits_j,drop=F]), W, Cholesky_C, pes[traits_j], chol_A_inv,tot_Y_prec[traits_j])		
	},mc.cores = ncores)
	theta = do.call(cbind,thetas)
	theta = theta[,order(unlist(unique_h2s_index))]
	theta
}
sample_MME_ZAZts = compiler::cmpfun(sample_MME_ZAZts)

sample_tot_prec = function(Y, tot_Y_prec_shape, tot_Y_prec_rate, Sigma_Choleskys,Sigma_Perm, h2s_index,ncores){
	n = nrow(Y)
	p = ncol(Y)

	Yp = Y
	if(!is.null(Sigma_Perm)){
		Yp = Sigma_Perm %*% Y
		Yp = matrix(Yp,nrow(Y))
	}
	unique_h2s = unique(h2s_index)
	unique_h2s_index = sapply(unique_h2s,function(x) which(h2s_index == x))	
	scores = mclapply(seq_along(unique_h2s),function(j){
		h2_index = unique_h2s[j]
		Cholesky_Sigma = Sigma_Choleskys[[unique_h2s[j]]]$Cholesky_Sigma
		colSums(solve(Cholesky_Sigma,Yp[,unique_h2s_index[[j]]],'L')^2)
	},mc.cores = ncores)
	scores = unlist(scores)[order(unlist(unique_h2s_index))]
	rgamma(p,shape = tot_Y_prec_shape + n/2, rate = tot_Y_prec_rate + 1/2*scores)
}
sample_tot_prec = compiler::cmpfun(sample_tot_prec)

log_prob_h2 = function(i, y_p, Sigma_Choleskys, n, tot_Y_prec, discrete_priors){
	Cholesky_Sigma = Sigma_Choleskys[[i]]$Cholesky_Sigma
	log_det_Sigma = Sigma_Choleskys[[i]]$log_det
	scores_2 = tot_Y_prec*colSums(solve(Cholesky_Sigma,y_p,'L')^2)

	log_p = -n/2 * log(2*pi) - 1/2*(log_det_Sigma - n*log(tot_Y_prec)) - 1/2 * scores_2 + log(discrete_priors)
	log_p
}
log_prob_h2 = compiler::cmpfun(log_prob_h2)

sample_h2s_discrete_MH = function(Y,tot_Y_prec, Sigma_Choleskys,Sigma_Perm,discrete_priors,h2_divisions,h2_index,step_size,ncores){
	# while this works, it is much slower than doing the full scan over all traits, at least for multiple traits
	# testing with p=100, solving the whole set takes ~4-5x solving just 1. And this method requires doing each trait separately
	# both methods are similarly easy to multiplex, so no advantage there either.
	n = nrow(Y)
	p = ncol(Y)
	discrete_bins = length(discrete_priors)

	Yp = Y
	if(!is.null(Sigma_Perm)){
		Yp = Sigma_Perm %*% Y
		Yp = matrix(Yp,nrow(Y))
	}

	h2s_index = mclapply(1:p,function(j) {
			## steps:
			# 1) calculate log_prob at current state
			# 2) propose new state
			# 3) calc log_prob at new state
			# 4) calc proposal prob from old to new
			# 5) calc proposal prob from new to old
			old_state <- h2_index[j]
			old_log_p <- log_prob_h2(old_state,Yp[,j],Sigma_Choleskys,n,tot_Y_prec[j],discrete_priors[old_state])

			candidate_new_states <- which(colSums((h2_divisions - c(h2_divisions[,old_state]))^2) < step_size)
			proposed_state <- sample(candidate_new_states,1)

			new_log_p <- log_prob_h2(proposed_state,Yp[,j],Sigma_Choleskys,n,tot_Y_prec[j],discrete_priors[proposed_state])
			candidate_state_from_new_state <- which(colSums((h2_divisions - c(h2_divisions[,proposed_state]))^2) < step_size)

			forward_prob = 1/length(candidate_new_states)
			if(old_state %in% candidate_state_from_new_state){
				back_prob = 1/length(candidate_state_from_new_state)
			} else{
				back_prob = 0
			}

			log_MH_ratio = new_log_p - old_log_p + log(forward_prob) - log(back_prob)

			if(log(runif(1)) < log_MH_ratio) {
				return(proposed_state)
			} else {
				return(old_state)
			}
	},mc.cores = ncores)
	do.call(c,h2s_index) 
}
sample_h2s_discrete_MH = compiler::cmpfun(sample_h2s_discrete_MH)

sample_h2s_discrete = function(Y,tot_Y_prec, Sigma_Choleskys,Sigma_Perm,discrete_priors,ncores){
	n = nrow(Y)
	p = ncol(Y)
	discrete_bins = length(discrete_priors)

	Yp = Y
	if(!is.null(Sigma_Perm)){
		Yp = Sigma_Perm %*% Y
		Yp = matrix(Yp,nrow(Y))
	}
	log_ps = mclapply(1:discrete_bins,function(i) {
		Cholesky_Sigma = Sigma_Choleskys[[i]]$Cholesky_Sigma
		log_det_Sigma = Sigma_Choleskys[[i]]$log_det
		scores_2 = tot_Y_prec*colSums(solve(Cholesky_Sigma,Yp,'L')^2)

		log_ps = -n/2 * log(2*pi) - 1/2*(log_det_Sigma - n*log(tot_Y_prec)) - 1/2 * scores_2 + log(discrete_priors[i])
		log_ps
	},mc.cores = ncores)
	if(length(log_ps) == 1) {
		log_ps = matrix(log_ps[[1]],nrow = p)
	} else{
		log_ps = do.call(cbind,log_ps)
	}
	h2s_index = sapply(1:p,function(j) {
		max_row = max(log_ps[j,])
		norm_factor = max_row+log(sum(exp(log_ps[j,]-max_row)))
		ps_j = exp(log_ps[j,] - norm_factor)
		sum(runif(1)>cumsum(ps_j))+1
	})
	return(h2s_index)
}
sample_h2s_discrete = compiler::cmpfun(sample_h2s_discrete)



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
		Posterior$resid_Y_prec[,sp_num] = tot_Y_prec/(1-colSums(resid_h2))
		Posterior$E_a_prec[,sp_num]     = tot_Y_prec/colSums(resid_h2)

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

    p = nrow(BSFG_state$current_state$Lambda)
    b = max(0,ncol(BSFG_state$current_state$B))
    n = nrow(BSFG_state$current_state$F)
    r = nrow(BSFG_state$current_state$E_a)
    n_RE = ncol(BSFG_state$current_state$resid_h2)
    
    Posterior = list(
		    Lambda        = matrix(0,nr=0,nc=0),
		    F_a           = matrix(0,nr=0,nc=0),
		    F             = matrix(0,nr=0,nc=0),
		    delta         = matrix(0,nr=0,nc=0),
            tot_F_prec    = matrix(0,nr=0,nc=0),
		    F_h2          = matrix(0,nr=0,nc=0),
		    tot_Y_prec    = matrix(0,nr = p,nc = 0),
		    resid_h2      = matrix(0,nr = p*n_RE,nc = 0),
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
				F_h2_index = c(F_h2_index,sample(1:ncol(h2_divisions),1))
				F_h2 = h2_divisions[,F_h2_index,drop=FALSE]
				tot_F_prec[k] = 1
				F_a = cbind(F_a,rnorm(r,0,sqrt(sum(F_h2[,k]))))
				F = cbind(F,rnorm(n,as.matrix(Z_all %*% F_a[,k]),sqrt(1-sum(F_h2[,k]))))
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
				F_h2 = F_h2[,nonred,drop=FALSE]
				F_h2_index = F_h2_index[nonred]
				tot_F_prec = tot_F_prec[nonred]
				F_a = F_a[,nonred]
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
		if(is.null(dim(current_state[[param]]))){
			current_state[[param]] = current_state[[param]][factor_order]
		} else if(dim(current_state[[param]])[2] == 1) {
			current_state[[param]] = current_state[[param]][factor_order,,drop=F]			
		} else {
			current_state[[param]] = current_state[[param]][,factor_order]			
		}
	}
	current_state$delta = c(current_state$tauh[1],current_state$tauh[-1]/current_state$tauh[-length(current_state$tauh)])
	BSFG_state$current_state = current_state

	# reorder Posterior
	Posterior = BSFG_state$Posterior
	if(ncol(Posterior$Lambda) == 0) return(BSFG_state)

	p = BSFG_state$run_parameters$setup$p
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
