
sample_MME_single_diagA = function(y,W,C,RinvSqW,prior_mean,prior_prec,chol_R){
	n = length(y)
	n_theta = length(prior_prec)
	theta_star = prior_mean + rnorm(n_theta,0,1/sqrt(prior_prec))
	e_star = chol_R %*% rnorm(n)
	W_theta_star = W %*% theta_star

	y_resid = y - W_theta_star - e_star@x
	WtRinvy = crossprod(RinvSqW,solve(t(chol_R),y_resid))

	theta_tilda = solve(C,WtRinvy)

	theta = theta_tilda@x + theta_star
	return(theta)
}
sample_MME_single_diagR = function(y,W,C,pe,prior_mean,chol_A){
	n = length(y)
	n_theta = nrow(chol_A)
	theta_star = chol_A %*% rnorm(n_theta)
	theta_star = theta_star@x + prior_mean
	e_star = rnorm(n)/sqrt(pe)
	W_theta_star = W %*% theta_star

	y_resid = y - W_theta_star@x - e_star
	WtRiy = crossprod(W,y_resid*pe)

	theta_tilda = solve(C,WtRiy)

	theta = theta_tilda@x + theta_star
	return(theta)
}

sample_MME_fixedEffects = function(Y,W,Sigmas, h2s_index, tot_Y_prec, prior_mean, prior_prec,ncores){
	require(parallel)
	# using method described in MCMC Course notes
	p = ncol(Y)
	n = nrow(Y)
	b = ncol(W)

	chunkSize = ceiling(p/ncores)
	res = mclapply(1:ceiling(p/chunkSize),function(chunk) {
		cols = 1:chunkSize + (chunk-1)*chunkSize
		cols = cols[cols <= p]
		thetas = do.call(cbind,lapply(cols,function(j) {
			chol_R = Sigmas[[h2s_index[j]]]$chol #/ sqrt(tot_Y_prec[j])
			chol_R@x = chol_R@x / sqrt(tot_Y_prec[j])
			RinvSqW = solve(t(chol_R),W)
			C = crossprod(RinvSqW)
			diag(C) = diag(C) + prior_prec[,j]
			theta_j = sample_MME_single_diagA(Y[,j], W, C, RinvSqW, prior_mean[,j],prior_prec[,j],chol_R)			
			theta_j
		}))
		thetas
	},mc.cores = ncores)
	res = do.call(cbind,res)
	res
}


sample_MME_ZAZts = function(Y, W, tot_Y_prec, prior_mean, randomEffect_Cs, Ai_mats, h2s, h2_index, chol_As,ncores){
	require(parallel)
	# using method described in MCMC Course notes
	p = ncol(Y)
	n = nrow(Y)
	b = ncol(W)
	pes = tot_Y_prec / (1-rowSums(h2s))

	chunkSize = ceiling(p/ncores)
	res = mclapply(1:ceiling(p/chunkSize),function(chunk) {
		cols = 1:chunkSize + (chunk-1)*chunkSize
		cols = cols[cols <= p]
		thetas = do.call(cbind,lapply(cols,function(j) {
			h2s_j = pmax(1e-10,h2s[j,])
			C = randomEffect_Cs[[h2_index[j]]]
			C@x = C@x * tot_Y_prec[j]
			chol_A = do.call(bdiag,lapply(1:ncol(h2s),function(i) {
				chol_Ai = chol_As[[i]]
				chol_Ai@x = chol_Ai@x *sqrt(h2s_j[i]/tot_Y_prec[j])
				chol_Ai
			}))
			theta_j = sample_MME_single_diagR(Y[,j], W, C, pes[j], prior_mean[,j],chol_A)
			theta_j
		}))
		thetas
	},mc.cores = ncores)
	res = do.call(cbind,res)
	res
}

sample_tot_prec = function(Y, tot_Y_prec_shape, tot_Y_prec_rate, Sigmas, h2s_index,ncores){
	n = nrow(Y)
	p = ncol(Y)

	chunkSize = ceiling(p/ncores)
	scores = mclapply(1:ceiling(p/chunkSize),function(chunk) {
		cols = 1:chunkSize + (chunk-1)*chunkSize
		cols = cols[cols <= p]
		sapply(cols,function(j) {
			chol_Sigma = Sigmas[[h2s_index[j]]]$chol
			sum(solve(t(chol_Sigma),Y[,j])^2)
		})
	},mc.cores = ncores)
	scores = do.call(c,scores)
	rgamma(p,shape = tot_Y_prec_shape + n/2, rate = tot_Y_prec_rate + 1/2*scores)
}

sample_h2s_discrete = function(Y,tot_Y_prec, Sigmas,h2_divisions,discrete_priors,ncores){
# 	Y = Y[,1:98]
# 	tot_Y_prec = tot_Y_prec[1:98]
	n = nrow(Y)
	p = ncol(Y)
	discrete_bins = length(discrete_priors)

	chunkSize = ceiling(discrete_bins/ncores)
	log_ps = mclapply(1:ceiling(discrete_bins/chunkSize),function(chunk) {
		rows = 1:chunkSize + (chunk-1)*chunkSize
		rows = rows[rows <= discrete_bins]
		sapply(rows,function(i) {
			chol_R = Sigmas[[i]]$chol
			det_R = Sigmas[[i]]$det
			# scores_2 = tot_Y_prec*diag(crossprod(Y,solve(R,Y)))
			scores_2 = tot_Y_prec*colSums(solve(t(chol_R),Y)^2)

			# Ri = solve(R)
			# scores_2b = tot_Y_prec*diag(crossprod(Y,Ri %*% Y))

			# chol_Ri = chol(ri)
			# scores_2c = tot_Y_prec*diag(crossprod(chol_Ri %*% Y))
			# scores_2d = tot_Y_prec*colSums((chol_Ri %*% Y)^2)
			# chol_R = chol(R)
			# scores_2e = tot_Y_prec*colSums(solve(chol_R, Y)^2)

			# microbenchmark(diag(crossprod(Y,solve(R,Y))),diag(crossprod(Y,Ri %*% Y)),diag(crossprod(chol_Ri %*% Y)),tot_Y_prec*colSums((chol_Ri %*% Y)^2),tot_Y_prec*colSums(solve(chol(R), Y)^2))

			log_ps = -n/2 * log(2*pi) - 1/2*(log(det_R) - n*log(tot_Y_prec)) - 1/2 * scores_2 + log(discrete_priors[i])
			log_ps
		})
	},mc.cores = ncores)
	log_ps = do.call(cbind,log_ps)
	h2s_index = sapply(1:p,function(j) {
		max_row = max(log_ps[j,])
		norm_factor = max_row+log(sum(exp(log_ps[j,]-max_row)))
		ps_j = exp(log_ps[j,] - norm_factor)
		sum(runif(1)>cumsum(ps_j))
	})
	# }))
	return(h2s_index)
}


save_posterior_samples = function( sp_num, current_state, Posterior) {
	# Posterior,Lambda,F,F_a,B,W,E_a,delta,F_h2,resid_Y_prec,E_a_prec,W_prec) {
	# save posteriors. Full traces are kept of the more interesting parameters.
	# Only the posterior means are kept of less interesting parameters. These
	# should be correctly calculated over several re-starts of the sampler.

	Posterior = with(current_state, {
		# transform variables so that the variance of each column of F is 1.
		F_var = 1/tot_F_prec
		F_a = sweep(F_a,2,sqrt(F_var),'/')
		F = sweep(F,2,sqrt(F_var),'/')
		F_h2 = rowSums(F_h2)
		resid_h2 = rowSums(resid_h2)
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

		# save B,U,W
		Posterior$B   = (Posterior$B*(sp_num-1) + B)/sp_num
		Posterior$E_a = (Posterior$E_a*(sp_num-1) + E_a)/sp_num
		Posterior
	})
	return(Posterior)
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
				F_h2[k,] = runif(ncol(F_h2))
				F_h2_index[k] = 1
				tot_F_prec[k] = 1
				F_a = cbind(F_a,rnorm(r,0,sqrt(sum(F_h2[k,]))))
				F = cbind(F,rnorm(n,as.matrix(Z_all %*% F_a[,k]),sqrt(1-sum(F_h2[k,]))))
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
				F_h2 = F_h2[nonred,,drop=FALSE]
				F_h2_index = F_h2_index[nonred]
				tot_F_prec = tot_F_prec[nonred]
				F_a = F_a[,nonred]
			}
		}
	}))
	current_state = current_state[current_state_members]

	return(current_state)
}
