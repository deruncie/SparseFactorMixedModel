#' Sample from Mixed Model equations given a diagonal K matrix
#'
#' Sample from Mixed Model equations given a diagonal K matrix with sparse but non-diagonal R
#'
#' This function draws a sample of the vector \eqn{theta} given the model:
#'     \deqn{y = W\theta + e}
#'     \deqn{\theta \sim N_b(\mu,K)} \deqn{e \sim N(0,R)} The algorithm follows that shown in the
#'     MCMCglmm course notes. This involves solving:
#'     \deqn{\tilde{\theta} = C^{-1}W'R^{-1}(1 - W\theta_{*} - e_{*})}
#'     \deqn{C = W'R^{-1}W + K^{-1}}
#'     \deqn{\theta_{*} \sim N(\mu,K), e_{*} \sim N(0,R)}
#'
#' To speed up the calculations, \eqn{R_hat = R*tot_Eta_prec = P' L L' P} is pre-calculated, and used
#' in place of \eqn{R^{-1}}.
#'
#' @param y vector of length n
#' @param W matrix of size n x b
#' @param Wp \eqn{R_perm * W} (see below)
#' @param prior_mean the vector \eqn{\mu}
#' @param prior_prec the diagonal of the matrix K^{-1}. All elements must be > 0
#' @param Cholesky_R Cholesky decomposition of the sparse matrix R calculated by
#'   \code{Cholesky(R_hat,perm=T,super=T)}. The factorization is  \eqn{P' L L' P}. R_hat here is
#'   \eqn{R*tot_Eta_prec}, i.e., the covariance up to normalization by the total variance.
#' @param chol_R the matrix \eqn{P' L} corresponding to the Cholesky_R decomposition above
#' @param R_Perm either \code{NULL} if \eqn{P} is diagonal, or the \eqn{P} matrix.
#' @param tot_Eta_prec the inverse of the total variance
#'
sample_MME_single_diagK = function(y, W, Wp, prior_mean,prior_prec,Cholesky_R,chol_R,R_Perm,tot_Eta_prec,randn_theta = NULL, randn_e = NULL) {
	n = length(y)
	n_theta = length(prior_prec)
	if(is.null(randn_theta)) {
	  randn_theta = rnorm(n_theta)
	}
	if(is.null(randn_e)) {
	  randn_e = rnorm(n)
	}
	theta_star = prior_mean + randn_theta/sqrt(prior_prec)
	e_star = (chol_R) %*% randn_e / sqrt(tot_Eta_prec)
	W_theta_star = W %*% theta_star
	if(is(W_theta_star,'Matrix')) W_theta_star = W_theta_star@x
	if(is(e_star,'Matrix')) e_star = e_star@x
	y_resid = y - W_theta_star - e_star
	if(!is.null(R_Perm)) {
		y_resid_p = R_Perm %*% y_resid
	} else{
		y_resid_p = y_resid
	}

	RinvSqW = solve(Cholesky_R,Wp,'L')
	WtRinvy = crossprod(RinvSqW, solve(Cholesky_R,y_resid_p,'L')) * tot_Eta_prec

	if(ncol(Wp) < nrow(Wp)) {
  	C = crossprod(RinvSqW) * tot_Eta_prec
  	diag(C) = diag(C) + prior_prec
  	theta_tilda = solve(C,WtRinvy)
	} else{
	  R = tcrossprod(chol_R)/tot_Eta_prec
	  AiU = t(Wp)/prior_prec
	  inner = AiU %*% solve(R + Wp %*% AiU,t(AiU))
	  WtRinvy = WtRinvy@x
	  theta_tilda = WtRinvy/prior_prec - (inner %*% WtRinvy)@x
	}
	if(is(theta_tilda,'Matrix')) theta_tilda = theta_tilda@x

	theta = theta_tilda + theta_star
	theta
}
sample_MME_single_diagK = compiler::cmpfun(sample_MME_single_diagK)

sample_MME_multiple_diagR = function(Y,W,Cholesky_C,pe,chol_K_inv,tot_Eta_prec, randn_theta = NULL, randn_e = NULL){
	n = nrow(Y)
	p = ncol(Y)
	n_theta = nrow(chol_K_inv)
	if(is.null(randn_theta)) {
	  randn_theta = matrix(rnorm(n_theta*p),ncol = p)
	}
	if(is.null(randn_e)) {
	  randn_e = matrix(rnorm(n*p), ncol = p)
	}
	theta_star = solve(chol_K_inv,randn_theta)
	e_star = sweep(randn_e,2,sqrt(pe),'/')
	W_theta_star = W %*% theta_star

	if(is(W_theta_star,'Matrix')) W_theta_star = W_theta_star@x
	if(is(theta_star,'Matrix')) theta_star = theta_star@x

	Y_resid = Y - W_theta_star - e_star
	WtRiy = crossprod(W,sweep(Y_resid,2,pe,'*'))

	theta_tilda = solve(Cholesky_C,WtRiy)

	theta = sweep(matrix(theta_tilda@x,nc=p),2,tot_Eta_prec,'/') + theta_star

	return(theta)
}
sample_MME_multiple_diagR = compiler::cmpfun(sample_MME_multiple_diagR)

sample_MME_fixedEffects = function(Y,W,Sigma_Choleskys, Sigma_Perm, h2s_index, tot_Eta_prec, prior_mean, prior_prec,ncores){
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

	# pre-sample z-scores because draws from parallel processes not consecutive
	randn_theta = matrix(rnorm(b*p),ncol = p)
	randn_e = matrix(rnorm(n*p),ncol = p)

	res = mclapply(1:p,function(j) {
		Cholesky_R = Sigma_Choleskys[[h2s_index[j]]]$Cholesky_Sigma
		chol_R = Sigma_Choleskys[[h2s_index[j]]]$chol_Sigma
		# recover()
		theta_j = sample_MME_single_diagK(Y[,j], W, Wp, prior_mean[,j],prior_prec[,j],Cholesky_R,chol_R,Sigma_Perm,tot_Eta_prec[j], randn_theta[,j],randn_e[,j])
		# chol_R_2 =as(chol_R,'dgCMatrix')
		# theta_j2 = sample_MME_single_diagK_c(Y[,j],W,prior_mean[,j],prior_prec[,j],chol_R_2,tot_Eta_prec[j], randn_theta[,j],randn_e[,j])
		# plot(theta_j,theta_j2);abline(0,1)
		# microbenchmark(sample_MME_single_diagK(Y[,j], W, Wp, prior_mean[,j],prior_prec[,j],Cholesky_R,chol_R,Sigma_Perm,tot_Eta_prec[j], randn_theta[,j],randn_e[,j]),
		#                sample_MME_single_diagK_c(Y[,j],W,prior_mean[,j],prior_prec[,j],chol_R_2,tot_Eta_prec[j], randn_theta[,j],randn_e[,j]))

	},mc.cores = ncores)
	res = do.call(cbind,res)
	res

	# chol_Rs = lapply(Sigma_Choleskys,function(x) as(x$chol_Sigma,'dgCMatrix'))
	# res2 = sample_MME_fixedEffects_c(Y,W,chol_Rs,h2s_index,tot_Eta_prec,prior_mean,prior_prec,randn_theta,randn_e,1)
	#
	# microbenchmark({
	#   res = mclapply(1:p,function(j) {
	#     Cholesky_R = Sigma_Choleskys[[h2s_index[j]]]$Cholesky_Sigma
	#     chol_R = Sigma_Choleskys[[h2s_index[j]]]$chol_Sigma
	#     # recover()
	#     theta_j = sample_MME_single_diagK(Y[,j], W, Wp, prior_mean[,j],prior_prec[,j],Cholesky_R,chol_R,Sigma_Perm,tot_Eta_prec[j], randn_theta[,j],randn_e[,j])
	#     # chol_R_2 =as(chol_R,'dgCMatrix')
	#     # theta_j2 = sample_MME_single_diagK_c(Y[,j],W,prior_mean[,j],prior_prec[,j],chol_R_2,tot_Eta_prec[j], randn_theta[,j],randn_e[,j])
	#     # plot(theta_j,theta_j2);abline(0,1)
	#     # microbenchmark(sample_MME_single_diagK(Y[,j], W, Wp, prior_mean[,j],prior_prec[,j],Cholesky_R,chol_R,Sigma_Perm,tot_Eta_prec[j], randn_theta[,j],randn_e[,j]),
	#     #                sample_MME_single_diagK_c(Y[,j],W,prior_mean[,j],prior_prec[,j],chol_R_2,tot_Eta_prec[j], randn_theta[,j],randn_e[,j]))
	#
	#   },mc.cores = ncores)
	#   res = do.call(cbind,res)
	# },
	# sample_MME_fixedEffects_c(Y,W,chol_Rs,h2s_index,tot_Eta_prec,prior_mean,prior_prec,randn_theta,randn_e,1)
	# )

}
sample_MME_fixedEffects = compiler::cmpfun(sample_MME_fixedEffects)

sample_MME_ZKZts = function(Y, W, tot_Eta_prec, randomEffect_C_Choleskys, h2s, h2s_index, chol_Ki_mats,ncores){
	# using method described in MCMC Course notes
	Y = as.matrix(Y)
	p = ncol(Y)
	n = nrow(Y)
	b = ncol(W)
	pes = tot_Eta_prec / (1-colSums(h2s))

	# pre-sample z-scores because draws from parallel processes not consecutive
	randn_theta = matrix(rnorm(b*p),ncol = p)
	randn_e = matrix(rnorm(n*p),ncol = p)

	unique_h2s = unique(h2s_index)
	unique_h2s_index = lapply(unique_h2s,function(x) which(h2s_index == x))
	thetas = mclapply(seq_along(unique_h2s),function(j){
		Cholesky_C = randomEffect_C_Choleskys[[unique_h2s[j]]]$Cholesky_C
		chol_K_inv = randomEffect_C_Choleskys[[unique_h2s[j]]]$chol_K_inv * sqrt(tot_Eta_prec[j])  # this is wrong!!!!!
		traits_j = unique_h2s_index[[j]]
		# recover()
		res1 = sample_MME_multiple_diagR(as.matrix(Y[,traits_j,drop=F]), W, Cholesky_C, pes[traits_j], chol_K_inv,tot_Eta_prec[traits_j], randn_theta[,traits_j,drop=F], randn_e[,traits_j,drop=F])
		chol_C = as(expand(Cholesky_C)$L,'dgCMatrix')
		chol_K_inv_2 = as(chol_K_inv,'dgCMatrix')
		print(paste('asdf',chol_K_inv_2[11,11]))
		# res2 = sample_MME_multiple_diagR_c(as.matrix(Y[,traits_j,drop=F]), W, chol_C, pes[traits_j], chol_K_inv_2,tot_Eta_prec[traits_j], randn_theta[,traits_j,drop=F], randn_e[,traits_j,drop=F])
		# res3 = sapply(traits_j,function(i) sample_MME_single_diagR_c(as.matrix(Y[,i,drop=F]), W, chol_C, pes[i], chol_K_inv_2,tot_Eta_prec[i], randn_theta[,i,drop=F], randn_e[,i,drop=F]))
		res4 = sample_MME_ZKZts_c(as.matrix(Y[,traits_j,drop=F]), W, tot_Eta_prec[traits_j], chol_Cs, chol_K_invs,h2s[,traits_j,drop=FALSE],h2s_index[traits_j],randn_theta[,traits_j,drop=F], randn_e[,traits_j,drop=F],1)
		# res5 = sapply(traits_j,function(i) sample_MME_ZKZts_c(as.matrix(Y[,i,drop=F]), W, tot_Eta_prec[i], chol_Cs, chol_K_invs,h2s[,i,drop=FALSE],h2s_index[i],randn_theta[,i,drop=F], randn_e[,i,drop=F],1))
		# plot(res3,res4);abline(0,1)
		# plot(res3,res5);abline(0,1)
		# plot(res4,res5);abline(0,1)
		# plot(res1,res2);abline(0,1)
		# plot(res1,res3);abline(0,1)
		# microbenchmark(
		#   sample_MME_multiple_diagR(as.matrix(Y[,traits_j,drop=F]), W, Cholesky_C, pes[traits_j], chol_K_inv,tot_Eta_prec[traits_j], randn_theta[,traits_j,drop=F], randn_e[,traits_j,drop=F]),
		#   sample_MME_multiple_diagR_c(as.matrix(Y[,traits_j,drop=F]), W, chol_C, pes[traits_j], chol_K_inv_2,tot_Eta_prec[traits_j], randn_theta[,traits_j,drop=F], randn_e[,traits_j,drop=F]),
		#   sapply(traits_j,function(j) sample_MME_multiple_diagR_c2(as.matrix(Y[,j,drop=F]), W, chol_C, pes[j], chol_K_inv_2,tot_Eta_prec[j], randn_theta[,j,drop=F], randn_e[,j,drop=F]))
		# )
		# colnames(res3) = traits_j
		# res3
		res1
	},mc.cores = ncores)
	theta = do.call(cbind,thetas)
	theta = theta[,order(unlist(unique_h2s_index))]
	theta

	sapply(chol_Cs[unique_h2s],function(x) x[11,11])


	chol_Cs = lapply(randomEffect_C_Choleskys,function(x) as(expand(x$Cholesky_C)$L,'dgCMatrix'))
	chol_K_invs = lapply(randomEffect_C_Choleskys,function(x) as(x$chol_K_inv,'dgCMatrix'))
	theta2 = sample_MME_ZKZts_c(Y,W,tot_Eta_prec,chol_Cs,chol_K_invs,h2s,h2s_index,randn_theta,randn_e,1)
  plot(theta,theta2);abline(0,1)

	microbenchmark({
	  thetas = mclapply(seq_along(unique_h2s),function(j){
	    Cholesky_C = randomEffect_C_Choleskys[[unique_h2s[j]]]$Cholesky_C
	    chol_K_inv = randomEffect_C_Choleskys[[unique_h2s[j]]]$chol_K_inv * sqrt(tot_Eta_prec[j])
	    traits_j = unique_h2s_index[[j]]
	    # recover()
	    res1 = sample_MME_multiple_diagR(as.matrix(Y[,traits_j,drop=F]), W, Cholesky_C, pes[traits_j], chol_K_inv,tot_Eta_prec[traits_j], randn_theta[,traits_j,drop=F], randn_e[,traits_j,drop=F])
	    # chol_C = as(expand(Cholesky_C)$L,'dgCMatrix')
	    # chol_K_inv_2 = as(chol_K_inv,'dgCMatrix')
	    # res2 = sample_MME_multiple_diagR_c(as.matrix(Y[,traits_j,drop=F]), W, chol_C, pes[traits_j], chol_K_inv_2,tot_Eta_prec[traits_j], randn_theta[,traits_j,drop=F], randn_e[,traits_j,drop=F])
	    # res3 = sapply(traits_j,function(j) sample_MME_multiple_diagR_c2(as.matrix(Y[,j,drop=F]), W, chol_C, pes[j], chol_K_inv_2,tot_Eta_prec[j], randn_theta[,j,drop=F], randn_e[,j,drop=F]))
	    # plot(res1,res2);abline(0,1)
	    # plot(res1,res3);abline(0,1)
	    # microbenchmark(
	    #   sample_MME_multiple_diagR(as.matrix(Y[,traits_j,drop=F]), W, Cholesky_C, pes[traits_j], chol_K_inv,tot_Eta_prec[traits_j], randn_theta[,traits_j,drop=F], randn_e[,traits_j,drop=F]),
	    #   sample_MME_multiple_diagR_c(as.matrix(Y[,traits_j,drop=F]), W, chol_C, pes[traits_j], chol_K_inv_2,tot_Eta_prec[traits_j], randn_theta[,traits_j,drop=F], randn_e[,traits_j,drop=F]),
	    #   sapply(traits_j,function(j) sample_MME_multiple_diagR_c2(as.matrix(Y[,j,drop=F]), W, chol_C, pes[j], chol_K_inv_2,tot_Eta_prec[j], randn_theta[,j,drop=F], randn_e[,j,drop=F]))
	    # )
	  },mc.cores = ncores)
	},
	sample_MME_ZKZts_c(Y,W,tot_Eta_prec,chol_Cs,chol_K_invs,h2s,h2s_index,randn_theta,randn_e,1)
	)
}
sample_MME_ZKZts = compiler::cmpfun(sample_MME_ZKZts)

sample_tot_prec = function(Y, tot_Eta_prec_shape, tot_Eta_prec_rate, Sigma_Choleskys,Sigma_Perm, h2s_index,ncores){
	n = nrow(Y)
	p = ncol(Y)

	Yp = Y
	if(!is.null(Sigma_Perm)){
		Yp = Sigma_Perm %*% Y
		Yp = matrix(Yp,nrow(Y))
	}
	unique_h2s = unique(h2s_index)
	unique_h2s_index = lapply(unique_h2s,function(x) which(h2s_index == x))
	scores = mclapply(seq_along(unique_h2s),function(j){
		h2_index = unique_h2s[j]
		Cholesky_Sigma = Sigma_Choleskys[[unique_h2s[j]]]$Cholesky_Sigma
		colSums(solve(Cholesky_Sigma,Yp[,unique_h2s_index[[j]]],'L')^2)
	},mc.cores = ncores)
	scores = unlist(scores)[order(unlist(unique_h2s_index))]
	matrix(rgamma(p,shape = tot_Eta_prec_shape + n/2, rate = tot_Eta_prec_rate + 1/2*scores),nrow=1)
}
sample_tot_prec = compiler::cmpfun(sample_tot_prec)

log_prob_h2 = function(i, y_p, Sigma_Choleskys, n, tot_Eta_prec, discrete_priors){
	Cholesky_Sigma = Sigma_Choleskys[[i]]$Cholesky_Sigma
	log_det_Sigma = Sigma_Choleskys[[i]]$log_det
	scores_2 = tot_Eta_prec*colSums(solve(Cholesky_Sigma,y_p,'L')^2)

	log_p = -n/2 * log(2*pi) - 1/2*(log_det_Sigma - n*log(tot_Eta_prec)) - 1/2 * scores_2 + log(discrete_priors)
	log_p
}
log_prob_h2 = compiler::cmpfun(log_prob_h2)

sample_h2s_discrete_MH = function(Y,tot_Eta_prec, Sigma_Choleskys,Sigma_Perm,discrete_priors,h2s_matrix,h2_index,step_size,ncores){
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

	# sample runif(p,0,1) before because parallel RNGs aren't consecutive.
	r_draws = runif(p)

	h2s_index = mclapply(1:p,function(j) {
			## steps:
			# 1) calculate log_prob at current state
			# 2) propose new state
			# 3) calc log_prob at new state
			# 4) calc proposal prob from old to new
			# 5) calc proposal prob from new to old
			old_state <- h2_index[j]
			old_log_p <- log_prob_h2(old_state,Yp[,j],Sigma_Choleskys,n,tot_Eta_prec[j],discrete_priors[old_state])

			candidate_new_states <- which(colSums((h2s_matrix - c(h2s_matrix[,old_state]))^2) < step_size)
			proposed_state <- sample(candidate_new_states,1)

			new_log_p <- log_prob_h2(proposed_state,Yp[,j],Sigma_Choleskys,n,tot_Eta_prec[j],discrete_priors[proposed_state])
			candidate_state_from_new_state <- which(colSums((h2s_matrix - c(h2s_matrix[,proposed_state]))^2) < step_size)

			forward_prob = 1/length(candidate_new_states)
			if(old_state %in% candidate_state_from_new_state){
				back_prob = 1/length(candidate_state_from_new_state)
			} else{
				back_prob = 0
			}

			log_MH_ratio = new_log_p - old_log_p + log(forward_prob) - log(back_prob)

			if(log(r_draws[j]) < log_MH_ratio) {
				return(proposed_state)
			} else {
				return(old_state)
			}
	},mc.cores = ncores)
	do.call(c,h2s_index)
}
sample_h2s_discrete_MH = compiler::cmpfun(sample_h2s_discrete_MH)

sample_h2s_discrete = function(Y,tot_Eta_prec, Sigma_Choleskys,Sigma_Perm,discrete_priors,ncores){
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
		scores_2 = tot_Eta_prec*colSums(solve(Cholesky_Sigma,Yp,'L')^2)

		log_ps = -n/2 * log(2*pi) - 1/2*(log_det_Sigma - n*log(tot_Eta_prec)) - 1/2 * scores_2 + log(discrete_priors[i])
		log_ps
	},mc.cores = ncores)
	if(length(log_ps) == 1) {
		log_ps = matrix(log_ps[[1]],ncol = p)
	} else{
		log_ps = do.call(rbind,log_ps)
	}
	h2s_index = sapply(1:p,function(j) {
		max_col = max(log_ps[,j])
		norm_factor = max_col+log(sum(exp(log_ps[,j]-max_col)))
		ps_j = exp(log_ps[,j] - norm_factor)
		sum(runif(1)>cumsum(ps_j))+1
	})
	return(h2s_index)
}
sample_h2s_discrete = compiler::cmpfun(sample_h2s_discrete)
