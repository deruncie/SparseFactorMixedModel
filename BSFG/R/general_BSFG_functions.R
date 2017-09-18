#' #' Sample from Mixed Model equations given a diagonal K matrix
#' #'
#' #' Sample from Mixed Model equations given a diagonal K matrix with sparse but non-diagonal R
#' #'
#' #' This function draws a sample of the vector \eqn{theta} given the model:
#' #'     \deqn{y = W\theta + e}
#' #'     \deqn{\theta \sim N_b(\mu,K)} \deqn{e \sim N(0,R)} The algorithm follows that shown in the
#' #'     MCMCglmm course notes. This involves solving:
#' #'     \deqn{\tilde{\theta} = C^{-1}W'R^{-1}(1 - W\theta_{*} - e_{*})}
#' #'     \deqn{C = W'R^{-1}W + K^{-1}}
#' #'     \deqn{\theta_{*} \sim N(\mu,K), e_{*} \sim N(0,R)}
#' #'
#' #' To speed up the calculations, \eqn{R_hat = R*tot_Eta_prec = P' L L' P} is pre-calculated, and used
#' #' in place of \eqn{R^{-1}}.
#' #'
#' #' @param y vector of length n
#' #' @param W matrix of size n x b
#' #' @param Wp \eqn{R_perm * W} (see below)
#' #' @param prior_mean the vector \eqn{\mu}
#' #' @param prior_prec the diagonal of the matrix K^{-1}. All elements must be > 0
#' #' @param Cholesky_R Cholesky decomposition of the sparse matrix R calculated by
#' #'   \code{Cholesky(R_hat,perm=T,super=T)}. The factorization is  \eqn{P' L L' P}. R_hat here is
#' #'   \eqn{R*tot_Eta_prec}, i.e., the covariance up to normalization by the total variance.
#' #' @param chol_R the matrix \eqn{P' L} corresponding to the Cholesky_R decomposition above
#' #' @param R_Perm either \code{NULL} if \eqn{P} is diagonal, or the \eqn{P} matrix.
#' #' @param tot_Eta_prec the inverse of the total variance
#' #'
#' sample_MME_fixedEffects = function(Y,W,Sigma_Choleskys,h2s_index, tot_Eta_prec, prior_mean, prior_prec,ncores){
#' 	# using method described in MCMC Course notes
#' 	p = ncol(Y)
#' 	n = nrow(Y)
#' 	b = ncol(W)
#'
#' 	# pre-sample z-scores because draws from parallel processes not consecutive
#' 	randn_theta = matrix(rnorm(b*p),ncol = p)
#' 	randn_e = matrix(rnorm(n*p),ncol = p)
#'   res = sample_MME_fixedEffects_c(Y,W,Sigma_Choleskys,h2s_index,tot_Eta_prec,prior_mean,prior_prec,randn_theta,randn_e,1)
#'
#' 	return(res)
#' }
#' sample_MME_fixedEffects_cis = function(Y,W,cis_genotypes,cis_effects_index, Sigma_Choleskys,h2s_index, tot_Eta_prec, prior_mean, prior_prec,ncores){
#'   # using method described in MCMC Course notes
#'   p = ncol(Y)
#'   n = nrow(Y)
#'   b = ncol(W)
#'
#'   # pre-sample z-scores because draws from parallel processes not consecutive
#'   randn_theta = matrix(rnorm(b*p),ncol = p)
#'   randn_e = matrix(rnorm(n*p),ncol = p)
#'   randn_cis = rnorm(cis_effects_index[length(cis_effects_index)]-1)
#'   res = sample_MME_fixedEffects_cis_c(Y,as.matrix(W),cis_genotypes,Sigma_Choleskys,h2s_index,tot_Eta_prec,prior_mean,prior_prec,randn_theta,randn_e,randn_cis,cis_effects_index-1,1)
#'
#'   return(res)
#' }
#'
#' sample_MME_ZKZts = function(Y, W, tot_Eta_prec, randomEffect_C_Choleskys, h2s, h2s_index, ncores){
#' 	# using method described in MCMC Course notes
#' 	Y = as.matrix(Y)
#' 	p = ncol(Y)
#' 	n = nrow(Y)
#' 	b = ncol(W)
#'
#' 	# pre-sample z-scores because draws from parallel processes not consecutive
#' 	randn_theta = matrix(rnorm(b*p),ncol = p)
#' 	randn_e = matrix(rnorm(n*p),ncol = p)
#'
#' 	theta = sample_MME_ZKZts_c(Y,W,tot_Eta_prec,randomEffect_C_Choleskys,h2s,h2s_index,randn_theta,randn_e,1)
#' }
#'
#'
#' generate_candidate_states = function(h2s_matrix,step_size){
#'   candidate_states = lapply(1:dim(h2s_matrix)[2],function(i){
#'     h2_dist = colSums(abs(h2s_matrix - c(h2s_matrix[,i])))
#'     which(h2_dist < step_size & h2_dist > 0)
#'   })
#' }
#'
#' sample_h2s_discrete_MH = function(Y,tot_Eta_prec, Sigma_Choleskys,discrete_priors,h2s_matrix,h2_index,step_size,grainSize){
#' 	# while this works, it is much slower than doing the full scan over all traits, at least for multiple traits
#' 	# testing with p=100, solving the whole set takes ~4-5x solving just 1. And this method requires doing each trait separately
#' 	# both methods are similarly easy to multiplex, so no advantage there either.
#' 	n = nrow(Y)
#' 	p = ncol(Y)
#' 	discrete_bins = length(discrete_priors)
#'
#'
#' 	# sample runif(p,0,1) before because parallel RNGs aren't consecutive.
#' 	r_draws = runif(p)
#' 	state_draws = runif(p)
#'
#' 	h2s_index = sample_h2s_discrete_MH_c(Y,tot_Eta_prec,discrete_priors,h2_index,h2s_matrix,Sigma_Choleskys,r_draws,state_draws,step_size,grainSize)+1
#'
#' 		return(h2s_index)
#' }
#'
#' sample_h2s_discrete = function(Y,tot_Eta_prec, Sigma_Choleskys,discrete_priors,ncores){
#' 	n = nrow(Y)
#' 	p = ncol(Y)
#'
#' 	log_ps = log_p_h2s(Y,tot_Eta_prec,Sigma_Choleskys,discrete_priors,1)
#' 	rs = runif(p)
#' 	h2s_index = sample_h2s(log_ps,rs,1)
#' 	return(h2s_index)
#' }
