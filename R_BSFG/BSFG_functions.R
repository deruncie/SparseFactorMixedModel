#library(R.matlab)

require(pracma)


#cholcov = function(X){
	# calculates a matrix U such that t(U) %*% U == X for X that is not PD
#	E = svd(X)
#	cols = E$d > 1e-14
#	U = E$u[,cols] %*% diag(sqrt(E$d[cols]))
#	return(t(U))
#}


cholcov = function(X){
  # calculates a matrix U such that t(U) %*% U == X for X that is not PD
  E = svd(X)
 # cols = E$d > 1e-14
  U = E$u %*% diag(sqrt(E$d))
  return(t(U))
}

sample_Lambda = function( Y_tilde,F,resid_Y_prec, E_a_prec,Plam,invert_aI_bZAZ ) {

	#Sample factor loadings Lambda while marginalizing over residual
	#genetic effects: Y - Z_2W = F*Lambda' + E, vec(E)~N(0,kron(Psi_E,In) + kron(Psi_U, ZAZ^T))
	#note: conditioning on F, but marginalizing over E_a.
	#sampling is done separately by trait because each column of Lambda is
	#independent in the conditional posterior
	#note: invert_aI_bZAZ has parameters that diagonalize aI + bZAZ for fast
	#inversion: inv(aI + bZAZ) = 1/b*U*diag(1./(s+a/b))*U'

	p = length(resid_Y_prec)
	k = ncol(F)

	U = invert_aI_bZAZ$U
	s = invert_aI_bZAZ$s

	FtU = t(F) %*% U
	UtY = t(U) %*% Y_tilde

	Zlams = matrix(rnorm(k*p),nr = k, nc = p)
	Lambda = matrix(0,nr = p,nc = k)

	for(j in 1:p) {
	  means = c()
	  Qlam = c()
	  for(pop in pops){
		  FUDi = E_a_prec[j] * sweep(FtU,2,(s + E_a_prec[j]/resid_Y_prec[j]),'/')
		  means = means + FUDi %*% UtY[,j]
		  Qlam = Qlam + FUDi %*% t(FtU) 
		}
	  Qlam = Qlam + diag(Plam[j,])

		# recover()
		Llam = t(chol(Qlam))
		vlam = forwardsolve(Llam,means)
		mlam = backsolve(t(Llam),vlam)
		ylam = backsolve(t(Llam),Zlams[,j])

		Lambda[j,] = ylam + mlam

	}

	return(Lambda)
}


sample_prec_discrete_conditional = function(Y,h2_divisions,h2_priors,invert_aI_bZAZ,res_prec) {
	#sample factor heritibilties conditional on a given residual precision
	#(res_precision)
	#prior given as relative weights on each of h2_divisions points. Doesn't
	#have to be normalized
	#samples conditional on F, marginalizes over F_a.
	#uses invert_aI_bZAZ.U and invert_aI_bZAZ.s to not have to invert aI + bZAZ
	#each iteration.
	#ident_prec is a in the above equation.

	# matlab = readMat('sample_prec_discrete_conditional_data.mat')
	# Y = matlab$Y
	# h2_divisions = c(matlab$h2.divisions)
	# h2_priors = c(matlab$h2.priors)
	# invert_aI_bZAZ = matlab$invert.aI.bZAZ
	# res_prec = matlab$res.prec
	# names(invert_aI_bZAZ) = c('U','s')

	# other code doesn't allow Prec == Inf. So should set prior h2_priors[1] = 0

	U = invert_aI_bZAZ$U
	s = invert_aI_bZAZ$s

	p = ncol(Y)
	n = nrow(Y)
	Trait_h2 = rep(0,p)

	log_ps = matrix(NA,p,h2_divisions)
	std_scores_b = t(Y) %*% U
	for(i in 1:h2_divisions) {
		h2 = (i-1)/(h2_divisions)
		if(h2 > 0){
			std_scores = sweep(sweep(std_scores_b,2,sqrt(s+(1-h2)/h2),'/'),1,sqrt(h2/(res_prec*(1-h2))),'/')
			det = colSums(log((s+(1-h2)/h2) %*% t(h2/(res_prec*(1-h2))))/2)
		} else {
			std_scores = sweep(t(Y),1,sqrt(1/res_prec),'/')
			det = n/2*log(1/res_prec)
		}
		log_ps[,i] = rowSums(dnorm(std_scores,log=T)) - det + log(h2_priors[i])
	}
	for(j in 1:p){
		norm_factor = max(log_ps[j,])+log(sum(exp(log_ps[j,]-max(log_ps[j,]))))
		ps_j = exp(log_ps[j,] - norm_factor)
		log_ps[j,] = ps_j
		Trait_h2[j] = sum(runif(1)>cumsum(ps_j))/(h2_divisions)
	}
	Prec = (res_prec*(1-Trait_h2))/Trait_h2

	return(Prec)
}


sample_h2s_discrete = function(F,h2_divisions,h2_priors,invert_aI_bZAZ){
	#sample factor heritibilties from a discrete set on [0,1)
	#prior places 50% of the weight at h2=0
	#samples conditional on F, marginalizes over F_a.
	#uses invert_aI_bZAZ.U and invert_aI_bZAZ.s to not have to invert aI + bZAZ
	#each iteration.

	# matlab = readMat('sample_h2s_discrete_data.mat')
	# for(i in 1:10) names(matlab) = sub('.','_',names(matlab),fixed=T)
	# F = matlab$F
	# h2_divisions = c(matlab$h2_divisions)
	# h2_priors = matlab$h2_priors
	# invert_aI_bZAZ = matlab$invert_aI_bZAZ
	# names(invert_aI_bZAZ) = c('U','s')

	U = invert_aI_bZAZ$U;
	s = invert_aI_bZAZ$s;

	k = ncol(F)
	F_h2 = rep(0,k)

	log_ps = matrix(0,k,h2_divisions)
	std_scores_b = t(F) %*% U

	for(i in 1:h2_divisions){
		h2 = (i-1)/(h2_divisions)
		if(h2 > 0) {
			std_scores = 1/sqrt(h2) * sweep(std_scores_b,2,sqrt(s+(1-h2)/h2),'/')
			det = sum(log((s+(1-h2)/h2)*h2)/2)
		} else {
			# std_scores = std_scores_b # note. This is the same with F or std_scores_b
			std_scores = t(F)
			det = 0
		}
		log_ps[,i] = rowSums(dnorm(std_scores,log=T)) - det + log(h2_priors[i])
	}
	for(j in 1:k) {
		norm_factor = max(log_ps[j,])+log(sum(exp(log_ps[j,]-max(log_ps[j,]))))
		ps_j = exp(log_ps[j,] - norm_factor)
		log_ps[j,] = ps_j
		F_h2[j] = sum(runif(1)>cumsum(ps_j))/(h2_divisions);
	}

	return(F_h2)
}


sample_means = function( Y_tilde, resid_Y_prec, E_a_prec, invert_aPXA_bDesignDesignT ) {
	#when used to sample [B;E_a]:
	# W - F*Lambda' = X*B + Z_1*E_a + E, vec(E)~N(0,kron(Psi_E,In)). 
	# Note: conditioning on F, Lambda and W.
	#The vector [b_j;E_{a_j}] is sampled simultaneously. Each trait is sampled separately because their
	#conditional posteriors factor into independent MVNs.
	#note:invert_aPXA_bDesignDesignT has parameters to diagonalize mixed model equations for fast inversion: 
	#inv(a*blkdiag(fixed_effects_prec*eye(b),Ainv) + b*[X Z_1]'[X Z_1]) = U*diag(1./(a.*s1+b.*s2))*U'
	#Design_U = [X Z_1]*U, which doesn't change each iteration. 
	
	# matlab = readMat('sample_means_data.mat')
	# for(i in 1:10) names(matlab) = sub('.','_',names(matlab),fixed=T)
	# Y_tilde = matlab$Y_tilde
	# resid_Y_prec = matlab$resid_Y_prec
	# E_a_prec = matlab$E_a_prec
	# invert_aPXA_bDesignDesignT = matlab$invert_aPXA_bDesignDesignT
	# names(invert_aPXA_bDesignDesignT) = c('U','s1','s2','Design_U')

	U = invert_aPXA_bDesignDesignT$U
	s1 = invert_aPXA_bDesignDesignT$s1
	s2 = invert_aPXA_bDesignDesignT$s2
	Design_U = invert_aPXA_bDesignDesignT$Design_U

	n = nrow(Y_tilde)
	p = length(resid_Y_prec)
	br = ncol(Design_U)

	means = sweep(t(Design_U) %*% Y_tilde,2,resid_Y_prec,'*')
	location_sample = matrix(0,nr = br,nc = p)
	Zlams = matrix(rnorm(br*p),nr = br, nc = p)

	# Zlams = matlab$Zlams
	for(j in 1:p) {
		d = s1*E_a_prec[j] + s2*resid_Y_prec[j]
		mlam = means[,j] / d
		location_sample[,j] = U %*% (mlam + Zlams[,j]/sqrt(d))
	}
	# location_sample = t(location_sample);

	return(location_sample)
}


sample_F_a = function(F,Z_1,F_h2,invert_aZZt_Ainv) {
	#samples genetic effects on factors (F_a) conditional on the factor scores F:
	# F_i = F_{a_i} + E_i, E_i~N(0,s2*(1-h2)*I) for each latent trait i
	# U_i = zeros(r,1) if h2_i = 0
	# it is assumed that s2 = 1 because this scaling factor is absorbed in
	# Lambda
	# invert_aZZt_Ainv has parameters to diagonalize a*Z_1*Z_1' + b*I for fast
	# inversion:
	
	U = invert_aZZt_Ainv$U
	s1 = invert_aZZt_Ainv$s1
	s2 = invert_aZZt_Ainv$s2

	k = ncol(F)
	r = ncol(Z_1)
	tau_e = 1/(1-F_h2)
	tau_u = 1/F_h2

	b = t(U) %*% t(Z_1) %*% sweep(F,2,tau_e,'*')
	z = matrix(rnorm(r*k),nr = r, nc = k)
	F_a = matrix(0,nr = r,nc = k)

	for(j in 1:k) {
		if(tau_e[j]==1) {
			F_a[,j] = 0
		} else if(tau_e[j] == Inf) {
			F_a[,j] = F[,j]
		} else {
			d = s2*tau_u[j] + s1*tau_e[j]
			mlam = b[,j] / d
			F_a[,j] = U %*% (mlam + z[,j]/sqrt(d))
		}
		if(tau_e[j] == 1) {
		   #  disp(F_a[,j])
		}
	}

	return(F_a)
}



sample_factors_scores = function( Y_tilde, Z_1,Lambda,resid_Y_prec,F_a,F_h2 ) {
#Sample factor scores given factor loadings (F_a), factor heritabilities (F_h2) and
#phenotype residuals

	# matlab = readMat('sample_Factor_scores_data.mat')
	# for(i in 1:10) names(matlab) = sub('.','_',names(matlab),fixed=T)
	# Y_tilde = matlab$Y_tilde
	# X = matlab$X
	# Z_1 = matlab$Z_1
	# Lambda = matlab$Lambda
	# resid_Y_prec = matlab$resid_Y_prec
	# F_b = matlab$F_b
	# F_a = matlab$F_a
	# F_h2 = c(matlab$F_h2)

	Lmsg = sweep(Lambda,1,resid_Y_prec,'*')
	tau_e = 1/(1-c(F_h2))
	S = chol(t(Lambda) %*% Lmsg + diag(tau_e))
	tS = t(S)

	# Meta = t(backsolve(S,t(Y_tilde %*% Lmsg + sweep(Z_1 %*% F_a,2,tau_e,'*'))))
	# F = t(forwardsolve(t(S),t(Meta + matrix(rnorm(prod(dim(Meta))),nr = nrow(Meta)))))
	Meta = t(forwardsolve(tS,t(Y_tilde %*% Lmsg + sweep(Z_1 %*% F_a,2,tau_e,'*'))))
	F = t(backsolve(S,t(Meta + matrix(rnorm(prod(dim(Meta))),nr = nrow(Meta)))))
	return(F)
}



sample_delta = function( delta,tauh,Lambda_prec,delta_1_shape,delta_1_rate,delta_2_shape,delta_2_rate,Lambda2 ) {
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
	
	return(delta)
}


update_k = function( F,Lambda,F_a,F_h2,Lambda_prec,Plam,delta,tauh,Z_W,Lambda_df,delta_2_shape,delta_2_rate,b0,b1,i,epsilon,prop ) {
# adapt the number of factors by dropping factors with only small loadings
# if they exist, or adding new factors sampled from the prior if all factors
# appear important. The rate of adaptation decreases through the chain,
# controlled by b0 and b1. Should work correctly over continuations of
# previously stopped chains.

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
			F_a = cbind(F_a,rnorm(r,0,sqrt(F_h2[k])))
			F = cbind(F,rnorm(n,Z_W %*% F_a[,k],sqrt(1-F_h2[k])))
		
			} else if(num > 0) { # drop redundant columns
			nonred = which(vec == 0) # non-redundant loadings columns
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
			F_a = F_a[,nonred]
		}
	}


	return(list(
		F           = F,
		Lambda      = Lambda,
		F_a         = F_a,
		F_h2        = F_h2,
		Lambda_prec = Lambda_prec,
		Plam        = Plam,
		delta       = delta,
		tauh         = tauh
		))
}


save_posterior_samples = function( sp_num,Posterior,Lambda,F,F_a,B,W,E_a,delta,F_h2,resid_Y_prec,E_a_prec,W_prec) {
	# save posteriors. Full traces are kept of the more interesting parameters.
	# Only the posterior means are kept of less interesting parameters. These
	# should be correctly calculated over several re-starts of the sampler.

	sp = ncol(Posterior$Lambda)
	#size(Posterior$Lambda,2)
  ncol(Posterior$Lambda)
	#save factor samples
	if(length(Lambda) > nrow(Posterior$Lambda)){
		# expand factor sample matrix if necessary
		Posterior$Lambda = rbind(Posterior$Lambda, 	matrix(0,nr = length(Lambda)-nrow(Posterior$Lambda),nc = sp))
		Posterior$F      = rbind(Posterior$F, 	   	matrix(0,nr = length(F)     -nrow(Posterior$F),		nc = sp))
		Posterior$F_a    = rbind(Posterior$F_a, 	matrix(0,nr = length(F_a) 	-nrow(Posterior$F_a),	nc = sp))
		Posterior$delta  = rbind(Posterior$delta, 	matrix(0,nr = length(delta) -nrow(Posterior$delta),	nc = sp))
		Posterior$F_h2   = rbind(Posterior$F_h2, 	matrix(0,nr = length(F_h2) 	-nrow(Posterior$F_h2),	nc = sp))
	}
	Posterior$Lambda[1:length(Lambda),sp_num] = c(Lambda)
	Posterior$F[1:length(F),sp_num]     = c(F)
	Posterior$F_a[1:length(F_a),sp_num] = c(F_a)
	Posterior$delta[1:length(delta),sp_num] = delta
	Posterior$F_h2[1:length(F_h2),sp_num] = F_h2

	Posterior$resid_Y_prec[,sp_num] = resid_Y_prec
	Posterior$E_a_prec[,sp_num]     = E_a_prec
	Posterior$W_prec[,sp_num]       = W_prec

	# save B,U,W
	Posterior$B   = (Posterior$B*(sp_num-1) + B)/sp_num
	Posterior$E_a = (Posterior$E_a*(sp_num-1) + E_a)/sp_num
	Posterior$W   = (Posterior$W*(sp_num-1) + W)/sp_num
	return(Posterior)
}
  
save_posterior_samples_comb = function( sp_num,Posterior,Lambda,F,F_a,B,W,E_a,delta,F_h2,resid_Y_prec,E_a_prec,W_prec) {
  # save posteriors. Full traces are kept of the more interesting parameters.
  # Only the posterior means are kept of less interesting parameters. These
  # should be correctly calculated over several re-starts of the sampler.
  
  # first unlist F, F_a, F_h2, Y, B, E_a
  nl = length(F_h2)
  F = unlist(do.call("rbind", F)) # n*k
  F_a = unlist(do.call("rbind", F_a)) #n*k
  F_h2 = Reduce("+", F_h2)/nl # k taking average
  B = Reduce("+",B)/nl # b*p taking average
  E_a = unlist(do.call("rbind", E_a)) # n*k
  W = Reduce("+",W)/nl  # k*p taking average
  resid_Y_prec = Reduce("+",resid_Y_prec)/nl #p taking average
  E_a_prec = Reduce("+",E_a_prec)/nl #p
  W_prec = Reduce("+",W_prec)/nl  #p
  
  #sp = ncol(Posterior$Lambda)
  sp = ncol(Posterior$Lambda)
  #size(Posterior$Lambda,2)
  #save factor samples
  # #if(length(Lambda) > nrow(Posterior[["Lambda"]])){  this will not happen since we stop updatek
  #   # expand factor sample matrix if necessary
  #   Posterior[["F"]]      = rbind(Posterior[["F"]], 	   	matrix(0,nr = length(F)     -nrow(Posterior[["F"]]),		nc = sp))
  #   Posterior[["F_a"]]    = rbind(Posterior[["F_a"]], 	matrix(0,nr = length(F_a) 	-nrow(Posterior[["F_a"]]),	nc = sp))
  #   Posterior[["delta"]]  = rbind(Posterior[["delta"]], 	matrix(0,nr = length(delta) -nrow(Posterior[["delta"]]),	nc = sp))
  #   Posterior[["F_h2"]]   = rbind(Posterior[["F_h2"]], 	matrix(0,nr = length(F_h2) 	-nrow(Posterior[["F_h2"]]),	nc = sp))
  # }
  Posterior$Lambda[1:length(Lambda),sp_num] = c(Lambda)
  Posterior$F[1:length(F),sp_num]     = c(F)
  Posterior$F_a[1:length(F_a),sp_num] = c(F_a)
  Posterior$delta[1:length(delta),sp_num] = delta
  Posterior$F_h2[1:length(F_h2),sp_num] = F_h2
  
  Posterior$resid_Y_prec[,sp_num] = resid_Y_prec
  Posterior$E_a_prec[,sp_num]     = E_a_prec
  Posterior$W_prec[,sp_num]       = W_prec
  
  # save B,U,W
  Posterior$B   = (Posterior$B*(sp_num-1) + B)/sp_num
  Posterior$E_a = (Posterior$E_a*(sp_num-1) + E_a)/sp_num
  Posterior$W   = (Posterior$W*(sp_num-1) + W)/sp_num
  
  return(Posterior)
}


save_posterior_samples_fixedlambda = function( j,Posterior,F,F_a,B,W,E_a,F_h2,resid_Y_prec,E_a_prec,W_prec) {
  # save posteriors. Full traces are kept of the more interesting parameters.
  # Only the posterior means are kept of less interesting parameters. These
  # should be correctly calculated over several re-starts of the sampler.
  
  sp = ncol(Posterior$F)
  #size(Posterior$Lambda,2)
  ncol(Posterior$F)
  #save factor samples
  if(length(F) > nrow(Posterior$F)){
    # expand factor sample matrix if necessary
    Posterior$F      = rbind(Posterior$F, 	   	matrix(0,nr = length(F)     -nrow(Posterior$F),		nc = sp))
    Posterior$F_a    = rbind(Posterior$F_a, 	matrix(0,nr = length(F_a) 	-nrow(Posterior$F_a),	nc = sp))
    Posterior$F_h2   = rbind(Posterior$F_h2, 	matrix(0,nr = length(F_h2) 	-nrow(Posterior$F_h2),	nc = sp))
  }
  Posterior$F[1:length(F),j]     = c(F)
  Posterior$F_a[1:length(F_a),j] = c(F_a)
  Posterior$F_h2[1:length(F_h2),j] = F_h2
  
  Posterior$resid_Y_prec[,j] = resid_Y_prec
  Posterior$E_a_prec[,j]     = E_a_prec
  Posterior$W_prec[,j]       = W_prec
  
  # save U,W

  Posterior$E_a = (Posterior$E_a*(j-1) + E_a)/j
  Posterior$W   = (Posterior$W*(j-1) + W)/j
  
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
    r2 = nrow(Posterior$W)
    
    Posterior = list(
		    Lambda        = matrix(0,nr=0,nc=0),
		    F_a           = matrix(0,nr=0,nc=0),
		    F             = matrix(0,nr=0,nc=0),
		    delta         = matrix(0,nr=0,nc=0),
		    F_h2          = matrix(0,nr=0,nc=0),
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

G_Matrix_Comp = function(BSFG_state){
  Posterior      = BSFG_state$Posterior
  traitnames     = BSFG_state$traitnames
  run_variables  = BSFG_state$run_variables
  sp_num = ncol(Posterior$Lambda) 
  n = run_variables$n
  p = run_variables$p
  k1 = nrow(Posterior$F_a)/n;
  k2 = nrow(Posterior$Lambda)/p;
  if (k2 >= k1){
    k = k1
  }else{
    k=k2
  }
  h2s = Posterior$F_h2[,1:sp_num]
  #G_Lambdas = array(0,dim = dim(Posterior$Lambda))
  #Lambda_est = matrix(0,p,k)
  G_est = E_est = matrix(0,p,p)
  #traces_G = matrix(,p*(p+1)/2,sp_num)
  #traces_G_cor = matrix(,p*(p+1)/2,sp_num)
  #traces_E = matrix(,p*(p+1)/2,sp_num)
  GMatrix = NULL
  for(j in 1:sp_num) {
    Lj = matrix(Posterior$Lambda[,j],p,k)
    h2j = Posterior$F_h2[,j]
    G_Lj = Lj %*%  diag(sqrt(h2j))
    #G_Lambdas[,j] = c(G_Lj)
    Gj = G_Lj %*%  t(G_Lj) + diag(1/Posterior$E_a_prec[,j])
    rownames(Gj) = traitnames
    #posterior mean
    GMatrix[[j]] = Gj
    G_est = G_est + Gj/sp_num
    #library(gdata)
    #traces_G[,j] = lowerTriangle(Gj,diag = TRUE)
    #traces_G_cor[,j] = lowerTriangle(CovToCor(Gj),diag = TRUE)
    
    #E_Lj = Lj  %*% diag(1-h2j) %*%  t(Lj) + diag(1/Posterior$resid_Y_prec[,j])
    #E_est = E_est + E_Lj/sp_num;
    #posterior mean for lambda
    #Lambda_est = Lambda_est + matrix(Posterior$Lambda[,j],p,k)/sp_num;
    #traces_E[,j] = lowerTriangle(E_Lj,diag = TRUE)
  }
  #G_Lambda = matrix(rowMeans(G_Lambdas),p,k)
  return(GMatrix)
}

G_Traces_Comp = function(BSFG_state){
  Posterior      = BSFG_state$Posterior
  traitnames     = BSFG_state$traitnames
  run_variables  = BSFG_state$run_variables
  n = run_variables$n
  sp_num = ncol(Posterior$F_a) 
  p = run_variables$p
  # the number of latent traits may be different for different model. If we fixed lambda
  # need to cut the redudant latent traits or add some latent traits
  k1 = nrow(Posterior$F_a)/n;
  k2 = nrow(Posterior$Lambda)/p;
  if (k2 >= k1){
    k = k1
  }else{
    k=k2
  }
  h2s = Posterior$F_h2[,1:sp_num]
  #G_Lambdas = array(0,dim = dim(Posterior$Lambda))
  Lambda_est = matrix(0,p,k)
  G_est = E_est = matrix(0,p,p)
  traces_G = matrix(,p*(p+1)/2,sp_num)
  #traces_G_cor = matrix(,p*(p+1)/2,sp_num)
  #traces_E = matrix(,p*(p+1)/2,sp_num)
  GMatrix = NULL
  for(j in 1:sp_num) {
    Lj = matrix(Posterior$Lambda[p*k,j],p,k)
    h2j = Posterior$F_h2[,j]
    G_Lj = Lj %*%  diag(sqrt(h2j))
    #G_Lambdas[,j] = c(G_Lj)
    Gj = G_Lj %*%  t(G_Lj) + diag(1/Posterior$E_a_prec[,j])
    rownames(Gj) = traitnames
    #posterior mean
    GMatrix[[j]] = Gj
    #G_est = G_est + Gj/sp_num
    library(gdata)
    traces_G[,j] = lowerTriangle(Gj,diag = TRUE)
    #traces_G_cor[,j] = lowerTriangle(CovToCor(Gj),diag = TRUE)
    
    #E_Lj = Lj  %*% diag(1-h2j) %*%  t(Lj) + diag(1/Posterior$resid_Y_prec[,j])
    #E_est = E_est + E_Lj/sp_num;
    #posterior mean for lambda
    #Lambda_est = Lambda_est + matrix(Posterior$Lambda[,j],p,k)/sp_num;
    #traces_E[,j] = lowerTriangle(E_Lj,diag = TRUE)
  }
  #G_Lambda = matrix(rowMeans(G_Lambdas),p,k)
  return(traces_G)
}