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
