save_posterior_samples_ipx = function( sp_num, current_state, Posterior) {
	# Posterior,Lambda,F,F_a,B,W,E_a,delta,F_h2,resid_Y_prec,E_a_prec,W_prec) {
	# save posteriors. Full traces are kept of the more interesting parameters.
	# Only the posterior means are kept of less interesting parameters. These
	# should be correctly calculated over several re-starts of the sampler.

	Posterior = with(current_state, {
		# transform variables so that the variance of each column of F is 1.
		F_var = 1/F_a_prec + 1/F_e_prec
		F_h2 = F_e_prec / (F_e_prec + F_a_prec)
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
		Posterior
	})
	return(Posterior)
}
             