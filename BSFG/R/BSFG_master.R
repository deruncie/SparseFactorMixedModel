# chooses which sampler to run, stores this in BSFG_state, calls sampler function

BSFG_init = function(Y, fixed, random, data, priors, run_parameters, A_mats = NULL, A_inv_mats = NULL,
                                    fixed_Factors = NULL, scaleY = TRUE, sampler = 'fast_BSFG',
                                    ncores = 1,simulation = F,setup = NULL,verbose=T)
{
	RE_names = rownames(attr(terms(random),'factors'))

	if(length(RE_names) > 1){
		print(sprintf('%d random effects. Using "general_BSFG" sampler',length(RE_names)))
		sampler = 'general_BSFG'
	}

	if(sampler == 'fast_BSFG'){
		BSFG_state = fast_BSFG_init(Y, fixed, random, data, priors, run_parameters, A_mats, A_inv_mats,
												fixed_Factors, scaleY, simulation, setup, verbose)
	} else{
		BSFG_state = general_BSFG_init(Y, fixed, random, data, priors, run_parameters, A_mats, A_inv_mats,
												fixed_Factors, scaleY, ncores, simulation, setup, verbose)
	}
	BSFG_state$run_parameters$sampler = sampler
	return(BSFG_state)
}

BSFG_sampler = function(BSFG_state,n_samples, ncores){
	sampler = BSFG_state$run_parameters$sampler
	if(sampler == 'fast_BSFG'){
		if(BSFG_state$run_parameters$verbose) print('fast_BSFG sampler')
		BSFG_state = fast_BSFG_sampler(BSFG_state,n_samples)
	} else{
		if(BSFG_state$run_parameters$verbose) print('general_BSFG sampler')
		BSFG_state = general_BSFG_sampler(BSFG_state,n_samples, ncores)
	}
	BSFG_state$run_parameters$sampler = sampler
	return(BSFG_state)
}
