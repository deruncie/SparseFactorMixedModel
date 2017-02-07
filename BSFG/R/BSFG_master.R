#' Initialize a BSFG model
#'
#' Sets up the BSFG model, selects starting values, and pre-calculates matrices for the GIBBS
#' sampler.
#'
#' @param Y a n x p matrix of data (n individuals x p traits)
#' @param fixed  formula for the fixed effects. \strong Note: currently only applies to residuals
#'   (not factors)
#' @param random formula for the random effects. \strong Note: currently should be of the form \code{~ x1
#'   + x2} (no interactions)
#' @param data data.frame with n rows containing columns corresponding to the fixed and random
#'   effects
#' @param priors list providing hyperparameters for the model priors. See details
#' @param run_parameters list providing various parameters for the model run
#' @param A_mats list of covariance matrices for random effects. If none provided (and none provided
#'   for A_inv_mats), assumed to be the identity.
#' @param A_inv_mats list of covariance matrices for random effects. If none provided (and none
#'   provided for A_mats), assumed to be the identity.
#' @param fixed_Factors currently no effect. To be used to specify fixed effect model for factors
#' @param scaleY should columns of Y be re-scaled to mean = 0, var = 1?
#' @param sampler ('fast_BSFG','general_BSFG'). Specifies which Gibbs sampler to use. If more than
#'   one random effect, defaults to general_BSFG, regardless of this parameter
#' @param ncores for \code{general_BSFG}, number of cores to use. \code{fast_BSFG} is parallelize by RcppParallel
#' @param simulation is this a simulation (and are known values included in setup?)?
#' @param setup a list of known values for Lambda (error_factor_lambda), h2, factor_h2s
#' @param verbose (T,F) should progress in initialization be reported?
#'
#' @return An object of class BSFG_state with components:
#' @return current_state: a list of parameters in the current iteration of the sampler
#' @return Posterior: a list of arrays of posterior samples
#' @return RNG current state of R's Random number generator (for re-starting chaings)
#' @return traitnames: vector of trait names (from colnames of Y)
#' @return run_parameters, run_variables, data_matrices, priors, simulation
#'
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
	class(BSFG_state) = append(class(BSFG_state),'BSFG_state')
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

summary.BSFG_state = function(BSFG_state){
  with(BSFG_state,{
    cat(
      c(sprintf('\n BSFG_state object for data of size %d x %d \n',nrow(data_matrices$Y),ncol(data_matrices$Y))),
      c(sprintf('Model dimensions: fixed = %d, random = %d \n',ncol(data_matrices$X),ncol(data_matrices$Z))),
      c(sprintf('Sampler: %s \n',run_parameters$sampler)),
      c(sprintf('Current iteration: %d, Posterior_samples: %d \n',current_state$nrun,Posterior$sp_num)),
      c(sprintf('Current factor dimension: %d factors \n',ncol(current_state$Lambda))),
      c(sprintf('Total time: %s \n\n',format(current_state$total_time)))
    )
  })
}

print.BSFG_state = function(BSFG_state){
  with(BSFG_state,{
    cat(
      c(sprintf('\n Current iteration: %d, Posterior_samples: %d \n',current_state$nrun,Posterior$sp_num)),
      c(sprintf('Current factor dimension: %d factors \n',ncol(current_state$Lambda))),
      c(sprintf('Total time: %s \n\n',format(current_state$total_time)))
    )
  })
}

plot.BSFG_state = function(BSFG_state){
  if(BSFG_state$simulation){
    draw_simulation_diagnostics(BSFG_state)
  } else{
    draw_results_diagnostics(BSFG_state)
  }
}
