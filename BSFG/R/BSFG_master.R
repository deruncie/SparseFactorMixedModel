#' Set BSFG run parameters
#'
#' Function to create run_parameters list for initializing BSFG model
#'
#' @param sampler specify the sampler to use. fast_BSFG is much faster, but only allows one random
#'   effect. If more are specified in \code{BSFG_init}, this is switched to general_BSFG.
#' @param fixed_factors should the fixed effect model be applied to factors as well as the factor
#'   residuals? If so, the first column of the design matrix is dropped.
#' @param simulaiton Is this a fit to simulated data? If so, a setup list will be expected providing
#'   the true values
#' @param scale_Y Should the Y values be centered and scaled? Recommend, except for simulated data.
#' @param b0 parameter of the \code{update_k} function. See Bhattacharya and Dunson 2011
#' @param b1 parameter of the \code{update_k} function. See Bhattacharya and Dunson 2011
#' @param epsilon parameter of the \code{update_k} function. Smallest \eqn{\lambda_{ij}} that is
#'   considered "large", signifying a factor should be kept
#' @param prop proportion of \eqn{\lambda{ij}} elements in a column of \eqn{\Lambda} that must be smaller than
#'   \code{epsilon} before factor is dropped
#' @param kinit initial number of factors
#' @param h2_divisions A scalar or vector of length equal to number of random effects, specifying
#'   the number of equally spaced divisions to devide each random effect variance component. All
#'   variance componets are scaled as percent_of_total. This number of divisions are created for
#'   each component, and then combined with \code{expand.grid()}. Then combinations of variances that sum to less than 1 are selected as valid.
#' @param burn burnin length of the MCMC chain
#' @param thin thinning rate of the MCMC chain
#' @seealso \code{\link{BSFG_init}}, \code{\link{sample_BSFG}}, \code{\link{print.BSFG_state}}
#'
BSFG_control = function(sampler = c('fast_BSFG','general_BSFG'),fixed_factors = c(T,F),simulation = c(F,T),scale_Y = c(T,F),
                        b0 = 1, b1 = 0.0005, epsilon = 1e-1, prop = 1.00,
                        k_init = 20, h2_divisions = 100,
                        burn = 100,
                        thin = 2) {

  all_args = formals()
  passed_args = as.list(match.call())[-1]
  all_args[names(passed_args)] = passed_args
  return(passed_args)
}

#' Initialize a BSFG model
#'
#' Sets up the BSFG model, selects starting values, and pre-calculates matrices for the GIBBS
#' sampler.
#'
#' The first step in fitting a BSFG model. This function setups up the model matrices based on the
#' fixed and random effect formulas provided, initializes all of model parameters, and then
#' pre-calculates a set of matrices and transformations to speed up each iteration of the Gibbs
#' sampler.
#'
#'
#' @param Y a n x p matrix of data (n individuals x p traits)
#' @param fixed  formula for the fixed effects
#' @param random formula for the random effects. Multiple random effects can be specified as
#'   ~random1+random2+..., but currently there is no support for interactions. Each random effect
#'   must correspond to a column of data (which should be a factor)
#' @param data data.frame with n rows containing columns corresponding to the fixed and random
#'   effects
#' @param priors list providing hyperparameters for the model priors. T
#' @param run_parameters list providing various parameters for the model run. See \link{BSFG_control}
#' @param A_mats list of covariance matrices for random effects. If none provided (and none provided
#'   for A_inv_mats), assumed to be the identity.
#' @param A_inv_mats list of covariance matrices for random effects. If none provided (and none
#'   provided for A_mats), assumed to be the identity.
#' @param ncores for \code{general_BSFG}, number of cores to use during initialization.
#' @param setup optional - a list of known values for Lambda (error_factor_lambda), h2, factor_h2s
#' @param verbose should progress in initialization be reported?
#'
#' @return An object of class BSFG_state with components:
#' @return current_state: a list of parameters in the current iteration of the sampler
#' @return Posterior: a list of arrays of posterior samples
#' @return RNG current state of R's Random number generator (for re-starting chaings)
#' @return traitnames: vector of trait names (from colnames of Y)
#' @return run_parameters, run_variables, data_matrices, priors, simulation: input data and
#'   parameters
#' @seealso \code{\link{BSFG_control}}, \code{\link{sample_BSFG}}, \code{\link{print.BSFG_state}}, \code{\link{plot.BSFG_state}}
#'
BSFG_init = function(Y, fixed, random, data, priors, run_parameters, A_mats = NULL, A_inv_mats = NULL,
                     sampler = c('fast_BSFG','general_BSFG'), ncores = 1,simulation = c(F,T),setup = NULL,verbose=T) {


  # ----------------------------- #
  # ---- build model matrices --- #
  # ----------------------------- #

	# model dimensions
	n = nrow(Y)
	p = ncol(Y)
	traitnames = colnames(Y)

	# missing data
	Y_missing = Matrix(is.na(Y))

	# scale Y
	if(run_parameters$scale_Y){
	  Mean_Y = colMeans(Y,na.rm=T)
	  VY = apply(Y,2,var,na.rm=T)
	  Y = sweep(Y,2,Mean_Y,'-')
	  Y = sweep(Y,2,sqrt(VY),'/')
	} else {
	  Mean_Y = rep(0,p)
	  VY = rep(1,p)
	}

	# build X from fixed model
	  # for residuals
	X = model.matrix(fixed,data)
	b = ncol(X)

	  # for factors.
	  #   If fixed_factors, then fixed effects are modeled for factors
	  #     in this case, the prior variance should be modeled
	  #   If fixed_factors == F, then an empty matrix
	if(run_parameters$fixed_factors){
	  X_F = X[,-1,drop=FALSE]
	} else{
	  X_F = matrix(0,nr = n, ncol = 0)
	}
	b_F = ncol(X_F)

	# build Z matrices from random model
	RE_names = rownames(attr(terms(random),'factors'))
	n_RE = length(RE_names)

	if(n_RE > 1 && run_parameters$sampler == 'fast_BSFG'){
	  print(sprintf('%d random effects. Using "general_BSFG" sampler',length(RE_names)))
	  run_parameters$sampler = 'general_BSFG'
	}

	Z_matrices = lapply(RE_names,function(re) {
	  Z = Matrix(model.matrix(formula(sprintf('~0 + %s',re)),data),sparse = TRUE)
	  Z[,paste0(re,levels(data[[re]]))]
	})
	names(Z_matrices) = RE_names
	r_RE = sapply(Z_matrices,function(x) ncol(x))

	Z = do.call(cbind,Z_matrices)

	# covariance matrices from random effects
	#   Three options for each random effect:
	#     1) pass an A matrix
	#     2) pass an A_inv matrix (say from a call to inverseA)
	#     3) NULL (or un-named). Will construct as an identity matrix
	#     options 1 and 2: the provided matrix must have rows/columns for each individual in data,
	#       but can have extra rows (and should for A_inv) if ancestors are known

	  # function to ensure that covariance matrices are sparse and symmetric
	fix_A = function(x) forceSymmetric(drop0(x,tol = 1e-10))

	  # construct A matrices for each random effect
	A_mats = lapply(RE_names,function(re) {
	  if(re %in% names(A_mats)) {
	    A = A_mats[[re]]
	  } else if(re %in% names(A_inv_mats)){
	    A = solve(A_inv_mats[[re]])
	    rownames(A) = rownames(A_inv_mats[[re]])
	  } else{
	    A = Diagonal(ncol(Z_matrices[[re]]))
	    rownames(A) = levels(data[[re]])
	  }
	  index = match(sub(re,'',colnames(Z_matrices[[re]])),rownames(A))  # A must have rownames
	  fix_A(A[index,index])
	})
	names(A_mats) = RE_names

	  # construct A_inverse matrices for each random effect
	Ai_mats = lapply(RE_names,function(re) {
	  if(re %in% names(A_inv_mats)){
	    Ai = A_inv_mats[[re]]
	    index = match(sub(re,'',colnames(Z_matrices[[re]])),rownames(Ai))  # Ai must have rownames
	    Ai = Ai[index,index]
	  } else {
	    Ai = solve(A_mats[[re]])
	  }
	  Ai = fix_A(Ai)
	  Ai
	})
	names(Ai_mats) = RE_names
	  # cholesky decompositions (L'L) of each A_inverse matrix
	chol_Ai_mats = lapply(Ai_mats,chol)

	  # table of possible h2s for each random effect
	  #   These are percentages of the total residual variance accounted for by each random effect
	  #   Each column is a set of percentages, the sum of which must be less than 1 (so that Ve is > 0)
	  # can specify different levels of granularity for each random effect
	h2_divisions = run_parameters$h2_divisions
	if(length(h2_divisions) < n_RE){
	  if(length(h2_divisions) != 1) stop('Must provide either 1 h2_divisions parameter, or 1 for each random effect')
	  h2_divisions = rep(ceiling(h2_divisions^(1/n_RE)),n_RE)  # evenly divide divisions among random effects
	}
	if(is.null(names(h2_divisions))) {
	  names(h2_divisions) = RE_names
	}
	h2s_matrix = expand.grid(lapply(RE_names,function(re) seq(0,1,length = h2_divisions[[re]]+1)))
	colnames(h2s_matrix) = RE_names
	h2s_matrix = t(h2s_matrix[rowSums(h2s_matrix) < 1,,drop=FALSE])
	colnames(h2s_matrix) = NULL

	data_matrices = list(
	  Y          = Y,
	  Y_missing  = Y_missing,
	  X          = X,
	  X_F        = X_F,
	  Z_matrices = Z_matrices,
	  Z          = Z,
	  h2s_matrix = h2s_matrix
	)

	run_variables = list(
	  p      = p,
	  n      = n,
	  r_RE   = r_RE,
	  b      = b,
	  b_F    = b_F,
	  Mean_Y = Mean_Y,
	  VY     = VY
	)


	# ----------------------------- #
	# ----- re-formulate priors --- #
	# ----------------------------- #
	if(any(sapply(priors, function(x) {try({return(x$nu <= 2)},silent=T);return(FALSE)}))) stop('priors nu must be > 2')
	  # fixed effects
	priors$fixed_prec_shape = with(priors$fixed_var,V * nu)
	priors$fixed_prec_rate  = with(priors$fixed_var,nu - 2)
	  # total precision
	priors$tot_Y_prec_shape = with(priors$tot_Y_var,V * nu)
	priors$tot_Y_prec_rate  = with(priors$tot_Y_var,nu - 2)
	priors$tot_F_prec_shape = with(priors$tot_F_var,V * nu)
	priors$tot_F_prec_rate  = with(priors$tot_F_var,nu - 2)
	  # delta: column shrinkage of Lambda
	priors$delta_1_shape    = with(priors$delta_1,V * nu)
	priors$delta_1_rate     = with(priors$delta_1,nu - 2)
	priors$delta_2_shape    = with(priors$delta_2,V * nu)
	priors$delta_2_rate     = with(priors$delta_2,nu - 2)


	# ----------------------------- #
	# -- create BSFG_state object - #
	# ----------------------------- #

	BSFG_state = list(
	  data_matrices  = data_matrices,
	  priors         = priors,
	  run_parameters = run_parameters,
	  run_variables  = run_variables,
	  setup          = setup,
	  traitnames     = traitnames
	)
	class(BSFG_state) = append(class(BSFG_state),c('BSFG_state',run_parameters$sampler))

	BSFG_state = initialize_BSFG(BSFG_state, A_mats, chol_Ai_mats,verbose)

	return(BSFG_state)
}


initialize_BSFG = function(BSFG_state,...){
  UseMethod("initialize_BSFG",BSFG_state)
}

#' Run BSFG Gibbs sampler
#'
#' Run MCMC chain for a specified number of iterations
sample_BSFG = function(BSFG_state,...){
  UseMethod("sample_BSFG",BSFG_state)
}

#' Print more detailed statistics on current BSFG state
#'
#' Print more detailed statistics on current BSFG state
#' @seealso \code{\link{BSFG_control}}, \code{\link{sample_BSFG}}, \code{\link{BSFG_init}},
#'   \code{\link{print.BSFG_state}}, \code{\link{plot.BSFG_state}}
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

#' Print statistics on current BSFG state
#'
#' Print statistics on current BSFG state
#' @seealso \code{\link{BSFG_control}}, \code{\link{sample_BSFG}}, \code{\link{BSFG_init}},
#'   \code{\link{summary.BSFG_state}}, \code{\link{plot.BSFG_state}}
print.BSFG_state = function(BSFG_state){
  with(BSFG_state,{
    cat(
      c(sprintf('\n Current iteration: %d, Posterior_samples: %d \n',current_state$nrun,Posterior$sp_num)),
      c(sprintf('Current factor dimension: %d factors \n',ncol(current_state$Lambda))),
      c(sprintf('Total time: %s \n\n',format(current_state$total_time)))
    )
  })
}

#' Make plots of current BSFG state
#'
#' Make plots of current BSFG state
#' @seealso \code{\link{BSFG_control}}, \code{\link{sample_BSFG}}, \code{\link{BSFG_init}},
#'   \code{\link{print.BSFG_state}}, \code{\link{summary.BSFG_state}}
plot.BSFG_state = function(BSFG_state){
  if(BSFG_state$run_parameters$simulation){
    plot_diagnostics_simulation(BSFG_state)
  } else{
    plot_diagnostics(BSFG_state)
  }
}
