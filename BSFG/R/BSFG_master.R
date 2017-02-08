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
BSFG_init = function(Y, fixed, random, data, priors, run_parameters, A_mats = NULL, A_inv_mats = NULL,sampler = 'fast_BSFG',
                                    ncores = 1,simulation = F,setup = NULL,verbose=T) {


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

	if(n_RE > 1 && sampler == 'fast_BSFG'){
	  print(sprintf('%d random effects. Using "general_BSFG" sampler',length(RE_names)))
	  sampler = 'general_BSFG'
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
	class(BSFG_state) = append(class(BSFG_state),c('BSFG_state',sampler))

	BSFG_state = initialize_BSFG(BSFG_state, A_mats, chol_Ai_mats,verbose)

	return(BSFG_state)
}


initialize_BSFG = function(BSFG_state,...){
  UseMethod("initialize_BSFG",BSFG_state)
}

sample_BSFG = function(BSFG_state,...){
  UseMethod("sample_BSFG",BSFG_state)
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
  if(BSFG_state$run_parameters$simulation){
    draw_simulation_diagnostics(BSFG_state)
  } else{
    draw_results_diagnostics(BSFG_state)
  }
}
