#' Set BSFG run parameters
#'
#' Function to create run_parameters list for initializing BSFG model
#'
#' @param sampler specify the sampler to use. fast_BSFG is much faster, but only allows one random
#'   effect. If more are specified in \code{BSFG_init}, this is switched to general_BSFG.
#' @param Posterior_folder folder to save posterior sample chunk files
#' @param fixed_factors should the fixed effect model be applied to factors as well as the factor
#'   residuals? If so, the first column of the design matrix is dropped.
#' @param simulation Is this a fit to simulated data? If so, a setup list will be expected providing
#'   the true values
#' @param save_Eta Should posterior samples of Eta be saved?
#' @param scale_Y Should the Y values be centered and scaled? Recommend, except for simulated data.
#' @param b0 parameter of the \code{update_k} function. See Bhattacharya and Dunson 2011
#' @param b1 parameter of the \code{update_k} function. See Bhattacharya and Dunson 2011
#' @param epsilon parameter of the \code{update_k} function. Smallest \eqn{\lambda_{ij}} that is
#'   considered "large", signifying a factor should be kept. See Bhattacharya and Dunson 2011
#' @param prop proportion of \eqn{\lambda{ij}} elements in a column of \eqn{\Lambda} that must be smaller than
#'   \code{epsilon} before factor is dropped. See Bhattacharya and Dunson 2011
#' @param kinit initial number of factors
#' @param h2_divisions A scalar or vector of length equal to number of random effects. In BSFG, random
#'   effects are re-scaled as percentages of the total variation. Then a discrete prior spanning [0,1)
#'   with \code{h2_divisions} equally spaced values is constructred for each variance component. If
#'   \code{h2_divisions} is a scalar, the prior for each variance component has this number of divisions.
#'   In the joint prior over all variance components, combinations of variance components with total variance != 1
#'   are assigned a prior of zero and ignored.
#' @param burn burnin length of the MCMC chain
#' @param thin thinning rate of the MCMC chain
#' @seealso \code{\link{BSFG_init}}, \code{\link{sample_BSFG}}, \code{\link{print.BSFG_state}}
#'
BSFG_control = function(sampler = c('fast_BSFG','general_BSFG'),Posterior_folder = 'Posterior',
                        fixed_factors = c(T,F),save_Eta = F, simulation = c(F,T),scale_Y = c(T,F),
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
#' @param model RHS of a model. The syntax is similar to \link{lmer}.
#'     Random effects are specified by (1+factor | group), with the left side of the '|' a design
#'     matrix, and the right side a random effect factor (group). For each random effect factor, a
#'     covariance matrix (\code{A_mats}) or precision matrix (\code{A_inv_mats}) can be provided.
#'     Unlike in \code{lmer}, each variable or covariate in the design matrix is given an
#'     independent random effect (ie no covariance among random effects is modeled), so two bars '||'
#'     gives an identical model to one bar.
#'     Note: the speed of the model will decrease dramatically with the number of random effects (multiplied
#'     by h2_divisions)
#' @param data data.frame with n rows containing columns corresponding to the fixed and random
#'   effects
#' @param priors list providing hyperparameters for the model priors. This must include: for \code{fixed_var},
#'     \code{tot_Eta_var}, and \code{tot_F_var}, a list of \code{V} and \code{nu}, which specify an inverse-gamma
#'     distribution as in MCMCglmm. For \code{delta_1} and \code{delta_2} a list of \code{shape} and \rate{rate} for
#'     a gamma distribution. For \code{Lambda_df}, the degrees of freedom of the implied t-distribution.
#'     Discrete priors on \code{h2_priors_factors} and \code{h2_priors_resids} are also required for \code{sample_BSFG},
#'     but can be appended to \code{BSFG_state$priors} after running \code{BSFG_init}. This is often easier if
#'     the total number of h2 divisions is not known beforehand (ie when multiple random effects are used).
#' @param run_parameters list providing various parameters for the model run. See \link{BSFG_control}
#' @param A_mats list of covariance matrices for random effects. If none provided (and none provided
#'   for A_inv_mats) for any of the random effects, A is assumed to be the identity.
#' @param A_inv_mats list of precision matrices for random effects. If none provided (and none
#'   provided for A_mats) for any of the random effects, A is assumed to be the identity.
#' @param data_model either a character vector identifying a provided data_model
#'     (ex. 'missing_data' calls \link{missing_data_model}),
#'     or a function modeled after \code{missing_data_model} (see code by typing this into the console)
#'     that draws posterior samples of Eta conditional on
#'     the observations (Y) and the current state of the BSFG model (current_state). The function should have
#'     the same form as \link{missing_data}, and must return Eta even if current_state is NULL
#' @param data_model_parameters a list of parameters necessary for computing \code{data_model}.
#'     Ex. Y_missing for the missing_data model.
#'     If \code{data_model == 'missing_data'}, Y_missing is computed automatically.
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
BSFG_init = function(Y, model, data, priors, run_parameters, A_mats = NULL, A_inv_mats = NULL,
                     data_model = 'missing_data', data_model_parameters = NULL,
                     sampler = c('fast_BSFG','general_BSFG'), ncores = detectCores(),simulation = c(F,T),setup = NULL,verbose=T) {

  # ----------------------------- #
  # ---- build model matrices --- #
  # ----------------------------- #

	# model dimensions
	n = nrow(data)
	p_Y = ncol(Y)
	traitnames = colnames(Y)

	Y_missing = Matrix(is.na(Y))

	# scale Y
	if(run_parameters$scale_Y){
	  Mean_Y = colMeans(Y,na.rm=T)
	  VY = apply(Y,2,var,na.rm=T)
	  Y = sweep(Y,2,Mean_Y,'-')
	  Y = sweep(Y,2,sqrt(VY),'/')
	} else {
	  Mean_Y = rep(0,p_Y)
	  VY = rep(1,p_Y)
	}

	# use data_model to get dimensions, names of Y
	if(is.character(data_model) && data_model == 'missing_data') {
	  data_model = missing_data_model
	  data_model_parameters = list(Y_missing = Matrix(is.na(Y)))
	}
	Eta = data_model(Y,data_model_parameters)$Eta
	p = ncol(Eta)
	traitnames = colnames(Y)

	# check that all terms in model are in data
	if(!all(all.vars(model) %in% colnames(data))) {
	  terms = all.vars(model)
	  missing_terms = terms[!terms %in% colnames(data)]
	  stop(sprintf('terms %s missing from data',paste(missing_terms,sep=', ')))
	}

	# build X from fixed model
	  # for residuals
	X = model.matrix(nobars(model),data)
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

	# use lme4 functions to parse random effects
	#  note: correlated random effects are not allowed. Will convert to un-correlated REs
	RE_terms = mkReTrms(findbars(model),data)  # extracts terms and builds Zt matrices

	Z_matrices = list()
	RE_names = c()
	RE_covs = c()
	for(i in 1:length(RE_terms$cnms)){
	  term = names(RE_terms$cnms)[i]
	  n_factors = length(RE_terms$cnms[[i]])  # number of factors for this grouping factor
	  RE_covs = c(RE_covs,rep(term,n_factors)) # name of covariance for this grouping factor (for A_mats)
	  if(sum(names(RE_terms$cnms) == term) > 1 || n_factors > 1) {
	    RE_names = c(RE_names,paste(term,RE_terms$cnms[[i]],sep='.'))  # name of variance component
	  } else{
	    RE_names = c(RE_names,term)
	  }
	  combined_Zt = RE_terms$Ztlist[[i]] # combined matrix of Zt matrices for this grouping factor
	  # split combined matrix into a list of Z matrices (transposed)
	  Zs_term = tapply(1:nrow(combined_Zt),gl(n_factors,1,nrow(combined_Zt),labels = RE_terms$cnms[[term]]),function(x) t(combined_Zt[x,]))
	  Z_matrices = c(Z_matrices,Zs_term)
	}
	names(Z_matrices) = RE_names
	Z = do.call(cbind,Z_matrices)

	n_RE = length(RE_names)
	r_RE = sapply(Z_matrices,function(x) ncol(x))

	if(n_RE > 1 && run_parameters$sampler == 'fast_BSFG'){
	  print(sprintf('%d random effects. Using \"general_BSFG\" sampler',length(RE_names)))
	  run_parameters$sampler = 'general_BSFG'
	}

	  # function to ensure that covariance matrices are sparse and symmetric
	fix_A = function(x) forceSymmetric(drop0(x,tol = 1e-10))

	  # construct A matrices for each random effect
	if(is.null(A_mats)) A_mats = list()
	for(i in 1:length(RE_names)){
	  re = RE_covs[i]
	  re_name = RE_names[i]
	  if(re %in% names(A_mats)) {
	    A = A_mats[[re]]
	  } else if(re %in% names(A_inv_mats)){
	    A_mats[[re]] = solve(A_inv_mats[[re]])
	    rownames(A_mats[[re]]) = rownames(A_inv_mats[[re]])
	    A = A_mats[[re]]
	  } else{
	    A_mats[[re]] = Diagonal(ncol(Z_matrices[[re_name]]))
	    rownames(A_mats[[re]]) = levels(data[[re]])
	    A = A_mats[[re]]
	  }
	  if(is.null(rownames(A))) stop('A must have rownames')
	  index = match(colnames(Z_matrices[[re_name]]),rownames(A)) # A must have rownames
	  if(any(is.na(index))) stop(sprintf('levels missing from covarince of random effect %s',re))
	  stopifnot(length(index) == ncol(Z_matrices[[re_name]]))
	  A_mats[[re_name]] = fix_A(A[index,index])
	}
	A_mats = A_mats[RE_names]

	  # construct A_inverse matrices for each random effect
	if(is.null(A_inv_mats)) A_inv_mats = list()
	Ai_mats = list()
	for(i in 1:length(RE_names)){
	  re = RE_covs[i]
	  re_name = RE_names[i]
	  if(re %in% names(A_inv_mats)){
	    Ai = A_inv_mats[[re]]
	  } else {
	    A_inv_mats[[re_name]] = solve(A_mats[[re_name]])
	    rownames(A_inv_mats[[re_name]]) = rownames(A_mats[[re_name]])
	    Ai = A_inv_mats[[re_name]]
	  }
	  if(is.null(rownames(Ai))) stop('A_inv must have rownames')
	  index = match(colnames(Z_matrices[[re_name]]),rownames(Ai)) # iA must have rownames
	  Ai_mats[[re_name]] = fix_A(Ai[index,index])
	}
	# names(Ai_mats) = RE_names
	  # cholesky decompositions (L'L) of each A_inverse matrix
	chol_Ai_mats = lapply(Ai_mats,chol)

	  # table of possible h2s for each random effect
	  #   These are percentages of the total residual variance accounted for by each random effect
	  #   Each column is a set of percentages, the sum of which must be less than 1 (so that Ve is > 0)
	  # can specify different levels of granularity for each random effect
	h2_divisions = run_parameters$h2_divisions
	if(length(h2_divisions) < n_RE){
	  if(length(h2_divisions) != 1) stop('Must provide either 1 h2_divisions parameter, or 1 for each random effect')
	  # h2_divisions = rep(ceiling(h2_divisions^(1/n_RE)),n_RE)  # evenly divide divisions among random effects
	  h2_divisions = rep(h2_divisions,n_RE)
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

	run_parameters$data_model = data_model
	run_parameters$data_model_parameters = data_model_parameters


	# ----------------------------- #
	# ----- re-formulate priors --- #
	# ----------------------------- #
	if(any(sapply(priors, function(x) {try({return(exists(x$nu) && x$nu <= 2)},silent=T);return(FALSE)}))) stop('priors nu must be > 2')
	  # fixed effects
	priors$fixed_prec_rate   = with(priors$fixed_var,V * nu)
	priors$fixed_prec_shape  = with(priors$fixed_var,nu - 1)
	  # total precision
	priors$tot_Eta_prec_rate   = with(priors$tot_Y_var,V * nu)
	priors$tot_Eta_prec_shape  = with(priors$tot_Y_var,nu - 1)
	priors$tot_F_prec_rate     = with(priors$tot_F_var,V * nu)
	priors$tot_F_prec_shape    = with(priors$tot_F_var,nu - 1)
	  # delta: column shrinkage of Lambda
	# priors$delta_1_rate    = with(priors$delta_1,V * nu)
	# priors$delta_1_shape   = with(priors$delta_1,nu - 1)
	# priors$delta_2_rate    = with(priors$delta_2,V * nu)
	# priors$delta_2_shape   = with(priors$delta_2,nu - 1)
	priors$delta_1_rate    = priors$delta_1$rate
	priors$delta_1_shape   = priors$delta_1$shape
	priors$delta_2_rate    = priors$delta_2$rate
	priors$delta_2_shape   = priors$delta_2$shape


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

	BSFG_state = initialize_BSFG(BSFG_state, A_mats, chol_Ai_mats,verbose=verbose,ncores=ncores)

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
      c(sprintf('Current iteration: %d, Posterior_samples: %d \n',current_state$nrun,Posterior$total_samples)),
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
      c(sprintf('\n Current iteration: %d, Posterior_samples: %d \n',current_state$nrun,Posterior$total_samples)),
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
plot.BSFG_state = function(BSFG_state,file = 'diagnostics_polts.pdf'){
  pdf(file)
  if(BSFG_state$run_parameters$simulation){
    plot_diagnostics_simulation(BSFG_state)
  } else{
    plot_diagnostics(BSFG_state)
  }
  dev.off()
}

#' Sample missing data
#'
#' Function to sample missing data given model parameters.
#'
#' This is also a template for \code{data_model} functions.
#'
#'    The function should draw a posterior sample for Eta (the (potentially) latent
#'    data on the linear scale) given the factor model and the observed data. It returns a list
#'    with Eta and any other variables associated with the data_model.
#'    These variables are added to current_state, but are not currently saved in Posterior
#'
#'    Initial values, hyperparameters, and any helper functions / parameters can be passed
#'    in \code{data_model_parameters}.
#'
#'    When run with an empty \code{current_state} and \code{data_matrices == NULL}, it should return
#'    a matrix Eta of the correct dimensions, but the values are unimportant.
#'
#' @param Y data matrix n_Y x p_Y
#' @param data_model_parameters List of parameters necessary for the data model.
#'      Here, a Matrix of coordinates of NAs in Y
#' @param current_state current_state from BSFG_state
#' @param data_matrices from BSFG_state
#' @return list of data_model variables including:
#' @return Eta matrix of latent data: n x p
missing_data_model = function(Y,data_model_parameters,current_state = list(),data_matrices = NULL){
  new_variables = c('Eta')
  data_model_state = with(c(data_model_parameters,data_matrices,current_state),{
        Eta = Y
        if(sum(Y_missing)  > 0){
          n = nrow(Y)
          p = ncol(Y)
          if(length(current_state) == 0) {
            Eta_mean = matrix(0,n,p)
            resids = matrix(rnorm(n*p),n,p)
          } else{
            Eta_mean = X %*% B + F %*% t(Lambda) + Z %*% E_a
            resid_Eta_prec = tot_Eta_prec / (1-resid_h2)
            resids = matrix(rnorm(p*n,0,sqrt(1/resid_Eta_prec)),nr = n,nc = p,byrow=T)
          }
          missing_indices = which(Y_missing)
          Eta[missing_indices] = Eta_mean[missing_indices] + resids[missing_indices]
        }
        return(list(Eta = as.matrix(Eta)))
    })
  return(data_model_state[new_variables])
}
