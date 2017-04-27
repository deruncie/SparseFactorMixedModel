#' Set BSFG run parameters
#'
#' Function to create run_parameters list for initializing BSFG model
#'
#' @param sampler specify the sampler to use. fast_BSFG is often faster, but only allows one random
#'   effect. If more are specified in \code{BSFG_init}, this is switched to general_BSFG.
#' @param Posterior_folder path to folder to save posterior samples. Samples of each parameter
#'     are saved in chuncks to limit memory requirements.
#' @param simulation Is this a fit to simulated data? If so, a setup list will be expected providing
#'   the true values
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
#' @param h2_step_size Either NULL, or a scaler in the range (0,1] giving specifying the range of h2 values for a Metropolis-Hastings
#'   update step for each h2 parameter vector. If NULL, h2's will be sampled based on the marginal probability
#'   over all possible h2 vectors. If a scalar, a Metropolis-Hastings update step will be used for each h2 vector.
#'   The trail value will be selected uniformly from all possible h2 vectors within this Euclidean distance from the current vector.
#' @param drop0_tol A scalar giving the a tolerance for the \code{drop0()} function that will be applied
#'     to various symmetric (possibly) sparse matrices to try to fix numerical errors and increase sparsity.
#' @param K_eigen_tol A scalar giving the minimum eigenvalue of a K matrix allowed. During pre-processing,
#'     eigenvalues of each K matrix will be calculated using \code{svd(K)}. Only eigenvectors of K with corresponding eigenvalues
#'     greater than this value will be kept. If smaller eigenvalues exist, the model will be transformed
#'     to reduce the rank of K, by multiplying Z by the remaining eigenvectors of K. This transformation
#'     is undone before posterior samples are recorded, so posterior samples of \code{U_F} and \code{U_R} are
#'     untransformed.
#' @param burn burnin length of the MCMC chain
#' @param thin thinning rate of the MCMC chain
#' @seealso \code{\link{BSFG_init}}, \code{\link{sample_BSFG}}, \code{\link{print.BSFG_state}}
#'
BSFG_control = function(sampler = c('fast_BSFG','general_BSFG'),Posterior_folder = 'Posterior',
                        simulation = c(F,T),scale_Y = c(T,F),
                        b0 = 1, b1 = 0.0005, epsilon = 1e-1, prop = 1.00,
                        k_init = 20, h2_divisions = 100, h2_step_size = NULL,
                        drop0_tol = 1e-14, K_eigen_tol = 1e-10,
                        burn = 100,
                        thin = 2) {

  all_args = lapply(formals(),function(x) eval(x)[1])
  passed_args = as.list(match.call())[-1]
  all_args[names(passed_args)] = passed_args
  return(all_args)
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
#'     covariance matrix (\code{K_mats}) or precision matrix (\code{K_inv_mats}) can be provided.
#'     Unlike in \code{lmer}, each variable or covariate in the design matrix is given an
#'     independent random effect (ie no covariance among random effects is modeled), so two bars '||'
#'     gives an identical model to one bar.
#'     Note: the speed of the model will decrease dramatically with the number of random effects (multiplied
#'     by h2_divisions)
#' @param data data.frame with n rows containing columns corresponding to the fixed and random
#'   effects
#' @param factor_model_fixed Fixed effect model formula specific to the latent factors. Optional. If NULL,
#'     the fixed effects for the latent factors will be the same as for Eta, except that the intercept
#'     will be removed. If a formula is provided (no random effects allowed), the model will be applied to
#'     each factor, again with the intercept dropped (or set to zero). Note: Random effects in \code{model}
#'     are applied to both Eta and F.
#' @param priors list providing hyperparameters for the model priors. This must include: for \code{fixed_var},
#'     \code{tot_Eta_var}, and \code{tot_F_var}, a list of \code{V} and \code{nu}, which specify an inverse-gamma
#'     distribution as in MCMCglmm. For \code{delta_1} and \code{delta_2} a list of \code{shape} and \code{rate} for
#'     a gamma distribution. For \code{Lambda_df}, the degrees of freedom of the implied t-distribution.
#'     Discrete priors on \code{h2_priors_factors} and \code{h2_priors_resids} are also required for \code{sample_BSFG},
#'     but can be appended to \code{BSFG_state$priors} after running \code{BSFG_init}. This is often easier if
#'     the total number of h2 divisions is not known beforehand (ie when multiple random effects are used).
#' @param run_parameters list providing various parameters for the model run. See \link{BSFG_control}
#' @param K_mats list of covariance matrices for random effects. If none provided (and none provided
#'   for K_inv_mats) for any of the random effects, K is assumed to be the identity.
#' @param K_inv_mats list of precision matrices for random effects. If none provided (and none
#'   provided for K_mats) for any of the random effects, K is assumed to be the identity.
#' @param data_model either a character vector identifying a provided data_model
#'     (ex. 'missing_data' calls \link{missing_data_model}),
#'     or a function modeled after \code{missing_data_model} (see code by typing this into the console)
#'     that draws posterior samples of Eta conditional on
#'     the observations (Y) and the current state of the BSFG model (current_state). The function should have
#'     the same form as \link{missing_data}, and must return Eta even if current_state is NULL
#' @param data_model_parameters a list of parameters necessary for computing \code{data_model}.
#'     Ex. Y_missing for the missing_data model.
#'     If \code{data_model == 'missing_data'}, Y_missing is computed automatically.
#' @param X_resid a design matrix. Alternative method for specifying fixed effects for \code{model}.
#'     Care must be taken with the intercept if both a fixed effect model is specified for in both
#'     \code{model} and \code{X_resid}.
#' @param X_factor a design matrix. Alternative method for specifying \code{factor_model_fixed}
#' @param cis_genotypes a list of design matrices of length \code{p} (ie number of columns of \code{Eta})
#'     This is used to specify trait-specific fixed effects, such a cis-genotypes
#' @param posteriorSample_params A character vector giving names of parameters to save all posterior samples
#' @param posteriorMean_params A character vector giving names of parameters to save only the posterior mean.
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
#' @seealso \code{\link{BSFG_control}}, \code{\link{sample_BSFG}}, \code{\link{print.BSFG_state}}, \code{\link{plot.BSFG_state}}#'
BSFG_init = function(Y, model, data, factor_model_fixed = NULL, priors, run_parameters, K_mats = NULL, K_inv_mats = NULL,
                     data_model = 'missing_data', data_model_parameters = NULL, X_resid = NULL, X_factor = NULL, cis_genotypes = NULL,
                     posteriorSample_params = c('Lambda','U_F','F','delta','tot_F_prec','F_h2','tot_Eta_prec','resid_h2', 'B', 'B_F','U_R','tau_B','tau_B_F','cis_effects'),
                     posteriorMean_params = c(),
                     ncores = detectCores(),setup = NULL,verbose=T) {

  # ----------------------------- #
  # ---- build model matrices --- #
  # ----------------------------- #

	# model dimensions
	n = nrow(data)
	p_Y = ncol(Y)
	traitnames = colnames(Y)

	Y_missing = Matrix(is.na(Y))

	# scale Y
	if(!is(Y,'matrix'))	Y = as.matrix(Y)
	if(run_parameters$scale_Y){
	  Mean_Y = colMeans(Y,na.rm=T)
	  VY = apply(Y,2,var,na.rm=T)
	  Y = sweep(Y,2,Mean_Y,'-')
	  Y = sweep(Y,2,sqrt(VY),'/')
	} else {
	  p_Y = dim(Y)[2]
	  Mean_Y = rep(0,p_Y)
	  VY = rep(1,p_Y)
	}

	# use data_model to get dimensions, names of Y
	if(is.character(data_model) && data_model == 'missing_data') {
	  data_model = missing_data_model
	  data_model_parameters = list(Y_missing = Matrix(is.na(Y)))
	}
	if(is.character(data_model) && data_model == 'voom_RNAseq') {
	  data_model = voom_model()
	  data_model_parameters$Y_std = Y * data_model_parameters$prec_Y
	}
	data_model_state = data_model(Y,data_model_parameters)
	Eta = data_model_state$state$Eta
	p = ncol(Eta)
	traitnames = colnames(Y)

	# if factor_model_fixed not specified, use fixed effects from model for both
	if(is.null(factor_model_fixed)) factor_model_fixed = nobars(model)

	# check that there are no random effects specifed in factor_model_fixed
	if(length(findbars(factor_model_fixed))) stop('Do not specify random effects in `factor_model_fixed`. Random effects for factors taken from `model`')

	# check that all terms in models are in data
	terms = c(all.vars(model),all.vars(factor_model_fixed))
	if(!all(terms %in% colnames(data))) {
	  missing_terms = terms[!terms %in% colnames(data)]
	  stop(sprintf('terms %s missing from data',paste(missing_terms,sep=', ')))
	}

	# -------- Fixed effects ---------- #

	# build X from fixed model
	  # for Eta
	X = model.matrix(nobars(model),data)
	X = cbind(X, X_resid) # add in X_resid if provided.
	b = ncol(X)

  # for F
	  # note we drop a column to force intercept to be zero
	X_F = model.matrix(factor_model_fixed,data)[,-1,drop = FALSE]
	X_F = cbind(X_F,X_factor)
	b_F = ncol(X_F)

	# -------- cis genotypes ---------- #
	if(is.null(cis_genotypes)){
	  cis_effects_index = NULL
	} else{
	  cis_effects_index = do.call(c,lapply(1:length(cis_genotypes),function(j) rep(j,ncol(cis_genotypes[[j]]))))
	}

	# -------- Random effects ---------- #
	# ensure that only K or K_inv provided for each random effect
	# check that all IDs from data are in rownames of K or K_inv
	# add IDs from K or K_inv to levels data[[re]] for IDs not in data

	RE_levels = list() # a list of levels for each of the random effects
	if(is.null(K_mats)) K_mats = list()
	for(re in names(K_mats)) {
	  # check that K is a matrix, then convert to Matrix
	  if(is.data.frame(K_mats[[re]])) K_mats[[re]] = as.matrix(K_mats[[re]])
	  if(is.matrix(K_mats[[re]])) K_mats[[re]] = Matrix(K_mats[[re]],sparse=T)
	  if(is.null(rownames(K_mats[[re]]))) stop(sprintf('K %s must have rownames',re))
	  if(re %in% names(K_inv_mats)) stop(sprintf('Both K and Kinv provided for %s. Please provide only one for each random effect',re))
	  RE_levels[[re]] = rownames(K_mats[[re]])
	}
	for(re in names(K_inv_mats)) {
	  # check that K_inv is a matrix, then convert to Matrix
	  if(is.data.frame(K_inv_mats[[re]])) K_inv_mats[[re]] = as.matrix(K_inv_mats[[re]])
	  if(is.matrix(K_inv_mats[[re]])) K_inv_mats[[re]] = Matrix(K_inv_mats[[re]],sparse=T)

	  if(is.null(rownames(K_inv_mats[[re]]))) stop(sprintf('K_inv %s must have rownames',re))
	  RE_levels[[re]] = rownames(K_inv_mats[[re]])
	}
	for(re in names(RE_levels)){
	  if(!re %in% colnames(data)) stop(sprintf('Column "%s" required in data',re))
	  data[[re]] = as.factor(data[[re]]) # ensure 'data[[re]]' is a factor
	  if(!all(levels(data[[re]]) %in% RE_levels[[re]])) stop(sprintf('Levels of random effect %s missing from provided %s',re,c('K_inv','K')[(re %in% names(K_mats))+1]))
	  data[[re]] = factor(data[[re]],levels = RE_levels[[re]]) # add levels to data[[re]]
	}

	# use lme4 functions to parse random effects
	#  note: correlated random effects are not allowed. Will convert to un-correlated REs
	RE_terms = mkReTrms(findbars(model),data,drop.unused.levels = FALSE)  # extracts terms and builds Zt matrices

	# construct Z matrices for each Random Effect
	Z_matrices = list()
	RE_names = c()
	RE_covs = c()
	for(i in 1:length(RE_terms$cnms)){
	  term = names(RE_terms$cnms)[i]
	  n_factors = length(RE_terms$cnms[[i]])  # number of factors for this grouping factor
	  RE_covs = c(RE_covs,rep(term,n_factors)) # name of covariance for this grouping factor (for K_mats)
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

	  # function to ensure that covariance matrices are sparse and symmetric
	fix_K = function(x) forceSymmetric(drop0(x,tol = run_parameters$drop0_tol))

	# construct K matrices for each random effect
	    # if K is PSD, find K = USVt and set K* = S and Z* of ZU
	svd_Ks = lapply(K_mats,function(K) svd(K))
	RE_U_matrices = list()

	for(i in 1:length(RE_names)){
	  re = RE_covs[i]
	  re_name = RE_names[i]
	  if(re %in% names(svd_Ks)) {
	    svd_K = svd_Ks[[re]]
	    r_eff = sum(svd_K$d > run_parameters$K_eigen_tol)  # truncate eigenvalues at this value
	    # if need to use reduced rank model, then use the svd of K in place of K and merge U into Z
	    # otherwise, use original K, set U = Diagonal(1,r)
	    if(r_eff < length(svd_K$d)) {
  	    K = Diagonal(r_eff,svd_K$d[1:r_eff])
  	    U = Matrix(svd_K$u[,1:r_eff])
	    } else{
	      K = K_mats[[re]]
	      U = Diagonal(nrow(K),1)
	    }
	  } else if(re %in% names(K_inv_mats)){
	    K_mats[[re]] = solve(K_inv_mats[[re]])
	    rownames(K_mats[[re]]) = rownames(K_inv_mats[[re]])
	    K = K_mats[[re]]
	    U = Diagonal(nrow(K),1)
	  } else{
	    K_mats[[re]] = Diagonal(ncol(Z_matrices[[re_name]]))
	    rownames(K_mats[[re]]) = levels(data[[re]])
	    K = K_mats[[re]]
	    U = Diagonal(nrow(K),1)
	  }
	  K_mats[[re_name]] = fix_K(K)
	  RE_U_matrices[[re_name]] = U
	  Z_matrices[[re_name]] = Z_matrices[[re_name]] %*% U
	}
	K_mats = K_mats[RE_names]
	# Fix Z_matrices based on PSD K's
	Z = do.call(cbind,Z_matrices[RE_names])
	Z = as(do.call(cbind,Z_matrices[RE_names]),'dgCMatrix')
	# The following matrix is used to transform random effects back to the original space had we sampoled from the original (PSD) K.
	if(length(RE_names) > 1) {
	  U_svd = do.call(bdiag,RE_U_matrices[RE_names])
	} else{
	  U_svd = RE_U_matrices[[1]]
	}
	r_RE = sapply(Z_matrices,function(x) ncol(x))  # re-calculate

	  # construct K_inverse matrices for each random effect
	if(is.null(K_inv_mats)) K_inv_mats = list()
	Ki_mats = list()
	for(i in 1:length(RE_names)){
	  re = RE_covs[i]
	  re_name = RE_names[i]
	  if(re %in% names(K_inv_mats)){
	    Ki = K_inv_mats[[re]]
	  } else {
	    K_inv_mats[[re_name]] = solve(K_mats[[re_name]])
	    rownames(K_inv_mats[[re_name]]) = rownames(K_mats[[re_name]])
	    Ki = K_inv_mats[[re_name]]
	  }
	  Ki_mats[[re_name]] = fix_K(Ki)
	}
	# names(Ki_mats) = RE_names
	  # cholesky decompositions (L'L) of each K_inverse matrix
	chol_Ki_mats = lapply(Ki_mats,chol)


	n_RE = length(RE_names)

	if(n_RE > 1 && run_parameters$sampler == 'fast_BSFG'){
	  cat(sprintf('%d random effects. Using \"general_BSFG\" sampler\n',length(RE_names)))
	  run_parameters$sampler = 'general_BSFG'
	}


	  # table of possible h2s for each random effect
	  #   These are percentages of the total residual variance accounted for by each random effect
	  #   Each column is a set of percentages, the sum of which must be less than 1 (so that Ve is > 0)
	  # can specify different levels of granularity for each random effect
	h2_divisions = run_parameters$h2_divisions
	if(length(h2_divisions) < n_RE){
	  if(length(h2_divisions) != 1) stop('Must provide either 1 h2_divisions parameter, or 1 for each random effect')
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
	  U_svd      = U_svd,  # matrix necessary to back-transform U_F and U_R (U*U_F and U*U_R) to get original random effects
	  h2s_matrix = h2s_matrix,
	  cis_genotypes = cis_genotypes,
	  cis_effects_index = cis_effects_index,
	  data       = data
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

	BSFG_state = initialize_BSFG(BSFG_state, K_mats, chol_Ki_mats,verbose=verbose,ncores=ncores)

	# Initialize Eta
	data_model_state = data_model(Y,data_model_parameters,BSFG_state)
	BSFG_state$current_state[names(data_model_state$state)] = data_model_state$state

	# ----------------------- #
	# -Initialize Posterior-- #
	# ----------------------- #
	Posterior = list(
	  posteriorSample_params = unique(c(posteriorSample_params,data_model_state$posteriorSample_params)),
	  posteriorMean_params = unique(c(posteriorMean_params,data_model_state$posteriorMean_params)),
	  total_samples = 0,
	  folder = run_parameters$Posterior_folder,
	  files = c()
	)
	Posterior = reset_Posterior(Posterior,BSFG_state)
	BSFG_state$Posterior = Posterior


	return(BSFG_state)
}


initialize_BSFG = function(BSFG_state,...){
  UseMethod("initialize_BSFG",BSFG_state)
}

#' Run BSFG Gibbs sampler
#'
#' Run MCMC chain for a specified number of iterations
#'
#' @param BSFG_state BSFG_state object of current chain
#' @param n_samples Number of iterations to add to the chain (not number of posterior samples to draw.
#'     This is determined by n_samples / thin)
#' @param grainSize Minimum size of sub-problems for dividing among processes. Sent to RcppParallel
sample_BSFG = function(BSFG_state,n_samples,grainSize = 1,...) {
  data_matrices  = BSFG_state$data_matrices
  priors         = BSFG_state$priors
  run_parameters = BSFG_state$run_parameters
  run_variables  = BSFG_state$run_variables

  # ----------------------------------------------- #
  # -----------Reset Global Random Number Stream--- #
  # ----------------------------------------------- #
  do.call("RNGkind",as.list(BSFG_state$RNG$RNGkind))  ## must be first!
  assign(".Random.seed", BSFG_state$RNG$Random.seed, .GlobalEnv)

  # ----------------------------------------------- #
  # ----------------Set up run--------------------- #
  # ----------------------------------------------- #
  save_freq    = run_parameters$save_freq
  burn         = run_parameters$burn
  thin         = run_parameters$thin
  start_i      = BSFG_state$current_state$nrun

  # ----------------------------------------------- #
  # ---Extend posterior matrices for new samples--- #
  # ----------------------------------------------- #

  sp = (start_i + n_samples - burn)/thin - BSFG_state$Posterior$total_samples
  BSFG_state$Posterior = expand_Posterior(BSFG_state$Posterior,max(0,sp))

  # ----------------------------------------------- #
  # --------------start gibbs sampling------------- #
  # ----------------------------------------------- #

  start_time = Sys.time()
  for(i in start_i+(1:n_samples)){
    BSFG_state$current_state$nrun = i

    # ----- Sample Lambda ---------------- #
    BSFG_state$current_state = sample_Lambda_B(BSFG_state,grainSize = grainSize,...)

    # BSFG_state$current_state$Lambda[1,2:min(5,BSFG_state$current_state$k)] = 1  # TEMPORARY!!!

    # ----- Sample other factor model parameters  ---------------- #
    BSFG_state$current_state = sample_latent_traits(BSFG_state,grainSize = grainSize,...)

    # -----Sample Lambda_prec ------------- #
    BSFG_state$current_state = sample_Lambda_prec(BSFG_state)

    # -----Sample prec_B ------------- #
    BSFG_state$current_state = sample_prec_B_ARD(BSFG_state)

    # ----- sample Eta ----- #
    data_model_state = run_parameters$data_model(data_matrices$Y,run_parameters$data_model_parameters,BSFG_state)$state
    BSFG_state$current_state[names(data_model_state)] = data_model_state

    # -- adapt number of factors to samples ---#
    if(i > 200 && runif(1) < with(BSFG_state$run_parameters,1/exp(b0 + b1*i))){  # adapt with decreasing probability per iteration
      BSFG_state$current_state = update_k(BSFG_state)
    }

    # -- save sampled values (after thinning) -- #
    if( (i-burn) %% thin == 0 && i > burn) {
      BSFG_state$Posterior = save_posterior_sample(BSFG_state)
    }
  }
  end_time = Sys.time()
  print(end_time - start_time)
  BSFG_state$current_state$total_time = BSFG_state$current_state$total_time + end_time - start_time

  # ----------------------------------------------- #
  # ------------Save state for restart------------- #
  # ----------------------------------------------- #

  current_state = BSFG_state$current_state
  save(current_state,file='current_state.RData')

  BSFG_state$RNG = list(
    Random.seed = .Random.seed,
    RNGkind = RNGkind()
  )
  return(BSFG_state)
}

sample_Lambda_B = function(BSFG_state,...){
  UseMethod("sample_Lambda_B",BSFG_state)
}

sample_latent_traits = function(BSFG_state,...){
  UseMethod("sample_latent_traits",BSFG_state)
}

sample_Lambda_prec = function(BSFG_state) {
  priors         = BSFG_state$priors
  run_variables  = BSFG_state$run_variables
  current_state  = BSFG_state$current_state

  current_state = with(c(priors,run_variables),within(current_state,{
		Lambda2 = Lambda^2
		Lambda_prec = matrix(rgamma(p*k,shape = (Lambda_df + 1)/2,rate = (Lambda_df + sweep(Lambda2,2,tauh,'*'))/2),nr = p,nc = k)

		# # trait one is special?
		# Lambda_prec[1,] = 1e-10

	 # # -----Sample delta, update tauh------ #
		shapes = c(delta_1_shape + 0.5*p*k,
		           delta_2_shape + 0.5*p*((k-1):1))
		times = 100
		randg_draws = matrix(rgamma(times*k,shape = shapes,rate = 1),nr=times,byrow=T)
		delta[] = sample_delta_c( delta,tauh,Lambda_prec,delta_1_shape,delta_1_rate,delta_2_shape,delta_2_rate,randg_draws,Lambda2)
		tauh[]  = matrix(cumprod(delta),nrow=1)

	 # # -----Update Plam-------------------- #
		Plam[] = sweep(Lambda_prec,2,tauh,'*')
  }))
  return(current_state)
}

sample_prec_B = function(BSFG_state){
  priors         = BSFG_state$priors
  run_variables  = BSFG_state$run_variables
  current_state  = BSFG_state$current_state

  current_state = with(c(priors,run_variables),within(current_state,{
    if(b > 1) {
      if(b_F == b-1) {  # assume that X_F == X, want same tau for both
        B2 = cbind(B[-1,,drop=FALSE],B_F)^2
      } else{
        B2 = B[-1,,drop=FALSE]^2
      }
      tau_B[1,-1] = rgamma(b-1, shape = fixed_prec_shape + ncol(B2)/2, rate = fixed_prec_rate + rowSums(B2)/2)
      prec_B = matrix(tau_B,nrow = b, ncol = p)
    }
    if(b_F > 1){
      if(b_F == (b-1)) {
        tau_B_F[1,] = tau_B[1,-1]
      } else{
        B_F2 = B_F^2
        tau_B_F[1,] = rgamma(b_F, shape = fixed_prec_shape + ncol(B_F2)/2, rate = fixed_prec_rate + rowSums(B_F2)/2)
      }
      prec_B_F = matrix(tau_B_F,nrow = b_F, ncol = k)
    }
  }))
  return(current_state)
}

sample_prec_B_ARD = function(BSFG_state){
  priors         = BSFG_state$priors
  run_variables  = BSFG_state$run_variables
  current_state  = BSFG_state$current_state
  current_state = with(c(priors,run_variables),within(current_state,{
    if(b > 1) {
      B2 = B^2
      tau_B[1,-1] = rgamma(b-1, shape = fixed_prec_shape + ncol(B2)/2, rate = fixed_prec_rate + rowSums((B2 * prec_B/c(tau_B))[-1,,drop=FALSE])/2)
      prec_B[-1,] = matrix(rgamma((b-1)*p,shape = (B_df + 1)/2,rate = (B_df + B2[-1,,drop=FALSE]*tau_B[-1])/2),nr = (b-1),nc = p)
      prec_B[-1,] = prec_B[-1,]*tau_B[-1]
    }
    if(b_F > 0) {
      B_F2 = B_F^2
      tau_B_F[1,] = rgamma(b_F, shape = fixed_prec_shape + ncol(B_F2)/2, rate = fixed_prec_rate + rowSums(B_F2 * prec_B_F/c(tau_B_F))/2)
      prec_B_F[] = matrix(rgamma(b_F*k,shape = (B_F_df + 1)/2,rate = (B_F_df + B_F2*c(tau_B_F))/2),nr = b_F,nc = k)
      prec_B_F[] = prec_B_F*c(tau_B_F)
    }
  }))
  return(current_state)
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
plot.BSFG_state = function(BSFG_state,file = 'diagnostics_plots.pdf'){
  pdf(file)
  if(BSFG_state$run_parameters$simulation){
    plot_diagnostics_simulation(BSFG_state)
  } else{
    plot_diagnostics(BSFG_state)
  }
  dev.off()
}
