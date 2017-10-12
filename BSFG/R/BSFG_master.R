  #' Set BSFG run parameters
#'
#' Function to create run_parameters list for initializing BSFG model
#'
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
#' @param k_init initial number of factors
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
#' @param delta_iteractions_factor Number of times to iterate through sample_delta per iteration of the other parameters
#' @param impute_missing If FALSE, missing data will be ignored in the sampler. If TRUE, missing data will be imputed during sampling.
#' @seealso \code{\link{BSFG_init}}, \code{\link{sample_BSFG}}, \code{\link{print.BSFG_state}}
#'
BSFG_control = function(
                        Posterior_folder = "Posterior",
                        simulation = c(F,T),scale_Y = c(T,F),
                        b0 = 1, b1 = 0.0005, epsilon = 1e-1, prop = 1.00,
                        k_init = 20, h2_divisions = 100, h2_step_size = NULL,
                        drop0_tol = 1e-14, K_eigen_tol = 1e-10,
                        burn = 100,thin = 2,
                        delta_iteractions_factor = 100,
                        impute_missing = FALSE,
                        ...
                        ) {

  formals_named = formals()
  formals_named = formals_named[names(formals_named) != '...']
  all_args = lapply(formals_named,function(x) eval(x)[1])
  passed_args = lapply(as.list(match.call())[-1],eval)
  if(any(names(passed_args) %in% names(formals_named) == F)){
    unused_names = names(passed_args)[names(passed_args) %in% names(formals_named) == F]
    warning(sprintf('No argument(s) named %s',paste(unused_names,sep=', ')))
  }
  all_args[names(passed_args)] = passed_args
  return(all_args)
}



#' Set BSFG priors
#'
#' Function to create list of priors for BSFG model.
#'
#' Default values are provided, but any can be replaced. Note: \code{h2_priors_resids} and
#'     \code{h2_priors_factors} can be set after calling \code{BSFG_init} and before \code{sample_BSFG}
#'     if that is easier.
#'
#' @param fixed_var List of parameters of inverse gamma distribution for fixed effects, specifically:
#'     \code{V} and \code{nu}, give shape = \code{nu/2} and scale = \code{nu*V/2}, so mean = \code{\frac{V*nu}{nu-2}}
#'     Will be applied to fixed effects of residuals (\code{B}) and factors (\code{B_F}).
#' @param fixed_resid_var If provided, overides \code{fixed_var} for fixed effects of residuals (\code{B}).
#' @param fixed_factors_var If provided, overides \code{fixed_var} for fixed effects of factors (\code{B_F}).
#' @param tot_Y_var List of parameters of inverse gamma distribution for residual variances. See \code{fixed_var}.
#' @param tot_F_var List of parameters of inverse gamma distribution for factor variances. See \code{fixed_var}.
#'     This parameter provides the parameter extension of Ghosh and Dunson (2009), but is removed
#'     from all Factor parameters before they are saved in Posterior
#' @param delta_1 List of parameters of inverse gamma distribution for \code{delta_1}. Specifically:
#'     \code{shape} and \code{rate}. This parameter is the column-shrinkage of the first factor.
#' @param delta_2 List of parameters of inverse gamma distribution for \code{delta_2 \dots delta_k}.
#'     Specifically: \code{shape} and \code{rate}.
#'     This is provides the additional column-shrinkage of higher-order columns.
#' @param Lambda_df Degrees of freedom of individual parameter shrinkage of Lambda from implied
#'     t-distribution
#' @param B_df Degrees of freedom of individual parameter shrinkage of B from implied
#'     t-distribution
#' @param B_F_df Degrees of freedom of individual parameter shrinkage of B_F from implied
#'     t-distribution
#' @param h2_priors_resids_fun function for that returns prior probability for a given value of h2
#'     for each random effect. Should take two argument - a vector \code{h2} values for each random effect,
#'     and \code{n} - the number of discrete levels of the prior.
#'     Alternatively, can be a scalar or vector of (relative) prior values for each value of the
#'     discrete prior.
#' @param h2_priors_factors_fun see \code{h2_priors_resids_fun}. Same, but for the h2s of the factors.
#'
#' @return a list with each of the prior components specified above.
#' @export
#'
#' @examples
BSFG_priors = function(
                        fixed_var = list(V = 1,     nu = 3),
                        fixed_resid_var = NULL,
                        fixed_factors_var = NULL,
                        QTL_resid_var = NULL,
                        QTL_factors_var = NULL,
                        tot_Y_var = list(V = 0.5,   nu = 3),
                        tot_F_var = list(V = 18/20, nu = 20),
                        h2_priors_resids_fun = function(h2s, n) 1,
                        h2_priors_factors_fun = function(h2s, n) 1,
                        Lambda_prior = list(
                                            sampler = sample_Lambda_prec_ARD,
                                            Lambda_df = 3,
                                            delta_1   = list(shape = 2.1,  rate = 1/20),
                                            delta_2   = list(shape = 3, rate = 1)
                          ),
                        B_prior = list(
                                        sampler = sample_B_prec_ARD,
                                        B_df      = 3,
                                        B_F_df    = 3
                                        )

                    ) {
  all_args = lapply(formals(),function(x) eval(x))
  passed_args = lapply(as.list(match.call())[-1],function(x) eval(x))
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
#' The model is specified as:
#'
#' y_i = g(eta_i)
#' Eta = rbind(eta_1,...,eta_n) = X*B + F*t(Lambda) + Z*U_R + E_R
#' F = X_F * B_F + Z*U_F + E_F
#'
#' For sampling, we reparameterize as:
#'
#' Qt*Eta = Qt*X*B + Qt*F*t(Lambda) + Qt*Z*U_R + Qt*E_R
#' Qt*F = Qt*X_F * B_F + Qt*Z*U_F + Qt*E_F
#'
#' We sample the quantities Qt*Eta, Qt*F, B, Lambda, U_R, U_F. We then back-calculate Eta and F.
#'
#' @param Y either a) a n x p matrix of data (n individuals x p traits), or b) a list describing
#'     the observation_model, data, and associated parameters. This list should contain:
#'     i) \code{observation_model}: a function modeled after \code{missing_observation_model}
#'         (see code by typing this into the console) that draws posterior samples of Eta conditional on
#'         the observations (Y) and the current state of the BSFG model (current_state).
#'         The function should have the same form as \link{missing_data},
#'         and must return Eta even if current_state is NULL.
#'     ii) \code{observations}: a data.frame containing the observaition-level data and associated covariates.
#'         This must include a column \code{ID} that is also present in \code{data}
#'     iii) any other parameters necessary for \code{observation_model}
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
#'     the fixed effects for the latent factors will be the same as for Eta, except that the columns will be centered.
#'     If a formula is provided (no random effects allowed), the model will be applied to
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
#' @param QTL_resid a design matrix. Intended for marker genotypes. Corresponding coefficients
#'     will be modeled with the \code{QTL_resid_var} priors. Applies to factor residuals \eqn{Y - F\Lambda^T}
#' @param QTL_factors a design matrix. Intended for marker genotypes. Corresponding coefficients
#'     will be modeled with the \code{QTL_factor_var} priors. Applies to factors \eqn{F}.
#' @param cis_genotypes a list of design matrices of length \code{p} (ie number of columns of \code{Eta})
#'     This is used to specify trait-specific fixed effects, such a cis-genotypes
#' @param posteriorSample_params A character vector giving names of parameters to save all posterior samples
#' @param posteriorMean_params A character vector giving names of parameters to save only the posterior mean.
#' @param ncores for \code{general_BSFG}, number of cores to use during initialization.
#' @param setup optional - a list of known values for Lambda (error_factor_lambda), h2, factor_h2s
#' @param verbose should progress in initialization be reported?
#' @param Sigma_Choleskys Pre-calculated matrices from \code{BSFG_state$run_variables$Sigma_Choleskys}
#'     can be provide directly. For general_BSFG sampler.
#' @param randomEffect_C_Choleskys Pre-calculated matrices from \code{BSFG_state$run_variables$randomEffect_C_Choleskys}
#'     can be provide directly. For general_BSFG sampler.
#' @param invert_aI_bZKZ Pre-calculated matrices from \code{BSFG_state$run_variables$invert_aI_bZKZ}
#'     can be provide directly. For fast_BSFG sampler.
#' @param invert_aZZt_Kinv Pre-calculated matrices from \code{BSFG_state$run_variables$invert_aZZt_Kinv}
#'     can be provide directly. For fast_BSFG sampler.#'
#' @return An object of class BSFG_state with components:
#' @return current_state: a list of parameters in the current iteration of the sampler
#' @return Posterior: a list of arrays of posterior samples
#' @return RNG current state of R's Random number generator (for re-starting chaings)
#' @return traitnames: vector of trait names (from colnames of Y)
#' @return run_parameters, run_variables, data_matrices, priors, simulation: input data and
#'   parameters
#' @seealso \code{\link{BSFG_control}}, \code{\link{sample_BSFG}}, \code{\link{print.BSFG_state}}, \code{\link{plot.BSFG_state}}#'
BSFG_init = function(Y, model, data, factor_model_fixed = NULL, priors = BSFG_priors(), run_parameters = BSFG_control(), K_mats = NULL, K_inv_mats = NULL,
                     QTL_resid = NULL, QTL_factors = NULL, cis_genotypes = NULL,
                     Sigma_Choleskys = NULL, randomEffect_C_Choleskys = NULL,
                     invert_aI_bZKZ = NULL, invert_aZZt_Kinv = NULL,
                     posteriorSample_params = c('Lambda','U_F','F','delta','tot_F_prec','F_h2','tot_Eta_prec','resid_h2', 'B', 'B_F','U_R','cis_effects'),
                     posteriorMean_params = c(),
                     ncores = detectCores(),setup = NULL,verbose=T) {

  # ----------------------------- #
  # ---- build model matrices --- #
  # ----------------------------- #

  if(is(Y,'list')){
    if(!'observation_model' %in% names(Y)) stop('observation_model not specified in Y')
    observation_model = Y$observation_model
    observation_model_parameters = Y[names(Y) != 'observation_model']
  } else{
    if(!is(Y,'matrix'))	Y = as.matrix(Y)
    observation_model = missing_data_model
    observation_model_parameters = list(
      Y = Y,
      scale_Y = run_parameters$scale_Y
    )
  }

  # initialize observation_model
  observation_model_parameters$observation_setup = observation_model(observation_model_parameters,list(data_matrices = list(data = data)))
  n = nrow(data)
  p = observation_model_parameters$observation_setup$p
  traitnames = observation_model_parameters$observation_setup$traitnames
  if(is.null(traitnames)) traitnames = paste('trait',1:p,sep='_')
  if(is.null(observation_model_parameters$observation_setup$Y_missing)) {
    observation_model_parameters$observation_setup$Y_missing = matrix(0,n,p)
  }
  if(!is(observation_model_parameters$observation_setup$Y_missing,'lgTMatrix')){
    observation_model_parameters$observation_setup$Y_missing = as(observation_model_parameters$observation_setup$Y_missing,'lgTMatrix')
  }
  Y_missing = observation_model_parameters$observation_setup$Y_missing

	# if factor_model_fixed not specified, use fixed effects from model for both
	if(is.null(factor_model_fixed)) {
	  factor_model_fixed = nobars(model)
	  same_fixed_model = TRUE  # keep track if the factors and residuals have same fixed effects
	} else{
	  same_fixed_model = FALSE
	}

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
	resid_intercept = FALSE     # is there an intercept? If so, don't penalize coefficient.
	QTL_columns_resid = NULL    # indexes of columns of X for QTLs. They get a different penalty

	X = model.matrix(nobars(model),data)
	linear_combos = caret::findLinearCombos(X)
	if(!is.null(linear_combos$remove)) {
	  cat(sprintf('dropping column(s) %s to make X_resid full rank\n',paste(linear_combos$remove,sep=',')))
	  X = X[,-linear_combos$remove]
	}
	if(ncol(X) > 0 && all(X[,1] == 1)) {
	  resid_intercept = TRUE
	}
	b = ncol(X)
	X = as(X,'dgCMatrix')
	if(any(is.na(X))) stop('Missing values in X_resid')

	QTL_resid_Z = QTL_resid_X = NULL
	if(!is.null(QTL_resid)){
	  if(class(QTL_resid) == 'list'){
	    QTL_terms = mkReTrms(findbars(QTL_resid$model),data,drop.unused.levels = FALSE)
	    QTL_resid_Z = as(t(QTL_terms$Zt),'dgCMatrix')
	    if(!all(colnames(QTL_factors_Z) %in% rownames(QTL_resid$X))) stop(sprintf('Missing %s from QTL_resid$X',names(QTL_terms$cnms)[1]))
	    QTL_resid_X = QTL_resid$X[colnames(QTL_factors_Z),]
	  } else{
	    if(nrow(QTL_resid) != nrow(data)) stop(sprintf('QTL_resid has wrong number of rows. Should be %d',nrow(data)))
	    QTL_resid_Z = as(diag(1,nrow(data)),'dgCMatrix')
	    if(is.data.frame(QTL_resid)) QTL_resid = as.matrix(QTL_resid)
	    QTL_resid_X = QTL_resid
	  }
	  same_fixed_model = FALSE
	  if(any(is.na(QTL_resid_X))) stop('Missing values in QTL_resid_X')
	  QTL_columns_resid = b+1:ncol(QTL_resid_X)
	  X = cbind(X,QTL_resid_Z %*% QTL_resid_X)
	  b = ncol(X)
	}

  # for F
	QTL_columns_factors = NULL
	  # note columns are centered, potentially resulting in zero-variance columns
	X_F = model.matrix(factor_model_fixed,data)
	linear_combos = caret::findLinearCombos(X_F)
	if(!is.null(linear_combos$remove)) {
	  cat(sprintf('dropping column(s) %s to make X_factor full rank\n',paste(linear_combos$remove,sep=',')))
	  X_F = X_F[,-linear_combos$remove]
	}
	b_F = ncol(X_F)
	X_F_zero_variance = apply(X_F,2,var) == 0
	X_F = as(X_F,'dgCMatrix')
	if(any(is.na(X_F))) stop('Missing values in X_F')


	QTL_factors_Z = QTL_factors_X = NULL
	if(!is.null(QTL_factors)){
	  if(class(QTL_factors) == 'list'){
	    QTL_terms = mkReTrms(findbars(QTL_factors$model),data,drop.unused.levels = FALSE)
	    QTL_factors_Z = as(t(QTL_terms$Zt),'dgCMatrix')
	    if(!all(colnames(QTL_factors_Z) %in% rownames(QTL_factors$X))) stop(sprintf('Missing %s from QTL_factors$X',names(QTL_terms$cnms)[1]))
	    QTL_factors_X = QTL_factors$X[colnames(QTL_factors_Z),]
	  } else{
	    if(nrow(QTL_factors) != nrow(data)) stop(sprintf('QTL_factors has wrong number of rows. Should be %d',nrow(data)))
	    QTL_factors_Z = as(diag(1,nrow(data)),'dgCMatrix')
	    if(is.data.frame(QTL_factors)) QTL_factors = as.matrix(QTL_factors)
	    QTL_factors_X = QTL_factors
	  }
	  same_fixed_model = FALSE
	  if(any(is.na(QTL_factors_X))) stop('Missing values in QTL_factors_X')
	  QTL_columns_factors = b_F+1:ncol(QTL_factors_X)
	  X_F = cbind(X_F,QTL_factors_Z %*% QTL_factors_X)
	  b_F = ncol(X_F)
	  X_F_zero_variance = apply(X_F,2,var) == 0
	  X_F = as(X_F,'dgCMatrix')
	}


	# -------- cis genotypes ---------- #
	if(is.null(cis_genotypes)){
	  n_cis_effects = NULL
	  cis_effects_index = NULL
	} else{
	  n_cis_effects = sapply(cis_genotypes,ncol)
	  cis_effects_index = c(0,cumsum(n_cis_effects))+1
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
	n_RE = length(RE_names)

	# find RE indices
	RE_lengths = sapply(Z_matrices,ncol)
	RE_starts = cumsum(c(0,RE_lengths)[1:n_RE])
	names(RE_starts) = RE_names
	RE_indices = lapply(RE_names,function(re) RE_starts[re] + 1:RE_lengths[re])
	names(RE_indices) = RE_names

	  # function to ensure that covariance matrices are sparse and symmetric
	fix_K = function(x) forceSymmetric(drop0(x,tol = run_parameters$drop0_tol))

	# construct K matrices for each random effect
	    # if K is PSD, find K = LDLt and set K* = D and Z* of ZL
	ldl_ks = lapply(K_mats,function(K) {
	  res = LDLt_sparse(as(K,'dgCMatrix')) # actually calculates K = PtLDLtP
	  if(is.character(validObject(res$L,test=TRUE)[1])) {
	    res = LDLt_notSparse(as.matrix(K))  # sparse sometimes fails with non-PD K
	  }
	  res
	})

	RE_L_matrices = list()

	for(i in 1:n_RE){
	  re = RE_covs[i]
	  re_name = RE_names[i]
	  if(re %in% names(ldl_ks)) {
	    ldl_k = ldl_ks[[re]]
	    large_d = ldl_k$d > run_parameters$K_eigen_tol
	    r_eff = sum(large_d)
	    # if need to use reduced rank model, then use D of K in place of K and merge L into Z
	    # otherwise, use original K, set L = Diagonal(1,r)
	    if(r_eff < length(ldl_k$d)) {
	      K = Diagonal(r_eff,ldl_k$d[large_d])
	      L = t(ldl_k$P) %*% ldl_k$L[,large_d]
	    } else{
	      K = K_mats[[re]]
	      L = Diagonal(nrow(K),1)
	    }
	  } else if(re %in% names(K_inv_mats)){
	    K_mats[[re]] = solve(K_inv_mats[[re]])
	    rownames(K_mats[[re]]) = rownames(K_inv_mats[[re]])
	    K = K_mats[[re]]
	    L = Diagonal(nrow(K),1)
	  } else{
	    K_mats[[re]] = Diagonal(ncol(Z_matrices[[re_name]]))
	    rownames(K_mats[[re]]) = levels(data[[re]])
	    K = K_mats[[re]]
	    L = Diagonal(nrow(K),1)
	  }
	  K_mats[[re_name]] = fix_K(K)
	  RE_L_matrices[[re_name]] = L
	  Z_matrices[[re_name]] = Z_matrices[[re_name]] %*% L
	}
	K_mats = K_mats[RE_names]
	# Fix Z_matrices based on PSD K's
	Z = do.call(cbind,Z_matrices[RE_names])
	Z = as(Z,'dgCMatrix')
	# The following matrix is used to transform random effects back to the original space had we sampled from the original (PSD) K.
	if(length(RE_names) > 1) {
	  RE_L = do.call(bdiag,RE_L_matrices[RE_names])
	} else{
	  RE_L = RE_L_matrices[[1]]
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


	# ------------------------------------ #
	# ----Precalculate ZKZts, chol_Ks ---- #
	# ------------------------------------ #

	# first, identify sets of traits with same pattern of missingness

# ideally, want to be able to restrict the number of sets. Should be possible to merge sets of columngs together.
	if(!run_parameters$impute_missing) {
    Y_col_obs = lapply(1:ncol(Y_missing),function(x) {
      obs = which(!Y_missing[,x],useNames=F)
      names(obs) = NULL
      obs
    })
    unique_Y_col_obs = unique(c(list(1:nrow(Y_missing)),Y_col_obs))
    unique_Y_col_obs_str = lapply(unique_Y_col_obs,paste,collapse='')
    Y_col_obs_index = sapply(Y_col_obs,function(x) which(unique_Y_col_obs_str == paste(x,collapse='')))

    Missing_data_map = lapply(seq_along(unique_Y_col_obs),function(i) {
      x = unique_Y_col_obs[[i]]
      if(length(x) == 0){
        stop('Columns of Y have 100% missing data! Please drop these columns.')
      }
      return(list(
        Y_obs = x,
        Y_cols = which(Y_col_obs_index == i)
      ))
    })
	} else{
	  Missing_data_map = list(list(
	    Y_obs = 1:n,
	    Y_cols = 1:p
	  ))
	}

  # now, for each set of columns, pre-calculate a set of matrices, etc
  # do calculations in several chunks
  group_size = 2*parallel::detectCores()
  n_groups = ceiling(ncol(h2s_matrix)/group_size)
  col_groups = tapply(1:ncol(h2s_matrix),gl(n_groups,group_size,ncol(h2s_matrix)),function(x) x)

  if(verbose) {
    print(sprintf("Pre-calculating random effect inverse matrices for %d groups of traits and %d sets of random effect weights", length(Missing_data_map), ncol(h2s_matrix)))
    pb = txtProgressBar(min=0,max = length(Missing_data_map)*  length(col_groups) * 2,style=3)
  }

  Qt_list = list()
  QtX_list = list()
  QtXF_list = list()
  QtZ_list = list()
  Qt_cis_genotypes_list = list()
  randomEffect_C_Choleskys_list = list()
  Sigma_Choleskys_list = list()

  for(set in seq_along(Missing_data_map)){
    x = Missing_data_map[[set]]$Y_obs
    cols = Missing_data_map[[set]]$Y_cols

    if(n_RE == 1){
      ZKZt = Z[x,,drop=FALSE] %*% K_mats[[1]] %*% t(Z[x,,drop=FALSE])
      result = svd(ZKZt)
      Qt = as(drop0(as(t(result$u),'dgCMatrix'),tol = run_parameters$drop0_tol),'dgCMatrix')
    } else{
      Qt = as(diag(1,length(x)),'dgCMatrix')
    }
    QtZ_matrices_set = lapply(Z_matrices,function(z) Qt %*% z[x,,drop=FALSE])
    QtZ_set = do.call(cbind,QtZ_matrices_set[RE_names])
    QtZ_set = as(QtZ_set,'dgCMatrix')
    QtX_set = Qt %**% X[x,,drop=FALSE]
    QtXF_set = Qt %**% X_F[x,,drop=FALSE]  # This one is only needed for set==1
    Qt_cis_genotypes_set = lapply(cis_genotypes[cols],function(X) Qt %**% X[x,])

    Qt_list[[set]]   = Qt
    QtX_list[[set]]  = QtX_set
    QtXF_list[[set]] = QtXF_set
    QtZ_list[[set]]  = QtZ_set
    Qt_cis_genotypes_list[[set]] = Qt_cis_genotypes_set

    ZKZts_set = list()
    for(i in 1:n_RE){
      ZKZts_set[[i]] = as(forceSymmetric(drop0(QtZ_matrices_set[[i]] %*% K_mats[[i]] %*% t(QtZ_matrices_set[[i]]),tol = run_parameters$drop0_tol)),'dgCMatrix')
    }


    Sigma_Choleskys_c_list = list()
    for(i in 1:length(col_groups)){
      Sigma_Choleskys_c_list[[i]] = new(Sigma_Cholesky_database,ZKZts_set,h2s_matrix[,col_groups[[i]],drop=FALSE],run_parameters$drop0_tol,1)
      if(verbose) setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
    }
    Sigma_Choleskys_list[[set]] = do.call(c,lapply(1:length(col_groups),function(j) {
      Sigma_Choleskys_c = Sigma_Choleskys_c_list[[j]]
      lapply(1:length(col_groups[[j]]),function(i) {
        list(log_det = Sigma_Choleskys_c$get_log_det(i),
             chol_Sigma = drop0(Sigma_Choleskys_c$get_chol_Sigma(i),tol = run_parameters$drop0_tol)
            )
      })
    }))

    ZtZ_set = as(forceSymmetric(drop0(crossprod(Z[x,]),tol = run_parameters$drop0_tol)),'dgCMatrix')

    randomEffect_C_Choleskys_c_list = list()
    for(i in 1:length(col_groups)){
      randomEffect_C_Choleskys_c_list[[i]] = new(randomEffect_C_Cholesky_database,lapply(chol_Ki_mats,function(x) as(x,'dgCMatrix')),h2s_matrix[,col_groups[[i]],drop=FALSE],ZtZ_set,run_parameters$drop0_tol,1)
      if(verbose) setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
    }
    randomEffect_C_Choleskys_list[[set]] = do.call(c,lapply(1:length(col_groups),function(j) {
      randomEffect_C_Choleskys_c = randomEffect_C_Choleskys_c_list[[j]]
      lapply(1:length(col_groups[[j]]),function(i) {
        list(chol_C     = randomEffect_C_Choleskys_c$get_chol_Ci(i),
             chol_K_inv = randomEffect_C_Choleskys_c$get_chol_K_inv_i(i)
        )
      })
    }))
  }
  if(verbose) close(pb)

  if(!is.null(QTL_factors_Z)){
    Qt1_QTL_Factors_Z = Qt_list[[1]] %*% QTL_factors_Z  # since QTL_factors_Z only used with complete data, only calculate this for set 1.
  } else{
    Qt1_QTL_Factors_Z = NULL
  }


  run_variables = list(
    p      = p,
    n      = n,
    r_RE   = r_RE,
    RE_names = RE_names,
    b      = b,
    b_F    = b_F,
    resid_intercept   = resid_intercept,
    same_fixed_model  = same_fixed_model,
    X_F_zero_variance = X_F_zero_variance,   # used to identify fixed effect coefficients that should be forced to zero
    n_cis_effects     = n_cis_effects,
    cis_effects_index = cis_effects_index,
    Qt_list    = Qt_list,
    QtX_list   = QtX_list,
    QtXF_list  = QtXF_list,
    QtZ_list   = QtZ_list,
    Qt_cis_genotypes_list = Qt_cis_genotypes_list,
    Qt1_QTL_Factors_Z = Qt1_QTL_Factors_Z,
    Missing_data_map      = Missing_data_map,
    Sigma_Choleskys_list          = Sigma_Choleskys_list,
    randomEffect_C_Choleskys_list = randomEffect_C_Choleskys_list
  )

	data_matrices = list(
	  X          = X,
	  X_F        = X_F,
	  Z_matrices = Z_matrices,
	  Z          = Z,
	  RE_L       = RE_L,  # matrix necessary to back-transform U_F and U_R (RE_L*U_F and RE_L*U_R) to get original random effects
	  RE_indices = RE_indices,
	  h2s_matrix = h2s_matrix,
	  cis_genotypes = cis_genotypes,
	  QTL_columns_resid = QTL_columns_resid,
	  QTL_columns_factors = QTL_columns_factors,
	  QTL_resid_Z = QTL_resid_Z,
	  QTL_resid_X = QTL_resid_X,
	  QTL_factors_Z = QTL_factors_Z,
	  QTL_factors_X = QTL_factors_X,
	  data       = data
	)

	run_parameters$observation_model = observation_model
	run_parameters$observation_model_parameters = observation_model_parameters
	run_parameters$traitnames = traitnames


	# ----------------------------- #
	# ----- re-formulate priors --- #
	# ----------------------------- #
	if(any(sapply(priors, function(x) {try({return(exists(x$nu) && x$nu <= 2)},silent=T);return(FALSE)}))) stop('priors nu must be > 2')
	  # fixed effects
	if(is.null(priors$fixed_resid_var)) priors$fixed_resid_var = priors$fixed_var
	if(is.null(priors$fixed_factors_var)) priors$fixed_factors_var = priors$fixed_var
	if(length(priors$fixed_resid_var$V == 1)) {
	  priors$fixed_resid_var$V = rep(priors$fixed_resid_var$V,b)
	  priors$fixed_resid_var$nu = rep(priors$fixed_resid_var$nu,b)
	}

	if(length(priors$fixed_factors_var$V) == 0){
	  priors$fixed_factors_var = priors$fixed_resid_var
	} else if(length(priors$fixed_factors_var$V) == 1){
	  priors$fixed_factors_var$V = rep(priors$fixed_factors_var$V,b)
	  priors$fixed_factors_var$nu = rep(priors$fixed_factors_var$nu,b)
	}

	if(length(QTL_columns_resid)>0){
	  if(length(priors$QTL_resid_var$V) == 1) {
	    priors$QTL_resid_var$V = rep(priors$QTL_resid_var$V,length(QTL_columns_resid))
	    priors$QTL_resid_var$nu = rep(priors$QTL_resid_var$nu,length(QTL_columns_resid))
	  }
	  priors$fixed_resid_var$V[QTL_columns_resid] = priors$QTL_resid_var$V
	  priors$fixed_resid_var$nu[QTL_columns_resid] = priors$QTL_resid_var$nu
	}
	if(length(QTL_columns_factors)>0){
	  if(length(priors$QTL_factors_var$V) == 1) {
	    priors$QTL_factors_var$V = rep(priors$QTL_factors_var$V,length(QTL_columns_factors))
	    priors$QTL_factors_var$nu = rep(priors$QTL_factors_var$nu,length(QTL_columns_factors))
	  }
	  priors$fixed_factors_var$V[QTL_columns_factors] = priors$QTL_factors_var$V
	  priors$fixed_factors_var$nu[QTL_columns_factors] = priors$QTL_factors_var$nu
	}

	priors$fixed_resid_prec_rate   = with(priors$fixed_resid_var,V * nu)
	priors$fixed_resid_prec_shape  = with(priors$fixed_resid_var,nu - 1)
	priors$fixed_factors_prec_rate   = with(priors$fixed_factors_var,V * nu)
	priors$fixed_factors_prec_shape  = with(priors$fixed_factors_var,nu - 1)

	  # total precision
	if(length(priors$tot_Y_var$V == 1)) {
	  priors$tot_Y_var$V = rep(priors$tot_Y_var$V,p)
	  priors$tot_Y_varnu = rep(priors$tot_Y_var$nu,p)
	}
	priors$tot_Eta_prec_rate   = with(priors$tot_Y_var,V * nu)
	priors$tot_Eta_prec_shape  = with(priors$tot_Y_var,nu - 1)
	priors$tot_F_prec_rate     = with(priors$tot_F_var,V * nu)
	priors$tot_F_prec_shape    = with(priors$tot_F_var,nu - 1)
	#   # delta: column shrinkage of Lambda
	# priors$delta_1_rate    = priors$delta_1$rate
	# priors$delta_1_shape   = priors$delta_1$shape
	# priors$delta_2_rate    = priors$delta_2$rate
	# priors$delta_2_shape   = priors$delta_2$shape

	# h2_priors_resids
	if(exists('h2_priors_resids',priors)) {
	  if(length(priors$h2_priors_resids) == 1) priors$h2_priors_resids = rep(priors$h2_priors_resids,ncol(h2s_matrix))
	  if(!length(priors$h2_priors_resids) == ncol(h2s_matrix)) stop('wrong length of priors$h2_priors_resids')
	} else{
	  if(!is(priors$h2_priors_resids_fun,'function')) stop('need to provide a priors$h2_priors_resids_fun() to specify discrete h2 prior for resids')
	  priors$h2_priors_resids = apply(h2s_matrix,2,priors$h2_priors_resids_fun,n = ncol(h2s_matrix))
	}
	priors$h2_priors_resids = priors$h2_priors_resids/sum(priors$h2_priors_resids)
	# h2_priors_factors
	if(exists('h2_priors_factors',priors)) {
	  if(length(priors$h2_priors_factors) == 1) priors$h2_priors_factors = rep(priors$h2_priors_factors,ncol(h2s_matrix))
	  if(!length(priors$h2_priors_factors) == ncol(h2s_matrix)) stop('wrong length of priors$h2_priors_factors')
	} else{
	  if(!is(priors$h2_priors_factors_fun,'function')) stop('need to provide a priors$h2_priors_factors_fun() to specify discrete h2 prior for factors')
	  priors$h2_priors_factors = apply(h2s_matrix,2,priors$h2_priors_factors_fun,n = ncol(h2s_matrix))
	}
	priors$h2_priors_factors = priors$h2_priors_factors/sum(priors$h2_priors_factors)

	# ----------------------------- #
	# -- create BSFG_state object - #
	# ----------------------------- #

	RNG = list(
	  Random.seed = .Random.seed,
	  RNGkind = RNGkind()
	)

	BSFG_state = list(
	  data_matrices  = data_matrices,
	  priors         = priors,
	  run_parameters = run_parameters,
	  run_variables  = run_variables,
	  RNG            = RNG,
	  setup          = setup
	)
	class(BSFG_state) = append('BSFG_state',class(BSFG_state))

	# ----------------------------- #
	# --- Initialize BSFG_state --- #
	# ----------------------------- #

	BSFG_state$Posterior = list(
	  posteriorSample_params = posteriorSample_params,
	  posteriorMean_params = posteriorMean_params,
	  total_samples = 0,
	  folder = run_parameters$Posterior_folder,
	  files = c()
	)

	BSFG_state = initialize_variables(BSFG_state)

	Posterior = reset_Posterior(BSFG_state$Posterior,BSFG_state)
	BSFG_state$Posterior = Posterior


	return(BSFG_state)
}

initialize_variables = function(BSFG_state,...){
  run_parameters = BSFG_state$run_parameters
  run_variables = BSFG_state$run_variables
  data_matrices = BSFG_state$data_matrices
  priors = BSFG_state$priors

  BSFG_state$current_state = with(c(run_parameters,run_variables,data_matrices,priors),{

    # Factors loadings:
    #  initial number of factors
    k = k_init

    Plam = matrix(1,p,k)

    # Lambda - factor loadings
    #   Prior: Normal distribution for each element.
    #       mu = 0
    #       sd = sqrt(1/Plam)
    Lambda = matrix(rnorm(p*k,0,sqrt(1/Plam)),nr = p,nc = k)
    rownames(Lambda) = traitnames

    # residuals
    # p-vector of factor precisions. Note - this is a 'redundant' parameter designed to give the Gibbs sampler more flexibility
    #  Prior: Gamma distribution for each element
    #       shape = tot_Eta_prec_shape
    #       rate = tot_Eta_prec_rate
    tot_Eta_prec = matrix(rgamma(p,shape = tot_Eta_prec_shape,rate = tot_Eta_prec_rate),nrow = 1)
    colnames(tot_Eta_prec) = traitnames

    # p-vector of factor precisions. Note - this is a 'redundant' parameter designed to give the Gibbs sampler more flexibility
    #  Prior: Gamma distribution for each element
    #       shape = tot_F_prec_shape
    #       rate = tot_F_prec_rate
    tot_F_prec = matrix(1,nrow=1,ncol=k)
    #with(priors,matrix(rgamma(k,shape = tot_F_prec_shape,rate = tot_F_prec_rate),nrow=1))

    # Factor scores:

    # Resid discrete variances
    # p-matrix of n_RE x p with
    resid_h2_index = sample(1:ncol(h2s_matrix),p,replace=T)
    resid_h2 = h2s_matrix[,resid_h2_index,drop=FALSE]

    # Factor discrete variances
    # k-matrix of n_RE x k with
    F_h2_index = sample(1:ncol(h2s_matrix),k,replace=T)
    F_h2 = h2s_matrix[,F_h2_index,drop=FALSE]

    U_F = lapply(RE_names,function(effect){
      matrix(rnorm(r_RE[effect] * k, 0, sqrt(F_h2[effect,] / tot_F_prec)),ncol = k, byrow = T)
    })
    names(U_F) = RE_names

    U_R = do.call(rbind,lapply(RE_names,function(effect){
      matrix(rnorm(r_RE[effect] * p, 0, sqrt(resid_h2[effect,] / tot_Eta_prec)),ncol = p, byrow = T)
    }))
    colnames(U_R) = traitnames
    rownames(U_R) = colnames(Z)

    # Fixed effects
    B = 0*matrix(rnorm(b*p), ncol = p)
    colnames(B) = traitnames

    # Factor fixed effects
    B_F = 0*matrix(rnorm(b_F * k),b_F,k)

    # cis effects
    cis_effects = matrix(rnorm(cis_effects_index[length(cis_effects_index)]-1,0,1),nrow=1)

    XB = as.matrix(X %*% B)

    F = X_F %*% B_F + matrix(rnorm(n * k, 0, sqrt((1-colSums(F_h2)) / tot_F_prec)),ncol = k, byrow = T)
    for(effect in RE_names) {
      F = F + Z_matrices[[effect]] %*% U_F[[effect]]
    }
    F = as.matrix(F)
    U_F = do.call(rbind,U_F)
    rownames(U_F) = colnames(Z)

    # ----------------------- #
    # ---Save initial values- #
    # ----------------------- #
    current_state = list(
      k              = k,
      Lambda         = Lambda,
      tot_F_prec     = tot_F_prec,
      F_h2_index     = F_h2_index,
      F_h2           = F_h2,
      U_F            = U_F,
      F              = F,
      tot_Eta_prec   = tot_Eta_prec,
      resid_h2_index = resid_h2_index,
      resid_h2       = resid_h2,
      U_R            = U_R,
      B              = B,
      B_F            = B_F,
      XB             = XB,
      cis_effects    = cis_effects,
      nrun           = 0,
      total_time     = 0
    )
    return(current_state)
  })

  # Initialize parameters for Lambda_prior and B_prior (may be model-specific)
  BSFG_state$current_state = BSFG_state$priors$Lambda_prior$sampler(BSFG_state)
  BSFG_state$current_state = BSFG_state$priors$B_prior$sampler(BSFG_state)

  # Initialize Eta
  observation_model_state = run_parameters$observation_model(run_parameters$observation_model_parameters,BSFG_state)
  BSFG_state$current_state[names(observation_model_state$state)] = observation_model_state$state

  BSFG_state$Posterior$posteriorSample_params = unique(c(BSFG_state$Posterior$posteriorSample_params,observation_model_state$posteriorSample_params))
  BSFG_state$Posterior$posteriorMean_params = unique(c(BSFG_state$Posterior$posteriorMean_params,observation_model_state$posteriorMean_params))

  return(BSFG_state)
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
  tempfile = 'diagnostics_plots_temp.pdf'
  pdf(tempfile)
  if(BSFG_state$run_parameters$simulation){
    plot_diagnostics_simulation(BSFG_state)
  } else{
    plot_diagnostics(BSFG_state)
  }
  dev.off()
  system(sprintf('mv %s %s',tempfile,file))
}
