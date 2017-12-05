  #' Set BSFG run parameters
#'
#' Function to create run_parameters list for initializing BSFG model
#'
#' @param Posterior_folder path to folder to save posterior samples. Samples of each parameter
#'     are saved in chuncks to limit memory requirements.
#' @param simulation Is this a fit to simulated data? If so, a setup list will be expected providing
#'   the true values
#' @param scale_Y Should the Y values be centered and scaled? Recommend, except for simulated data.
#' @param lambda_propto_Vp Should the prior for lambda include tot_Eta_prec?
#' @param cauchy_sigma_tot Should the prior on sigma_j be half-cauchy?
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
#' @param num_NA_groups If 0, all NAs will be imputed during sampling. If Inf, all NAs will be marginalized over.
#'     If in (0,Inf), up to this many groups of columns will be separately sampled.
#'     The minimum number of NAs in each column not in one of these groups will be imputed.
#' @param svd_K If TRUE, and a 1-random effect model is specified, the the diagonalization of ZKZt is accomplished using this algorithm:
#'     https://math.stackexchange.com/questions/67231/singular-value-decomposition-of-product-of-matrices which doesn't require forming ZKTt.
#'     If FALSE, the SVD of ZKZt is calculated directly. TRUE is generally faster if the same genomes are repeated several times.
#' @seealso \code{\link{BSFG_init}}, \code{\link{sample_BSFG}}, \code{\link{print.BSFG_state}}
#'
BSFG_control = function(
                        Posterior_folder = "Posterior",
                        simulation = c(F,T),scale_Y = c(T,F),
                        lambda_propto_Vp = TRUE,cauchy_sigma_tot = FALSE,
                        b0 = 1, b1 = 0.0005, epsilon = 1e-1, prop = 1.00,
                        k_init = 20, h2_divisions = 100, h2_step_size = NULL,
                        drop0_tol = 1e-14, K_eigen_tol = 1e-10,
                        burn = 100,thin = 2,
                        num_NA_groups = Inf,
                        svd_K = TRUE,
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
#' @param tot_Y_var List of parameters of inverse gamma distribution for residual variances, specifically:
#'     \code{V} and \code{nu}, give shape = \code{nu/2} and scale = \code{nu*V/2}, so mean = \code{\frac{V*nu}{nu-2}}
#' @param tot_F_var List of parameters of inverse gamma distribution for factor variances. See \code{tot_Y_var}.
#'     This parameter provides the parameter extension of Ghosh and Dunson (2009), but is removed
#'     from all Factor parameters before they are saved in Posterior
#' @param h2_priors_resids_fun function for that returns prior probability for a given value of h2
#'     for each random effect. Should take two argument - a vector \code{h2} values for each random effect,
#'     and \code{n} - the number of discrete levels of the prior.
#'     Alternatively, can be a scalar or vector of (relative) prior values for each value of the
#'     discrete prior.
#' @param h2_priors_factors_fun see \code{h2_priors_resids_fun}. Same, but for the h2s of the factors.
#' @param Lambda_prior A list with elements:
#'     1) \code{sampler}: a function that draws samples of the precision matrix for Lambda. Ex: \code{sample_Lambda_prec_ARD}; 2)
#'     any other hyperparameters and control parameters for \code{sampler}
#' @param B_prior A list with elements:
#'     1) \code{sampler}: a function that draws samples of the precision matrix for B and B_F Ex: \code{sample_B_prec_ARD}; 2)
#'     any other hyperparameters and control parameters for \code{sampler}
#' @param QTL_prior A list with elements:
#'     1) \code{sampler}: a function that draws samples of the precision matrix for B_QTL and B_QTL_F Ex: \code{sample_QTL_prec_horseshoe}; 2)
#'     any other hyperparameters and control parameters for \code{sampler}
#'
#' @return a list with each of the prior components specified above.
#' @export
#'
#' @examples
BSFG_priors = function(
                        tot_Y_var = list(V = 0.5,   nu = 3),
                        tot_F_var = list(V = 18/20, nu = 20),
                        h2_priors_resids_fun = function(h2s, n) 1,
                        h2_priors_factors_fun = function(h2s, n) 1,
                        Lambda_prior = list(
                          sampler = sample_Lambda_prec_ARD,
                          Lambda_df = 3,
                          delta_1   = list(shape = 2,  rate = 1),
                          delta_2   = list(shape = 3, rate = 1),
                          delta_iteractions_factor = 100
                        ),
                        B_prior = list(
                          sampler = sample_B_prec_ARD,
                          global   = list(V = 1,nu = 3),
                          global_F = list(V = 1,nu = 3),
                          B_df      = 3,
                          B_F_df    = 3
                        ),
                        QTL_prior = list(
                          sampler = sample_QTL_prec_horseshoe,
                          separate_QTL_shrinkage=T,
                          cauchy_iteractions_factor = 10
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
#' Qt*Eta = Qt*X*B + Qt*F*t(Lambda) + Qt*ZL*U_R + Qt*E_R
#' Qt*F = Qt*X_F * B_F + Qt*ZL*U_F + Qt*E_F
#'
#' where LTL = K and ZL = Z*L
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
                     posteriorSample_params = c('Lambda','U_F','F','delta','tot_F_prec','F_h2','tot_Eta_prec','resid_h2', 'B', 'B_F','B_QTL','B_QTL_F','U_R','cis_effects'),
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
    if(nrow(Y) != nrow(data)) stop('Y and data have different numbers of rows')
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
	if(any(is.na(X))) stop('Missing values in X_resid')

  # for F
	# if factor_model_fixed not specified, use fixed effects from model for both
	if(is.null(factor_model_fixed)) {
	  factor_model_fixed = nobars(model)
	}

	# check that there are no random effects specifed in factor_model_fixed
	if(length(findbars(factor_model_fixed))) stop('Do not specify random effects in `factor_model_fixed`. Random effects for factors taken from `model`')

	X_F = model.matrix(factor_model_fixed,data)
	linear_combos = caret::findLinearCombos(X_F)
	if(!is.null(linear_combos$remove)) {
	  cat(sprintf('dropping column(s) %s to make X_factor full rank\n',paste(linear_combos$remove,sep=',')))
	  X_F = X_F[,-linear_combos$remove]
	}
	X_F = sweep(X_F,2,colMeans(X_F),'-') # note columns are centered, potentially resulting in zero-variance columns
	b_F = ncol(X_F)
	X_F_zero_variance = apply(X_F,2,var) == 0
	X_Fm = X_F
	if(any(is.na(X_F))) stop('Missing values in X_F')

	# identify fixed effects present in both X and X_F
	fixed_effects_common = matrix(0,nr=2,nc=0)  # 2 x b_c matrix with row1 indexes of X and row2 corresponding indexes of X_F
	non_null_X = which(apply(X,2,var)>0)
	non_null_XF = which(!X_F_zero_variance)
	if(length(non_null_X) > 0 && length(non_null_XF) > 0){
  	cor_X = abs(cor(X[,non_null_X,drop=FALSE],X_F[,non_null_XF,drop=FALSE]))
  	for(j in 1:nrow(cor_X)){
  	  if(any(cor_X[j,] == 1)){
        fixed_effects_common = cbind(fixed_effects_common,c(non_null_X[j],non_null_XF[which(cor_X[j,] == 1)[1]]))
      }
  	}
	}
	fixed_effects_only_resid = NULL
	fixed_effects_only_factors = NULL
	if(ncol(fixed_effects_common) < ncol(X)){
	  fixed_effects_only_resid = (1:ncol(X))[-fixed_effects_common[1,]]
	}
	if(ncol(fixed_effects_common) < ncol(X_F)){
	  fixed_effects_only_factors = (1:ncol(X_F))[-fixed_effects_common[2,]]
	}

	X = as(X,'dgCMatrix')
	X_F = as(X_F,'dgCMatrix')

	# -------- QTL effects ---------- #
	QTL_resid_Z = QTL_resid_X = NULL
	b_QTL = 0
	if(!is.null(QTL_resid)){
	  if(class(QTL_resid) == 'list'){
	    QTL_model = findbars(QTL_resid$model)
	    if(length(QTL_model) == 0) stop('no grouping factors found in QTL_resid model')
	    if(length(QTL_model) > 1)  stop('more than one grouping factor found in QTL_resid model')

	    group_model = as.formula(paste0('~',as.character(QTL_model[[1]][2])))
	    group_mm = model.matrix(group_model,data)

	    groupIDs = levels(as.factor(data[[as.character(QTL_model[[1]][[3]])]]))

	    QTL_terms = mkReTrms(QTL_model,data,drop.unused.levels = FALSE)
	    QTL_resid_Z = as(t(QTL_terms$Zt),'dgCMatrix')
	    QTL_resid_Z = QTL_resid_Z[,c(matrix(1:ncol(QTL_resid_Z),ncol = ncol(group_mm),byrow=T))]
	    if(!all(colnames(QTL_resid_Z) %in% rownames(QTL_resid$X))) stop(sprintf('Missing %s from QTL_resid$X',names(QTL_terms$cnms)[1]))
	    if(!all(colnames(group_mm) == QTL_terms$cnms[[1]])) stop("QTL_resid model didn't parse correctly. \nYou may have to create the X matrix yourself and use ~(1|group) as the model")
	    if(!all(colnames(QTL_resid_Z)[1:length(groupIDs)] == groupIDs)) stop("QTL_resid model didn't parse correctly. \nYou may have to create the X matrix yourself and use ~(1|group) as the model")

	    QTL_resid_X = do.call(bdiag,lapply(1:ncol(group_mm),function(i) QTL_resid$X[groupIDs,]))
	    QTL_resid_X = as.matrix(QTL_resid_X)

	    if(is.null(colnames(QTL_resid$X))) {
	      colnames(QTL_resid_X) = rep(paste0('m',1:ncol(QTL_resid$X)),ncol(group_mm))
	    } else {
	      colnames(QTL_resid_X) = rep(colnames(QTL_resid$X),ncol(group_mm))
	    }
	    colnames(QTL_resid_X) = paste(rep(colnames(group_mm),each = ncol(QTL_resid$X)),colnames(QTL_resid_X),sep='.')
	  } else{
	    if(nrow(QTL_resid) != nrow(data)) stop(sprintf('QTL_resid has wrong number of rows. Should be %d',nrow(data)))
	    QTL_resid_Z = as(diag(1,nrow(data)),'dgCMatrix')
	    if(is.data.frame(QTL_resid)) QTL_resid = as.matrix(QTL_resid)
	    QTL_resid_X = QTL_resid
	  }
	  b_QTL = ncol(QTL_resid_X)
	}

	QTL_factors_Z = QTL_factors_X = NULL
	b_QTL_F = 0
	if(!is.null(QTL_factors)){
	  if(class(QTL_factors) == 'list'){
	    QTL_model = findbars(QTL_factors$model)
	    if(length(QTL_model) == 0) stop('no grouping factors found in QTL_factors model')
	    if(length(QTL_model) > 1)  stop('more than one grouping factor found in QTL_factors model')

	    group_model = as.formula(paste0('~',as.character(QTL_model[[1]][2])))
	    group_mm = model.matrix(group_model,data)

	    groupIDs = levels(as.factor(data[[as.character(QTL_model[[1]][[3]])]]))

	    QTL_terms = mkReTrms(QTL_model,data,drop.unused.levels = FALSE)
	    QTL_factors_Z = as(t(QTL_terms$Zt),'dgCMatrix')
	    QTL_factors_Z = QTL_factors_Z[,c(matrix(1:ncol(QTL_factors_Z),ncol = ncol(group_mm),byrow=T))]
	    if(!all(colnames(QTL_factors_Z) %in% rownames(QTL_factors$X))) stop(sprintf('Missing %s from QTL_factors$X',names(QTL_terms$cnms)[1]))
	    if(!all(colnames(group_mm) == QTL_terms$cnms[[1]])) stop("QTL_factors model didn't parse correctly. \nYou may have to create the X matrix yourself and use ~(1|group) as the model")
	    if(!all(colnames(QTL_factors_Z)[1:length(groupIDs)] %in% groupIDs)) stop("QTL_factors model didn't parse correctly. \nYou may have to create the X matrix yourself and use ~(1|group) as the model")

	    QTL_factors_X = do.call(bdiag,lapply(1:ncol(group_mm),function(i) QTL_factors$X[groupIDs,]))
	    QTL_factors_X = as.matrix(QTL_factors_X)

	    if(is.null(colnames(QTL_factors$X))) {
	      colnames(QTL_factors_X) = rep(paste0('m',1:ncol(QTL_factors$X)),ncol(group_mm))
	    } else {
	      colnames(QTL_factors_X) = rep(colnames(QTL_factors$X),ncol(group_mm))
	    }
	    colnames(QTL_factors_X) = paste(rep(colnames(group_mm),each = ncol(QTL_factors$X)),colnames(QTL_factors_X),sep='.')
	  } else{
	    if(nrow(QTL_factors) != nrow(data)) stop(sprintf('QTL_factors has wrong number of rows. Should be %d',nrow(data)))
	    QTL_factors_Z = as(diag(1,nrow(data)),'dgCMatrix')
	    if(is.data.frame(QTL_factors)) QTL_factors = as.matrix(QTL_factors)
	    QTL_factors_X = QTL_factors
	  }
	  b_QTL_F = ncol(QTL_factors_X)
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

	# form Z
	Z = do.call(cbind,Z_matrices[RE_names])
	Z = as(Z,'dgCMatrix')

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
	ZL_matrices = list()

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
	    rownames(K_mats[[re]]) = levels(as.factor(data[[re]]))
	    K = K_mats[[re]]
	    L = Diagonal(nrow(K),1)
	  }
	  rownames(L) = paste(levels(as.factor(data[[re]])),re_name,sep='::')
	  K_mats[[re_name]] = fix_K(K)
	  RE_L_matrices[[re_name]] = L
	  ZL_matrices[[re_name]] = Z_matrices[[re_name]] %*% L
	}
	K_mats = K_mats[RE_names]
	# Construct ZL based on PSD K's
	ZL = do.call(cbind,ZL_matrices[RE_names])
	ZL = as(ZL,'dgCMatrix')
	# The following matrix is used to transform random effects back to the original space had we sampled from the original (PSD) K.
	if(length(RE_names) > 1) {
	  RE_L = do.call(bdiag,RE_L_matrices[RE_names])
	  rownames(RE_L) = do.call(c,lapply(RE_L_matrices,rownames))
	} else{
	  RE_L = RE_L_matrices[[1]]
	}
	r_RE = sapply(ZL_matrices,function(x) ncol(x))  # re-calculate

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
	if(run_parameters$num_NA_groups > 0) {
	  # columns with same patterns of missing data
	  Y_missing_mat = as.matrix(Y_missing)
    Y_col_obs = lapply(1:ncol(Y_missing_mat),function(x) {
      obs = which(!Y_missing_mat[,x],useNames=F)
      names(obs) = NULL
      obs
    })
    non_missing_rows = unname(which(rowSums(!Y_missing_mat)>0))
    unique_Y_col_obs = unique(c(list(non_missing_rows),Y_col_obs))
    unique_Y_col_obs_str = lapply(unique_Y_col_obs,paste,collapse='')
    Y_col_obs_index = sapply(Y_col_obs,function(x) which(unique_Y_col_obs_str == paste(x,collapse='')))

    if(length(unique_Y_col_obs) > run_parameters$num_NA_groups){
      col_counts = sort(tapply(Y_col_obs_index,Y_col_obs_index,length),decreasing = T)
      biggest_cols = as.numeric(names(col_counts))[1:run_parameters$num_NA_groups]
      biggest_cols = unique(c(1,biggest_cols))
      unique_Y_col_obs = unique_Y_col_obs[biggest_cols]
      Y_col_obs_index_new = rep(NA,length(Y_col_obs))
      for(i in 1:length(Y_col_obs_index)){
        if(Y_col_obs_index[i] %in% biggest_cols){
          Y_col_obs_index_new[i] = match(Y_col_obs_index[i],biggest_cols)
        } else{
          diffs = sapply(unique_Y_col_obs,function(x) sum(Y_col_obs[[i]] %in% x == F) + sum(x %in% Y_col_obs[[i]] == F))
          Y_col_obs_index_new[i] = order(diffs)[1]
          Y_missing_mat[,i] = T
          Y_missing_mat[unique_Y_col_obs[[Y_col_obs_index_new[i]]],i] = F
        }
      }
      Y_col_obs_index = Y_col_obs_index_new
    }

    Missing_data_map = lapply(seq_along(unique_Y_col_obs),function(i) {
      x = unique_Y_col_obs[[i]]
      return(list(
        Y_obs = x,
        Y_cols = which(Y_col_obs_index == i)
      ))
    })

    # rows with same patterns of missing data
    Y_row_obs = lapply(1:nrow(Y_missing_mat),function(x) {
      obs = which(!Y_missing_mat[x,],useNames=F)
      names(obs) = NULL
      obs
    })
    non_missing_cols = unname(which(colSums(!Y_missing_mat)>0))
    unique_Y_row_obs = unique(c(list(non_missing_cols),Y_row_obs))
    unique_Y_row_obs_str = lapply(unique_Y_row_obs,paste,collapse='')
    Y_row_obs_index = sapply(Y_row_obs,function(x) which(unique_Y_row_obs_str == paste(x,collapse='')))

    Missing_row_data_map = lapply(seq_along(unique_Y_row_obs),function(i) {
      x = unique_Y_row_obs[[i]]
      return(list(
        Y_cols = x,
        Y_obs = which(Y_row_obs_index == i)
      ))
    })
	} else{
	  Missing_data_map = list(list(
	    Y_obs = 1:n,
	    Y_cols = 1:p
	  ))
	  Missing_row_data_map = list(list(
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
  Qt_QTL_resid_Z_list = list()
  QtZL_list = list()
  Qt_cis_genotypes_list = list()
  randomEffect_C_Choleskys_list = list()
  Sigma_Choleskys_list = list()

  if(n_RE == 1 && run_parameters$svd_K == TRUE){
    svd_K1 = svd(K_mats[[1]])
  }

  for(set in seq_along(Missing_data_map)){
    x = Missing_data_map[[set]]$Y_obs
    cols = Missing_data_map[[set]]$Y_cols
    if(length(x) == 0) next

    if(n_RE == 1){
      if(run_parameters$svd_K == TRUE) {
        # a faster way of taking the SVD of ZLKZLt, particularly if ncol(ZL) < nrow(ZL). Probably no benefit if ncol(K) > nrow(ZL)
        qr_ZU = qr(ZL[x,,drop=FALSE] %*% svd_K1$u)
        R_ZU = drop0(qr.R(qr_ZU,complete=F),tol=run_parameters$drop0_tol)
        Q_ZU = drop0(qr.Q(qr_ZU,complete=T),tol=run_parameters$drop0_tol)
        RKRt = R_ZU %*% diag(svd_K1$d) %*% t(R_ZU)
        svd_RKRt = svd(RKRt)
        RKRt_U = svd_RKRt$u
        if(ncol(Q_ZU) > ncol(RKRt_U)) RKRt_U = bdiag(RKRt_U,diag(1,ncol(Q_ZU)-ncol(RKRt_U)))
        Qt = t(Q_ZU %*% RKRt_U)
      } else{
        ZKZt = ZL[x,,drop=FALSE] %*% K_mats[[1]] %*% t(ZL[x,,drop=FALSE])
        result = svd(ZKZt)
        Qt = t(result$u)
      }
      Qt = as(drop0(as(Qt,'dgCMatrix'),tol = run_parameters$drop0_tol),'dgCMatrix')
    } else{
      Qt = as(diag(1,length(x)),'dgCMatrix')
    }
    QtZL_matrices_set = lapply(ZL_matrices,function(ZL) Qt %*% ZL[x,,drop=FALSE])
    QtZL_set = do.call(cbind,QtZL_matrices_set[RE_names])
    QtZL_set = as(QtZL_set,'dgCMatrix')
    QtX_set = Qt %**% X[x,,drop=FALSE]
    if(!is.null(QTL_resid_Z)){
      Qt_QTL_resid_Z_set = as(drop0(Qt %*% QTL_resid_Z[x,,drop=FALSE],tol=run_parameters$drop0_tol),'dgCMatrix')
    } else{
      Qt_QTL_resid_Z_set = NULL
    }
    Qt_cis_genotypes_set = lapply(cis_genotypes[cols],function(X) Qt %**% X[x,,drop=FALSE])

    Qt_list[[set]]   = Qt
    QtX_list[[set]]  = QtX_set
    Qt_QTL_resid_Z_list[[set]] = Qt_QTL_resid_Z_set
    QtZL_list[[set]]  = QtZL_set
    Qt_cis_genotypes_list[[set]] = Qt_cis_genotypes_set

    ZKZts_set = list()
    for(i in 1:n_RE){
      ZKZts_set[[i]] = as(forceSymmetric(drop0(QtZL_matrices_set[[i]] %*% K_mats[[i]] %*% t(QtZL_matrices_set[[i]]),tol = run_parameters$drop0_tol)),'dgCMatrix')
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

    ZtZ_set = as(forceSymmetric(drop0(crossprod(ZL[x,]),tol = run_parameters$drop0_tol)),'dgCMatrix')

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

  # Qt matrices for factors are only used with row set 1
  x = Missing_data_map[[1]]$Y_obs
  Qt1_XF = Qt_list[[1]] %**% X_F[x,,drop=FALSE]
  if(!is.null(QTL_factors_Z)){
    Qt1_QTL_Factors_Z = as(drop0(Qt_list[[1]] %*% QTL_factors_Z[x,,drop=FALSE],tol = run_parameters$drop0_tol),'dgCMatrix')
  } else{
    Qt1_QTL_Factors_Z = NULL
  }


  run_variables = list(
    p      = p,
    n      = n,
    r_RE   = r_RE,
    RE_names = RE_names,
    b       = b,
    b_F     = b_F,
    b_QTL   = b_QTL,
    b_QTL_F = b_QTL_F,
    fixed_effects_common = fixed_effects_common,
    fixed_effects_only_resid = fixed_effects_only_resid,
    fixed_effects_only_factors = fixed_effects_only_factors,
    resid_intercept   = resid_intercept,
    X_F_zero_variance = X_F_zero_variance,   # used to identify fixed effect coefficients that should be forced to zero
    n_cis_effects     = n_cis_effects,
    cis_effects_index = cis_effects_index,
    Qt_list    = Qt_list,
    QtX_list   = QtX_list,
    Qt_QTL_resid_Z_list = Qt_QTL_resid_Z_list,
    QtZL_list   = QtZL_list,
    Qt_cis_genotypes_list = Qt_cis_genotypes_list,
    Qt1_XF  = Qt1_XF,
    Qt1_QTL_Factors_Z = Qt1_QTL_Factors_Z,
    Missing_data_map      = Missing_data_map,
    Missing_row_data_map  = Missing_row_data_map,
    Sigma_Choleskys_list          = Sigma_Choleskys_list,
    randomEffect_C_Choleskys_list = randomEffect_C_Choleskys_list
  )

	data_matrices = list(
	  X           = X,
	  X_F         = X_F,
	  Z_matrices  = Z_matrices,
	  Z           = Z,
	  ZL_matrices = ZL_matrices,
	  ZL          = ZL,
	  RE_L        = RE_L,  # matrix necessary to back-transform U_F and U_R (RE_L*U_F and RE_L*U_R) to get original random effects
	  RE_indices  = RE_indices,
	  h2s_matrix  = h2s_matrix,
	  cis_genotypes = cis_genotypes,
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
	# total precision
	if(length(priors$tot_Y_var$V == 1)) {
	  priors$tot_Y_var$V = rep(priors$tot_Y_var$V,p)
	  priors$tot_Y_varnu = rep(priors$tot_Y_var$nu,p)
	}
	priors$tot_Eta_prec_rate   = with(priors$tot_Y_var,V * nu)
	priors$tot_Eta_prec_shape  = with(priors$tot_Y_var,nu - 1)
	priors$tot_F_prec_rate     = with(priors$tot_F_var,V * nu)
	priors$tot_F_prec_shape    = with(priors$tot_F_var,nu - 1)

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
    rownames(U_R) = colnames(ZL)

    # Fixed effects
    B = 0*matrix(rnorm(b*p), ncol = p)
    colnames(B) = traitnames

    # Factor fixed effects
    B_F = 0*matrix(rnorm(b_F * k),b_F,k)

    # QTL effects
    B_QTL = 0*rstdnorm_mat(b_QTL,p)
    B_QTL_F = 0*rstdnorm_mat(b_QTL_F,k)

    # cis effects
    cis_effects = matrix(rnorm(cis_effects_index[length(cis_effects_index)]-1,0,1),nrow=1)

    XB = X %**% B
    if(b_QTL > 0) XB = XB + QTL_resid_Z %**% (QTL_resid_X %**% B_QTL)

    F = X_F %*% B_F + matrix(rnorm(n * k, 0, sqrt((1-colSums(F_h2)) / tot_F_prec)),ncol = k, byrow = T)
    for(effect in RE_names) {
      F = F + ZL_matrices[[effect]] %*% U_F[[effect]]
    }
    F = as.matrix(F)
    U_F = do.call(rbind,U_F)
    rownames(U_F) = colnames(ZL)

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
      B_QTL          = B_QTL,
      B_QTL_F        = B_QTL_F,
      XB             = XB,
      cis_effects    = cis_effects,
      nrun           = 0,
      total_time     = 0
    )
    return(current_state)
  })

  # Initialize parameters for Lambda_prior, B_prior, and QTL_prior (may be model-specific)
  BSFG_state$current_state = BSFG_state$priors$Lambda_prior$sampler(BSFG_state)
  BSFG_state$current_state = BSFG_state$priors$B_prior$sampler(BSFG_state)
  BSFG_state$current_state = BSFG_state$priors$QTL_prior$sampler(BSFG_state)

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
      c(sprintf('Model dimensions: fixed = %d, random = %d \n',ncol(data_matrices$X),ncol(data_matrices$ZL))),
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
