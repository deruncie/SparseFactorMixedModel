
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
#' @param BSFG_state a BSFG_state object. Generally, only current_state and data_matrices is used. If
#'    empty, will return a default set of parameters of the appropriate size for model initialization.
#' @return list of data_model variables including:
#' @return state a list of parameters associated with the data_model. Here, only the matrix Eta
#' @return sample_params a list of parameter names to record posterior samples of
#' @return posteriorMean_params a list of parameters to record the posterior mean, but not save individual
#'     posterior samples
missing_data_model = function(Y,data_model_parameters,BSFG_state = list()){
  current_state = BSFG_state$current_state
  data_matrices = BSFG_state$data_matrices

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
  return(list(state = data_model_state[new_variables],
              sample_params = c(),
              posteriorMean_params = c('Eta')
  ))
}

#' Sample Eta given a voom observation model
#'
#' \code{voom} or \code{voomWithQualityWeights} (limma) should be run on RNAseq data,
#'    resulting in logCPM (E) and inverseWeights (weights)
#'    for each observation. The observations (E) should be passed as Y, and the weights as
#'    data_model_parameters$prec_Y.
#'
#' When running \code{voom}, a fully-specified fixed effect model for the data should be specified.
#'
#' @inheritParams missing_data_model
#' @inheritSection @value
#' @references Law, C. W., Chen, Y., Shi, W., & Smyth, G. K. (2014).
#'     voom: precision weights unlock linear model analysis tools for RNA-seq read counts.
#'     Genome Biology, 15(2), R29. http://doi.org/10.1186/gb-2004-5-10-r80
voom_model = function(Y,data_model_parameters,BSFG_state = list()){
  current_state = BSFG_state$current_state
  data_matrices = BSFG_state$data_matrices

  new_variables = c('Eta')
  data_model_state = with(c(data_model_parameters,data_matrices,current_state),{
    Eta = Y
    if(length(current_state) > 0){
      n = nrow(Y)
      p = ncol(Y)
      Eta_mean = X %*% B + F %*% t(Lambda) + Z %*% E_a
      resid_Eta_prec = tot_Eta_prec / (1-resid_h2)
      prec = sweep(prec_Y,2,resid_Eta_prec,'+')
      # Y_std = Y * prec_Y. This is calculated once in the initialization.
      Eta_hat = (Y_std + sweep(Eta_mean,2,resid_Eta_prec,'*')) / prec
      resid = matrix(rnorm(n*p),n,p) / sqrt(prec)
      Eta = Eta_hat + resid
    }
    return(list(Eta = as.matrix(Eta)))
  })
  return(list(state = data_model_state[new_variables],
              sample_params = c(),
              posteriorMean_params = c('Eta')
             )
         )
}


#' Sample Eta given B-splines individual-level model
#'
#' **Not yet functional**.
#' Eta is a matrix of individual-level parmaters for a B-spline with
#'    equivalent knots over all individuals
#'
#' This function should pre-calculate design matrices for each individual (assuming they are all
#'     unique), as well as SVDs for efficient repeated solvings of the posterior distributions
#'     as resid_Y_prec, Eta_mean and resid_Eta_prec are updated each iteraction.
#'
#' @param Y a vector of observation. Pre-scaled and centered if desired.
#' @param data_model_parameters a list including a data.frame with columns: \code{ID} and \code{covariate}
#'     and the variables df, knots, degree and intercept. Empty values will use the defaults
#'     for \code{bs()}. The data.frame should be ordered by ID.
bs_model = function(Y,data_model_parameters,BSFG_state = list()){
  current_state = BSFG_state$current_state
  data_matrices = BSFG_state$data_matrices

  new_variables = c('Eta')
  data_model_state = with(c(data_model_parameters,data_matrices,current_state),{

    if(length(current_state) == 0){
      if(!exists(df)) df = NULL
      if(!exists(knots)) knots = NULL
      if(!exists(degree)) degree = NULL
      if(!exists(intercept)) intercept = NULL
      coefficients = bs(covariate,df = df,knots=knots,degree=degree,intercept = intercept)

      if(!exists(ncores)) ncores = 1

      model.matrices = tapply(1:nrow(coefficients),ID,function(x) {
        X = Matrix(coefficients[x,])
        XtX = crossprod(X)
        Diagonal(XtX) = rep(1,ncol(X))
        Cholesky_XtX = Cholesky(XtX)
        return(list(X=X,XtX = XtX))
        # svd_XtX = svd(crossprod(X))
        # return(list(U = Matrix, s = svd_XtX$d))
      })

      n = length(model.matrices)
      p = ncol(coefficients)

      Eta_mean = matrix(0,n,p)
      resid_Eta_prec = matrix(0,1,p)
      resid_Y_prec = 1  # only a single precision parameter for the data_model?
    } else{
      Eta_mean = X %*% B + F %*% t(Lambda) + Z %*% E_a
      resid_Eta_prec = tot_Eta_prec / (1-resid_h2)
      if(!exists(resid_Y_prec)) resid_Y_prec = matrix(1)
    }

    Eta = do.call(rbind,mclapply(1:n,function(i) {
      with(model.matrices[[i]],{
        ## should be able to use sample_MME_multiple_diagR, need to work out maths
        # Sigm_inv = XtX*precY;
        # Diagonal(Sigma_inv) = Diagonal(Sigma_inv) + resid_Y_prec
        # chol_Sigma_inv = update(Cholesky_XtX,Sigma_inv)
        # Eta_hat = backXtY*precY + sweep(Eta_mean,2,resid_Eta_prec,'*')
      })
    },mc.cores = ncores))

    Y_fitted = do.call(c,mclapply(1:n,function(i) {
      model.matrices[[i]] %*% t(Eta[i,])
    },mc.cores = ncores))

    Y_tilde = Y - Y_fitted

    resid_Y_prec = rgamma(1,shape = resid_Y_prec_shape + 0.5*n*p, rate = resid_Y_prec_rate + 0.5*sum(Y_tilde^2))

    return(list(Eta = Eta, resid_Y_prec = resid_Y_prec, model.matrices = model.matrices))
  })
  return(list(state = data_model_state,
              sample_params = c(),
              posteriorMean_params = c('Eta','resid_Y_prec')
  )
  )
}
