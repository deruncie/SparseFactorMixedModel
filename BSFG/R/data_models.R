
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
#' @return posteriorSample_params a list of parameter names to record posterior samples of
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
        Eta_mean = X %*% B + F %*% t(Lambda) + Z %*% U_R
        resid_Eta_prec = tot_Eta_prec / (1-resid_h2)
        resids = matrix(rnorm(p*n,0,sqrt(1/resid_Eta_prec)),nr = n,nc = p,byrow=T)
      }
      missing_indices = which(Y_missing)
      Eta[missing_indices] = Eta_mean[missing_indices] + resids[missing_indices]
    }
    return(list(Eta = as.matrix(Eta)))
  })
  return(list(state = data_model_state[new_variables],
              posteriorSample_params = c(),
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
      Eta_mean = X %*% B + F %*% t(Lambda) + Z %*% U_R
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
              posteriorSample_params = c(),
              posteriorMean_params = c('Eta')
             )
         )
}


#' Sample Eta given B-splines individual-level model
#'
#' Eta is a matrix of individual-level parmaters for a B-spline with
#'    equivalent knots over all individuals
#'
#' This function should pre-calculate design matrices for each individual (assuming they are all
#'     unique), as well as SVDs for efficient repeated solvings of the posterior distributions
#'     as resid_Y_prec, Eta_mean and resid_Eta_prec are updated each iteraction.
#'
#' @param Y a vector of observation. Pre-scaled and centered if desired.
#' @param data_model_parameters a list including:
#'     \code{observations}, a data.frame with columns: \code{ID}, and \code{covariate}, ordered by ID.
#'     and the variables df, knots, degree and intercept.
#'     Empty values will use the defaults for \code{bs()}.
bs_model = function(Y,data_model_parameters,BSFG_state = list()){
  current_state = BSFG_state$current_state
  data_matrices = BSFG_state$data_matrices

  data_model_state = with(c(data_model_parameters,data_matrices,current_state),{

    if(!exists('model_matrices')){
      if(!exists('df') || !is.numeric(df)) df = NULL
      if(!exists('knots') || !is.numeric(knots)) knots = NULL
      if(!exists('degree') || !is.numeric(degree)) degree = 3
      if(!exists('intercept') || !is.numeric(intercept)) intercept = FALSE
      covariate = observations$covariate
      coefficients = splines::bs(covariate,df = df,knots=knots,degree=degree,intercept = intercept)

      model_matrices = tapply(1:nrow(observations),observations$ID,function(x) {
        list(
          X = coefficients[x,,drop=FALSE],
          y = Y[x,,drop=FALSE]
        )
      })
    }

    n = length(model_matrices)
    p = ncol(coefficients)

    if(length(current_state) == 0){
      Eta_mean = matrix(0,n,p)
      resid_Eta_prec = matrix(1e-10,1,p)
      resid_Y_prec = matrix(1)  # only a single precision parameter for the data_model?
    } else{
      Eta_mean = X %*% B + F %*% t(Lambda) + Z %*% U_R
      resid_Eta_prec = tot_Eta_prec / (1-resid_h2)
      if(!exists('resid_Y_prec')) resid_Y_prec = matrix(1)
    }

    if(!exists('ncores')) ncores = 1
    Eta = do.call(cbind,mclapply(1:n,function(i) {
      X = model_matrices[[i]]$X
      y = model_matrices[[i]]$y
      eta_i = sample_coefs_parallel_sparse_c(y,X,matrix(0),resid_Y_prec,rep(1,nrow(X)),matrix(Eta_mean[i,]),t(resid_Eta_prec),1)
    },mc.cores = ncores))
    Eta = t(Eta)

    Y_fitted = do.call(c,mclapply(1:n,function(i) {
      y_fitted = model_matrices[[i]]$X %*% Eta[i,]
      if(is(y_fitted,'Matrix')) y_fitted = y_fitted@x
      y_fitted
    },mc.cores = ncores))

    Y_tilde = Y - Y_fitted

    resid_Y_prec = matrix(rgamma(1,shape = resid_Y_prec_shape + 0.5*n*p, rate = resid_Y_prec_rate + 0.5*sum(Y_tilde^2)))

    return(list(Eta = Eta, resid_Y_prec = resid_Y_prec, model_matrices = model_matrices, coefficients = coefficients))
  })
  return(list(state = data_model_state,
              posteriorSample_params = c(),
              posteriorMean_params = c('Eta','resid_Y_prec')
  )
  )
}

#' Sample Eta given B-splines individual-level model with binomial observations
#'
#' Eta is a matrix of individual-level parmaters for a B-spline with
#'    equivalent knots over all individuals
#'
#' This function should pre-calculate design matrices for each individual (assuming they are all
#'     unique)
#'
#' @param Y a vector of observation. Pre-scaled and centered if desired.
#' @param data_model_parameters a list including:
#'     \code{observations}, a data.frame with columns: \code{ID}, \code{N} and \code{covariate}, ordered by ID.
#'     and the variables df, knots, degree and intercept.
#'     Empty values will use the defaults for \code{bs()}.
bs_binomial_model = function(Y,data_model_parameters,BSFG_state = list()){
  current_state = BSFG_state$current_state
  data_matrices = BSFG_state$data_matrices

  data_model_state = with(c(data_model_parameters,data_matrices,current_state),{

  log_binom = function(beta,X,y,N,mu, sigma2){
    Xbeta = X %*% beta
    eXbeta = exp(Xbeta)
    p = eXbeta/(1+eXbeta)
    return( sum(y*log(p) + (N-y)*log(1-p))
            - sum((beta - mu)^2)/(2*sigma2)
    )
  }

  if(!exists('model_matrices')){
    if(!exists('df') || !is.numeric(df)) df = NULL
    if(!exists('knots') || !is.numeric(knots)) knots = NULL
    if(!exists('degree') || !is.numeric(degree)) degree = 3
    if(!exists('intercept') || !is.numeric(intercept)) intercept = FALSE
    covariate = observations$covariate
    coefficients = splines::bs(covariate,df = df,knots=knots,degree=degree,intercept = intercept)

    model_matrices = tapply(1:nrow(observations),observations$ID,function(x) {
      list(
        X = Matrix(coefficients[x,]),
        y = Y[x,],
        N = observations$N[x]
      )
    })
  }

  n = length(model_matrices)
  p = ncol(coefficients)

  if(length(current_state) == 0){
    Eta_mean = matrix(0,n,p)
    resid_Eta_prec = matrix(1e-10,1,p)
  } else{
    Eta_mean = X %*% B + F %*% t(Lambda) + Z %*% U_R
    resid_Eta_prec = tot_Eta_prec / (1-resid_h2)
  }

  if(!exists(ncores)) ncores = 1
  Eta = do.call(rbind,mclapply(1:n,function(i) {
    X = model_matrices[[i]]$X
    y = model_matrices[[i]]$y
    N = model_matrices[[i]]$N
    eta_i = MfUSampler::MfU.Sample(Eta[i,], f=log_binom, uni.sampler="slice", X=X, y=y,N=N)
  },mc.cores = ncores))

    return(list(Eta = Eta, model_matrices = model_matrices, coefficients = coefficients))
  })
  return(list(state = data_model_state,
              posteriorSample_params = c(),
              posteriorMean_params = c('Eta')
  )
  )
}

probe_gene_model = function(Y,data_model_parameters,BSFG_state = list()){
  current_state = BSFG_state$current_state
  data_matrices = BSFG_state$data_matrices
  new_variables = c('Eta','resid_Y_prec','mu_probe')
  data_model_state = with(c(data_model_parameters,data_matrices,current_state),{
    p = nrow(Z_Y)
    n = nrow(Y)

    if(!exists('Eta') || !all(dim(Eta) == c(n,p))) Eta = matrix(0,n,p)
    if(!exists('resid_Y_prec')) resid_Y_prec = rep(1,ncol(Y))

    Y_tilde = Y - Eta %*% Z_Y
    mu_probe = colMeans(Y_tilde) + rnorm(ncol(Y))/sqrt(n*resid_Y_prec)

    Y_tilde = sweep(Y,2,mu_probe,'-')
    Y_std = sweep(Y_tilde,2,resid_Y_prec,'*')

    if(!exists('tot_Eta_prec')){
      sum_Y_prec = Z_Y %*% resid_Y_prec
      prec = sum_Y_prec@x
      post_mean = sweep(Y_std %*% t(Z_Y),2,prec,'/')
      resid = sweep(matrix(rnorm(n*p),n,p),2,sqrt(prec),'/')
    } else{
      resid_Eta_prec = tot_Eta_prec / (1-resid_h2)
      Eta_mean = X %*% B + F %*% t(Lambda) + Z %*% E_a

      Eta_mean_std = sweep(Eta_mean,2,resid_Eta_prec,'*')
      sum_Y_prec = Z_Y %*% resid_Y_prec
      prec = sum_Y_prec@x + resid_Eta_prec

      post_mean = sweep(Y_std %*% t(Z_Y) + Eta_mean_std,2,prec,'/')
      resid = sweep(matrix(rnorm(n*p),n,p),2,sqrt(prec),'/')
    }
    Eta = as.matrix(post_mean) + resid

    Y_tilde = Y_tilde - Eta %*% Z_Y
    resid_Y_prec = rgamma(ncol(Y),shape = resid_Y_prec_shape + 0.5*nrow(Y),
                          rate = resid_Y_prec_rate + 0.5 * colSums(Y_tilde^2))
    return(list(Eta = Eta, mu_probe = mu_probe, resid_Y_prec = resid_Y_prec))
  })
  return(list(state = data_model_state[new_variables],
              posteriorSample_params = c(),
              posteriorMean_params = c('Eta','resid_Y_prec','mu_probe')
  )
  )
}


#' Sample Eta given gene-specific cis-eQTL genotypes
#'
#' Eta is a matrix of gene expression residuals after accounting for cis-eQTL
#'
#' This function should pre-calculate design matrices for each individual (assuming they are all
#'     unique), as well as SVDs for efficient repeated solvings of the posterior distributions
#'     as resid_Y_prec, Eta_mean and resid_Eta_prec are updated each iteraction.
#'
#' @param Y a vector of observation. Pre-scaled and centered if desired.
#' @param data_model_parameters a list of model matrices for the cis-genotypes of each gene.
#'      If named, the names should correspond to the column names of Y
# cis_eQTL_model = function(Y,data_model_parameters,BSFG_state = list()){
#   current_state = BSFG_state$current_state
#   data_matrices = BSFG_state$data_matrices
#
#   new_variables = c('Eta')
#   data_model_state = with(data_matrices,current_state),{
#
#     if(length(current_state) == 0){
#       # prep model matrices
#       if(!is.null(names(data_model_parameters))) {
#         if(is.null(colnames(Y)) & length(data_model_parameters) == ncol(Y)) colnames(Y) = colnames(data_model_parameters)
#         if(!all(colnames(Y) %in% colnames(data_model_parameters))) stop('column names of Y and cis-genotype names do not match')
#       } else{
#         if(is.null(colnames(Y))){
#           if(!length(data_model_parameters) == ncol(Y)) stop('wrong number of cis-genotypes provided')
#           colnames(Y) = names(data_model_parameters) = 1:ncol(Y)
#         }
#       }
#       model_matrices = lapply(colnames(Y),function(i) {
#         X = model.matrix(~0+data_model_parameters[[i]])
#         XtX = crossprod(X)
#         XtXiXt = solve(XtX) %*% t(X)
#         chol_XtX = chol(XtX)
#         return(list(X = X, XtXiXt = XtXiXt, chol_XtX = chol_XtX))
#       })
#       names(model_matrices) = colnames(Y)
#
#       cis_to_gene = names(model_matrices)[sapply(model_matrices,function(x) ncol(chol_XtX))]
#
#       p = ncol(Y)
#       n = nrow(Y)
#
#       Eta = Eta_mean = matrix(0,n,p)
#       resid_Eta_prec = matrix(0,1,p)
#       resid_Y_prec = matrix(rep(1,p),dimnames = list(NULL,colnames(Y)))  # only a single precision parameter for the data_model?
#     } else{
#       Eta_mean = X %*% B + F %*% t(Lambda) + Z %*% U_R
#       resid_Eta_prec = tot_Eta_prec / (1-resid_h2)
#       if(!exists(resid_Y_prec)) resid_Y_prec = matrix(rep(1,p),dimnames = list(NULL,colnames(Y)))
#     }
#
#     # sample cis_effects conditional on Eta
#     Y_tilde_Eta = Y - Eta
#     cis_effects = lapply(colnames(Y),function(i) {
#       cis_mean = model_matrices[[i]]$XtXiXt %*% Y_tilde_Eta[,i]
#       cis_effect = cis_mean + backsolve(model_matrices[[i]]$chol_XtX,rnorm(length(cis_mean))) * sqrt(resid_Y_prec[i])
#       cis_effect
#     })
#     cis_means = sapply(colnames(Y),function(i) model_matrices[[i]]$X %*% cis_effects[[i]])
#
#     # sample Eta conditional on cis_effects
#     Y_tilde_cis = Y - cis_means
#     Y_std = sweep(Y_tilde_cis,2,resid_Y_prec,'*')
#
#     Eta_mean_std = sweep(Eta_mean,2,resid_Eta_prec,'*')
#     prec = resid_Y_prec + resid_Eta_prec
#
#     post_mean = sweep(Y_std + Eta_mean_std,2,prec,'/')
#     resid = sweep(matrix(rnorm(n*p),n,p),2,sqrt(prec),'/')
#     Eta = as.matrix(post_mean) + resid
#   }
#
#   Y_tilde = Y_tilde - Eta %*% Z_Y
#   resid_Y_prec = rgamma(ncol(Y),shape = resid_Y_prec_shape + 0.5*nrow(Y),
#                         rate = resid_Y_prec_rate + 0.5 * colSums(Y_tilde^2))
#   return(list(Eta = Eta, mu_probe = mu_probe, resid_Y_prec = resid_Y_prec))
# })
#
#
#
#     if(!exists(ncores)) ncores = 1
#     Eta = do.call(rbind,mclapply(1:n,function(i) {
#       X = model_matrices[[i]]$X
#       y = model_matrices[[i]]$y
#       Cholesky_R = model_matrices[[i]]$Cholesky_R
#       chol_R = model_matrices[[i]]$chol_R
#       eta_i = sample_MME_single_diagA(y,X,X,Eta_mean[i,],resid_Eta_prec,Cholesky_R,chol_R,R_Perm = NULL,resid_Y_prec)
#     },mc.cores = ncores))
#
#     Y_fitted = do.call(c,mclapply(1:n,function(i) {
#       y_fitted = model_matrices[[i]]$X %*% t(Eta[i,])
#       if(is(y_fitted,'Matrix')) y_fitted = y_fitted@x
#       y_fitted
#     },mc.cores = ncores))
#
#     Y_tilde = Y - Y_fitted
#
#     resid_Y_prec = rgamma(1,shape = resid_Y_prec_shape + 0.5*n*p, rate = resid_Y_prec_rate + 0.5*sum(Y_tilde^2))
#
#     return(list(Eta = Eta, resid_Y_prec = resid_Y_prec, model.matrices = model.matrices, coefficients = coefficients))
#   })
#   return(list(state = data_model_state,
#               posteriorSample_params = c(),
#               posteriorMean_params = c('Eta','resid_Y_prec')
#   )
#   )
# }


