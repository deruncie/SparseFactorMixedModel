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
#'    in \code{observation_model_parameters}.
#'
#'    When run with an empty \code{current_state} and \code{data_matrices == NULL}, it should return
#'    a matrix Eta of the correct dimensions, but the values are unimportant.
#'
#' @param Y data matrix n_Y x p_Y
#' @param observation_model_parameters List of parameters necessary for the data model.
#'      Here, a Matrix of coordinates of NAs in Y
#' @param BSFG_state a BSFG_state object. Generally, only current_state and data_matrices is used. If
#'    empty, will return a default set of parameters of the appropriate size for model initialization.
#' @return list of data_model variables including:
#' @return state a list of parameters associated with the data_model. Here, only the matrix Eta
#' @return posteriorSample_params a list of parameter names to record posterior samples of
#' @return posteriorMean_params a list of parameters to record the posterior mean, but not save individual
#'     posterior samples
missing_data_model = function(observation_model_parameters,BSFG_state = list()){
  current_state = BSFG_state$current_state
  data_matrices = BSFG_state$data_matrices

  if(!'observation_setup' %in% names(observation_model_parameters)) {
    observation_setup = with(c(observation_model_parameters,data_matrices,current_state),{
      if(scale_Y){
        Mean_Y = colMeans(Y,na.rm=T)
        VY = apply(Y,2,var,na.rm=T)
        Y = sweep(Y,2,Mean_Y,'-')
        Y = sweep(Y,2,sqrt(VY),'/')
      } else {
        p_Y = dim(Y)[2]
        Mean_Y = rep(0,p_Y)
        VY = rep(1,p_Y)
      }
      Y_missing = as(is.na(Y),'lgTMatrix')# un-compressed logical sparse matrix
      return(list(
        n = nrow(Y),
        p = ncol(Y),
        traitnames = colnames(Y),
        Mean_Y = Mean_Y,
        VY = VY,
        Y_missing = Y_missing,
        n_missing = sum(Y_missing),
        missing_indices = which(Y_missing)
      ))
    })
    return(observation_setup)
  }

  observation_model_state = with(c(observation_model_parameters,observation_model_parameters$observation_setup,data_matrices,current_state),{
    Eta = Y
    Eta_mean = matrix(0,0,0)
    if(n_missing > 0){
      n = nrow(Y)
      p = ncol(Y)
      if(length(current_state) == 0) {
        Eta_mean = matrix(0,n,p)
        resids = rnorm(n_missing)
      } else{
        Eta_mean = XB + F %*% t(Lambda) + toDense(Z %*% U_R)
        resid_Eta_prec = tot_Eta_prec / (1-colSums(resid_h2))
        resids = rnorm(n_missing,0,sqrt(1/resid_Eta_prec[Y_missing@j+1]))  # sample resids from normal distribution with appropriate variance
      }
      Eta[missing_indices] = Eta_mean[missing_indices] + resids
    }
    return(list(Eta = Eta,Eta_mean = Eta_mean,Y_missing = Y_missing))
  })
  return(list(state = observation_model_state,
              posteriorSample_params = c('Eta'),
              posteriorMean_params = c('Eta_mean')
  ))
}



#' Sample Eta given a voom observation model
#'
#' \code{voom} or \code{voomWithQualityWeights} (limma) should be run on RNAseq data,
#'    resulting in logCPM (E) and inverseWeights (weights)
#'    for each observation. The observations (E) should be passed as Y, and the weights as
#'    observation_model_parameters$prec_Y.
#'
#' When running \code{voom}, a fully-specified fixed effect model for the data should be specified.
#'
#' @inheritParams missing_data_model
#' @inheritSection @value
#' @references Law, C. W., Chen, Y., Shi, W., & Smyth, G. K. (2014).
#'     voom: precision weights unlock linear model analysis tools for RNA-seq read counts.
#'     Genome Biology, 15(2), R29. http://doi.org/10.1186/gb-2004-5-10-r80
voom_model = function(observation_model_parameters,BSFG_state = list()){
  current_state = BSFG_state$current_state
  data_matrices = BSFG_state$data_matrices

  new_variables = c('Eta')
  observation_model_state = with(c(observation_model_parameters,data_matrices,current_state),{
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
  return(list(state = observation_model_state[new_variables],
              posteriorSample_params = c(),
              posteriorMean_params = c('Eta')
             )
         )
}


#' Sample Eta given regression-splines individual-level model
#'
#' Eta is a matrix of individual-level parmaters for a regression-spline with
#'    equivalent knots over all individuals
#'
#' This function should pre-calculate design matrices for each individual (assuming they are all
#'     unique). During sampling, should sample regression coefficients \code{Eta} given the factor model state.
#'
#' @param observation_model_parameters a list including:
#'     1) \code{observations}, a data.frame with observation-level data including columns \code{ID} and \code{Y}
#'     3) \code{individual_model} the model that should be applied to the data for each ID.
#' @param BSFG_state The current \code{BSFG_state}.
#'     For initialization, can be a list with \code{data_matrices = list(data=data)}
#'     Where \code{data} contains \code{ID} for ordering.
regression_model = function(observation_model_parameters,BSFG_state = list()){
  current_state = BSFG_state$current_state
  data_matrices = BSFG_state$data_matrices

  if(!'observation_setup' %in% names(observation_model_parameters)) {
    observation_setup = with(c(observation_model_parameters,data_matrices,current_state),{
      if(!'ID' %in% colnames(data)) stop('ID column required in data')
      if(!length(unique(data$ID)) == nrow(data)) stop('duplicate IDs in data')

      # ensure ID is a character for matching
      observations$ID = as.character(observations$ID)
      data$ID = as.character(data$ID)
      observations = subset(observations,ID %in% data$ID)

      # extract model Terms and Y matrix
      mf = model.frame(individual_model,observations)
      mm = model.matrix(individual_model,observations)
      traits = all.vars(update(individual_model,'~.0'))
      Y = as.matrix(observations[,traits,drop=FALSE])
      Terms = delete.response(terms(mf))

      # smooth bs() terms
      mm_terms = attr(mm,'assign')
      bs_terms = grep('bs(',attr(Terms,'term.labels'),fixed=T)
      mm_rot_list = list()
      for(i in 0:max(attr(mm,'assign'))){
        U = diag(1,sum(mm_terms==i))
        colnames(U) = colnames(mm)[mm_terms == i]
        if(i %in% bs_terms && (!'do_not_penalize_bs' %in% ls() || !do_not_penalize_bs)) {
          mm_cols = mm_terms == i

          # # construct a differences matrix to smooth the B-spline
          # lower-diagonal
          U = matrix(0,sum(mm_cols),sum(mm_cols))
          diag(U) = 1
          U[lower.tri(U)] = 1

          colnames(U) = paste0(colnames(mm)[mm_terms == i][1:ncol(U)],'_differences')
        }
        mm_rot_list[[i+1]] = U
      }
      mm_rotation = as.matrix(do.call(bdiag,mm_rot_list))
      colnames(mm_rotation) = do.call(c,lapply(mm_rot_list,colnames))
      mm = mm %*% mm_rotation
      n_terms = ncol(mm)

      id_index = tapply(1:nrow(observations),observations$ID,function(x) x)
      model_matrices = lapply(data$ID,function(id) {
        x = id_index[[id]]
        if(length(x) > 0){
          X = mm[x,,drop=FALSE]
        } else{
          x = matrix(0,0,0)
          X = matrix(0,0,n_terms)
        }
        # keep only columns of X which are non-zero, but record which columns these were of original X.
        nonZero_cols_X = which(colSums(X!=0) > 0)
        X = X[,nonZero_cols_X,drop=FALSE]
        list(
          X = X,
          y = Y[x,,drop=FALSE],
          position = x,
          nonZero_cols_X = nonZero_cols_X
        )
      })
      names(model_matrices) = data$ID

      n = length(model_matrices)
      n_traits = ncol(model_matrices[[1]]$y) # number of traits
      p_trait = n_terms  # number of coefficients per trait
      p = p_trait * n_traits # total number of coefficients
      traits = colnames(model_matrices[[1]]$y)
      traitnames = paste(rep(traits,each = p_trait),rep(colnames(mm),length(traits)),sep='::')
      Eta_row_names = data$ID
      Y_missing = t(sapply(model_matrices,function(x) rep(!(seq_len(n_terms) %in% x$nonZero_cols_X),n_traits)))

      observation_setup = list(
        Y = Y,
        n = length(model_matrices),
        p = p,
        Terms = Terms,
        n_traits = n_traits,
        traitnames = traitnames,
        Eta_row_names = Eta_row_names,
        model_matrices = model_matrices,
        mm_rotation = mm_rotation,
        Y_missing = Y_missing
      )
      return(observation_setup)
    })
    return(observation_setup)
  }

  observation_model_state = with(c(observation_model_parameters,observation_model_parameters$observation_setup,data_matrices,current_state),{
    if(!'var_Eta' %in% ls()) var_Eta = rep(1,p)
    if(length(current_state) == 0){
      Eta_mean = matrix(0,n,p)
      resid_Eta_prec = matrix(1,1,p)
      resid_Y_prec = matrix(rgamma(n_traits,shape = resid_Y_prec_shape,rate = resid_Y_prec_rate),nr=1) # only a single precision parameter for the data_model?
    } else{
      Eta_mean = X %**% B + F %*% t(Lambda) + Z %**% U_R
      resid_Eta_prec = tot_Eta_prec / (1-colSums(resid_h2))
      if(!'resid_Y_prec' %in% ls()) resid_Y_prec = matrix(rgamma(n_traits,shape = resid_Y_prec_shape,rate = resid_Y_prec_rate),nr=1) # only a single precision parameter for the data_model?
    }

    # re-scale Eta_mean and resid_Eta_prec
    Eta_mean = sweep(Eta_mean,2,sqrt(var_Eta),'*')
    resid_Eta_prec[] = resid_Eta_prec / var_Eta
    coefs = sample_coefs_set_c(model_matrices,resid_Y_prec, t(Eta_mean),matrix(resid_Eta_prec,length(resid_Eta_prec),n),1)

    Y_fitted = get_fitted_set_c(model_matrices,coefs,1)
    Eta = t(coefs)
    colnames(Eta) = traitnames
    rownames(Eta) = Eta_row_names

    # un-scale Eta
    Eta = sweep(Eta,2,sqrt(var_Eta),'/')

    Y_tilde = Y - Y_fitted

    resid_Y_prec = matrix(rgamma(n_traits,shape = resid_Y_prec_shape + 0.5*nrow(Y_tilde), rate = resid_Y_prec_rate + 0.5*colSums(Y_tilde^2)),nr=1)

    return(list(Eta = Eta, resid_Y_prec = resid_Y_prec, Y_fitted=Y_fitted, Y=Y, var_Eta = var_Eta))
  })
  return(list(state = observation_model_state,
              posteriorSample_params = c('Y_fitted','Eta','resid_Y_prec'),
              posteriorMean_params = c()
  )
  )
}


#' Sample cis_eQTL coefficients
#'
#' @param observation_model_parameters list with:
#'     \code{Y} gene expression data matrix n x p
#'     \code{cis_genotypes} a list of design matrices of length \code{p} (ie number of columns of \code{Eta})
#'     This is used to specify trait-specific fixed effects, such a cis-genotypes
#'
#' @param BSFG_state
#'
#' @return list including:
#'     \code{state} parameters to add to \code{current_state}
#'     \code{posteriorSample_params} character vector of parameter names to include in the list of parameters to record posterior samples
#'     \code{posteriorMean_params} character vector of parameter names to include in the list of parameters to record posterior mean
cis_eQTL_model = function(observation_model_parameters,BSFG_state = list()){
  current_state = BSFG_state$current_state
  data_matrices = BSFG_state$data_matrices

  observation_model_state = with(c(observation_model_parameters,data_matrices,current_state),{
    n = nrow(Y)
    p = ncol(Y)
    Ut_cis = as(diag(1,n),'dgCMatrix')
    s_cis = rep(1,n)
    if(!exists('Eta')) {
      Eta_mean = matrix(0,n,p)
    } else{
      Eta_mean = XB + F %*% t(Lambda) + Z %*% U_R
    }
    if(!exists('tot_Eta_prec')) {
      resid_Eta_prec = matrix(1,1,p)
    } else{
      resid_Eta_prec = tot_Eta_prec / (1-colSums(resid_h2))
    }

    Y_tilde = Y - Eta_mean
    if(!exists('ncores')) ncores = 1
    cis_effects_list = mclapply(1:p,function(j) {
      cis_X_j = cis_genotypes[[j]]
      b_j = ncol(cis_X_j)
      if(b_j > 0) {
        prior_mean = matrix(0,b_j,1)
        prior_prec = matrix(1e-10,b_j,1)
        prior_prec[apply(cis_X_j,2,var)==0] = 1e10
        randn_theta = matrix(rnorm(b_j),b_j,1)
        randn_e = matrix(rnorm(n),n,1)
        coefs_j = sample_coefs_parallel_sparse_c_Eigen(Ut_cis,Y_tilde[,j],cis_X_j,
                                           0, resid_Eta_prec[,j],
                                           s_cis,prior_mean,prior_prec,
                                           randn_theta,randn_e,
                                           1)
        return(coefs_j)
      }
      return(NULL)
    },mc.cores = ncores)

    cis_fitted = do.call(cbind,mclapply(1:p,function(j) {
      cis_X_j = cis_genotypes[[j]]
      if(ncol(cis_X_j) > 0) {
        return(cis_X_j %*% cis_effects_list[[j]])
      }
      return(rep(0,n))
    },mc.cores = ncores))

    Eta = Y - cis_fitted
    cis_effects = matrix(do.call(c,cis_effects_list),nrow=1)
    return(list(Eta = Eta, cis_effects = cis_effects))
  })
  return(list(state = observation_model_state,
              posteriorSample_params = c('Eta','cis_effects'),
              posteriorMean_params = c()
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
#' @param observation_model_parameters a list including:
#'     \code{observations}, a data.frame with columns: \code{ID}, \code{N} and \code{covariate}, ordered by ID.
#'     and the variables df, knots, degree and intercept.
#'     Empty values will use the defaults for \code{bs()}.
bs_binomial_model = function(observation_model_parameters,BSFG_state = list()){
  current_state = BSFG_state$current_state
  data_matrices = BSFG_state$data_matrices

  observation_model_state = with(c(observation_model_parameters,data_matrices,current_state),{

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
  return(list(state = observation_model_state,
              posteriorSample_params = c(),
              posteriorMean_params = c('Eta')
  )
  )
}

probe_gene_model = function(observation_model_parameters,BSFG_state = list()){
  current_state = BSFG_state$current_state
  data_matrices = BSFG_state$data_matrices
  new_variables = c('Eta','resid_Y_prec','mu_probe')
  observation_model_state = with(c(observation_model_parameters,data_matrices,current_state),{
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
  return(list(state = observation_model_state[new_variables],
              posteriorSample_params = c(),
              posteriorMean_params = c('Eta','resid_Y_prec','mu_probe')
  )
  )
}



#' Sample Eta given a list of different observation_models
#'
#'
#' @param observation_model_parameters a list of observation_model lists
#' @param BSFG_state the current BSFG_state object. Can be just a list with \code{data}
#'
#' @return
#' @export
#'
#' @examples
combined_model = function(observation_model_parameters,BSFG_state = list()){

  ## NOTE: Somehow need to pass a subset of the factor model to each sub-model, or maybe calculate Eta_mean once and pass a subset to each model.

  sub_models = observation_model_parameters$sub_models
  n_models = length(sub_models)
  if(is.null(names(sub_models))) names(sub_models) = 1:n_models

  # run each of the separate observation models
  results = lapply(names(sub_models),function(i) {
    Y = sub_models[[i]]
    if(!'observation_model' %in% names(Y)) stop(sprintf('observation_model not specified in model %s',i))
    observation_model = Y$observation_model

    # adjust model-specific parameter names
    if(length(names(BSFG_state$current_state))>0){
      names(BSFG_state$current_state) = sub(sprintf('.%s',i),'',names(BSFG_state$current_state))
    }

    # run observation_model
    new_state = observation_model(Y[names(Y) != 'observation_model'],BSFG_state)

    # fix names of parameters
    names(new_state$state) = paste(names(new_state$state),i,sep='.')
    if(length(new_state$posteriorSample_params)>0) new_state$posteriorSample_params = paste(new_state$posteriorSample_params,i,sep='.')
    if(length(new_state$posteriorMean_params)>0) new_state$posteriorMean_params = paste(new_state$posteriorMean_params,i,sep='.')

    # return new_state
    new_state
  })
  names(results) = names(sub_models)

  # merge Etas
  Eta = do.call(cbind,lapply(names(results),function(x) results[[x]]$state[[paste('Eta',x,sep='.')]]))

  # merge Y_missing
  Y_missing = do.call(cbind,lapply(names(results),function(x) results[[x]]$state[[paste('Y_missing',x,sep='.')]]))

  # make combined new_state
  new_state = do.call('c',lapply(results,function(x) x$state))
  names(new_state) = do.call('c',lapply(results,function(x) names(x$state)))
  new_state$Eta = Eta
  new_state$Y_missing = Y_missing

  posteriorSample_params = do.call('c',lapply(results,function(x) x$posteriorSample_params))
  posteriorMean_params = do.call('c',lapply(results,function(x) x$posteriorMean_params))

  return(list(state = new_state,
              posteriorSample_params = posteriorSample_params,
              posteriorMean_params = posteriorMean_params
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
#' @param observation_model_parameters a list of model matrices for the cis-genotypes of each gene.
#'      If named, the names should correspond to the column names of Y
# cis_eQTL_model = function(Y,observation_model_parameters,BSFG_state = list()){
#   current_state = BSFG_state$current_state
#   data_matrices = BSFG_state$data_matrices
#
#   new_variables = c('Eta')
#   observation_model_state = with(data_matrices,current_state),{
#
#     if(length(current_state) == 0){
#       # prep model matrices
#       if(!is.null(names(observation_model_parameters))) {
#         if(is.null(colnames(Y)) & length(observation_model_parameters) == ncol(Y)) colnames(Y) = colnames(observation_model_parameters)
#         if(!all(colnames(Y) %in% colnames(observation_model_parameters))) stop('column names of Y and cis-genotype names do not match')
#       } else{
#         if(is.null(colnames(Y))){
#           if(!length(observation_model_parameters) == ncol(Y)) stop('wrong number of cis-genotypes provided')
#           colnames(Y) = names(observation_model_parameters) = 1:ncol(Y)
#         }
#       }
#       model_matrices = lapply(colnames(Y),function(i) {
#         X = model.matrix(~0+observation_model_parameters[[i]])
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
#   return(list(state = observation_model_state,
#               posteriorSample_params = c(),
#               posteriorMean_params = c('Eta','resid_Y_prec')
#   )
#   )
# }


