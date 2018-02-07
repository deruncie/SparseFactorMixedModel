load_simulation_data = function(file = NULL){
  if(is.null(file)){
    if(file.exists('../setup.RData')) {
      load('../setup.RData')
      for(i in 1:10) names(setup) = sub('.','_',names(setup),fixed=T)
    } else{
      if (!requireNamespace("R.matlab", quietly = TRUE)) {
        stop("R.matlab needed to load setup.mat. Please install it.",
             call. = FALSE)
      }
      setup = readMat('../setup.mat')
      for(i in 1:10) names(setup) = sub('.','_',names(setup),fixed=T)
    }
  } else{
    load(file)
  }
  Y = setup$Y
  K = setup$A
  r = dim(K)[1]
  n = nrow(Y)
  if(dim(setup$X)[1] != n) setup$X = t(setup$X)
  data = data.frame(animal = gl(r,n/r))
  rownames(K) = data$animal
  if(is.null(colnames(Y))) colnames(Y) = paste('Trait',1:ncol(Y),sep='_')
  return(list(Y = Y, data = data, K_mats = list(animal = K),setup = setup))
}



#' Multiply matrices
#'
#' Multiplies two matrices. Like \code{%*%}, but always returns a base matrix, not a Matrix, even when
#' one or both of the two matrices are of class matrix.
#'
#' Uses RcppEigen when one or both matrices are of class dgCMatrix
#'
#' @param X1 matrix-like object
#' @param X2 matrix-like object
#'
#' @return matrix
#' @export
#'
#' @examples
`%**%` = function(X1,X2){
  if(is.null(X1)) return(X2)
  if(inherits(X1,'dgCMatrix') && inherits(X2,'matrix')) return(SxD(X1,X2))
  if(inherits(X1,'dgCMatrix') && inherits(X2,'dgCMatrix')) return(SxS(X1,X2))
  if(inherits(X1,'matrix') && inherits(X2,'matrix')) return(X1 %*% X2)
  if(inherits(X1,'data.frame') || inherits(X2,'data.frame')) return(as.matrix(X1) %**% as.matrix(X2))
  return(as.matrix(X1 %*% X2))
}


#' Modification to \code{image()}.
#'
#' Like `image,CHMfactor-method`, but with better defaults for correlation matrices especially.
#'
#' The default is to have the colors evenly spaced, instead of having the full range <0 and >0.
#'
#' @param x
#' @param zlim
#' @param breaks
#' @param colors
#'
#' @return
#' @export
#'
#' @examples
Image = function(X,dimnames=TRUE,...) {
  require(ggplot2)
  X = as.matrix(X)
  X[] = as.numeric(X)
  if(!dimnames) rownames(X) <- colnames(X) <- NULL
  if(length(unique(rownames(X))) < nrow(X)) rownames(X) <- NULL
  if(length(unique(colnames(X))) < ncol(X)) colnames(X) <- NULL

  X_tall = reshape2::melt(X)
  colnames(X_tall) = c('Var1','Var2','value')
  if(!is.null(colnames(X))) X_tall$Var2 = factor(X_tall$Var2,levels = colnames(X))
  if(!is.null(rownames(X))) X_tall$Var1 = factor(X_tall$Var1,levels = rev(rownames(X)))
  X_tall$value = as.numeric(X_tall$value)
  X_tall = X_tall[X_tall$value != 0,,drop=FALSE] # make it sparse again
  p <- ggplot(X_tall,aes(x=Var2,y=Var1,fill=value)) + geom_tile(alpha = 0.8) + xlab('') + ylab('') + theme_minimal() + expand_limits(fill=0)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  if(is.numeric(X_tall$Var2)) p <- p + xlim(0,ncol(X)+1)
  if(is.numeric(X_tall$Var1)) p <- p + ylim(nrow(X)+1,0)
  if(length(unique(X_tall$value))>2) p <- p + scale_fill_gradient2(na.value = "grey90",...)
  # print(p)
  p
}
Image2=function(x,zlim = NULL,breaks=20,colors = c('blue','white','red'),colorkey = TRUE,aspect=NULL,...){
  # if zlim not passed and the range of the data is outside of (-1,1), expands zlim range
  if(missing(zlim)){
    if(all(x[!is.na(x)]>0)) {
      zlim = c(0,max(x,na.rm=T))
    } else{
      zlim = c(-1,1)*max(abs(x),na.rm=T)
    }
  }
  zlim[is.na(zlim)] = 0
  if(missing(aspect)){
    if(ncol(x) > (2*nrow(x))) aspect = 1/1.5
    if(nrow(x) > (2*ncol(x))) aspect = 1.5
  }
  at = seq(zlim[1]-1e-10,zlim[2]+1e-10,length=breaks)
  colors = colorRampPalette(colors)(breaks)
  image(Matrix(x),at=at,col.regions=colors,colorkey=colorkey,aspect=aspect,...)
}

#' Checks if factors (columns of Lambda) can be safely dropped
#'
#' Factors are dropped if \code{prop} faction of the \eqn{\lambda_{ij}} elements are less than
#' \code{epsilon}. Checking happens stochastically at a decreasing rate during the chain controlled
#' by b0 and b1. If no factors can be dropped, then a new one is appended. The rate of adaptation
#' decreases through the chain, controlled by b0 and b1. Should work correctly over continuations of
#' previously stopped chains.
#'
#' @seealso \code{\link{sample_BSFG}}, \code{\link{plot.BSFG_state}}
# update_k = function( current_state, priors,run_parameters,data_matrices) {
update_k = function( BSFG_state) {
  data_matrices  = BSFG_state$data_matrices
  priors         = BSFG_state$priors
  run_parameters = BSFG_state$run_parameters
  current_state  = BSFG_state$current_state

  current_state_members = names(current_state)
  current_state = with(c(priors,run_parameters,data_matrices),within(current_state, {
    i = nrun
    n = nrow(F)
    k = ncol(Lambda)
    r = nrow(U_F)
    p = nrow(Lambda)
    b_F = ncol(X_F)

    # to determine if Lambda is kept, first remove empirical variance of F.
    Lambda_star = sweep(Lambda,2,sqrt(1/colMeans(F^2)),'*')
    lind = colMeans(abs(Lambda_star) < epsilon)    # proportion of elements in each column less than eps in magnitude
    vec = lind >= prop
    num = sum(vec)       # number of redundant columns

    if(num == 0 && all(lind < 0.995) && k < 2*p) { #add a column
      k             = k+1
      Lambda_prec   = cbind(Lambda_prec,rgamma(p,shape = Lambda_df/2, rate = Lambda_df/2))
      delta         = cbind(delta,rgamma(1,shape = delta_2_shape,rate = delta_2_rate))
      tauh          = matrix(cumprod(delta),nrow = 1)
      Plam          = sweep(Lambda_prec,2,tauh,'*')
      Lambda        = cbind(Lambda,rnorm(p,0,sqrt(1/Plam[,k])))
      F_h2_index    = c(F_h2_index,sample(1:ncol(h2s_matrix),1))
      F_h2          = h2s_matrix[,F_h2_index,drop=FALSE]
      tot_F_prec    = cbind(tot_F_prec,1)
      U_F           = cbind(U_F,rnorm(r,0,sqrt(sum(F_h2[,k]))))
      B_F           = cbind(B_F,0*rnorm(b_F,0,1))
      prec_B_F      = cbind(prec_B_F,c(tau_B_F))
      F             = cbind(F,rnorm(n,as.matrix(X_F %*% B_F[,k]) + as.matrix(Z %*% U_F[,k]),sqrt(1-sum(F_h2[,k]))))
    } else if(num > 0) { # drop redundant columns
      nonred = which(vec == 0) # non-redundant loadings columns
      while(length(nonred) < 2) {
        nonred = c(nonred,which(vec != 0)[1])
        vec[nonred[length(nonred)]] = 0
      }
      k = length(nonred)
      Lambda = Lambda[,nonred,drop=FALSE]
      Lambda_prec = Lambda_prec[,nonred,drop=FALSE]
      F = F[,nonred,drop=FALSE]
      for(red in which(vec == 1)){
        if(red == length(vec)) next
        # combine deltas so that the shrinkage of kept columns doesnt
        # decrease after dropping redundant columns
        delta[red+1] = delta[red+1]*delta[red]
      }
      if(is.null(dim(delta))) recover()
      delta = delta[,nonred,drop=FALSE]
      tauh = matrix(cumprod(delta),nrow=1)
      Plam = sweep(Lambda_prec,2,tauh,'*')
      F_h2 = F_h2[,nonred,drop=FALSE]
      F_h2_index = F_h2_index[nonred]
      tot_F_prec = tot_F_prec[,nonred,drop=FALSE]
      U_F = U_F[,nonred,drop=FALSE]
      B_F = B_F[,nonred,drop=FALSE]
      prec_B_F = prec_B_F[,nonred,drop=FALSE]
    }
  }))
  current_state = current_state[current_state_members]

  return(current_state)
}

get_meff = function(BSFG_state){
  get_k = function(n,sigma2,tau2,lambda2,s2=1){
    1/(1+n/sigma2 * tau2 * s2 * lambda2)
  }
  current_state = BSFG_state$current_state
  current_state = within(current_state, {
                         tauh = tauh / colMeans(F^2)
                         p = nrow(Lambda)
                         n = nrow(F)
                         })

  ki = with(current_state,get_k(n,1/t(tot_Eta_prec[rep(1,k),]),
             (Lambda_omega2[1]/tauh[rep(1,p),]),
             Lambda_phi2)
  )
  meff = colSums(1-ki)
  meff
}

#' Re-orders factors in decreasing order of magnitude
#'
#' @seealso \code{\link{sample_BSFG}}, \code{\link{plot.BSFG_state}}
reorder_factors = function(BSFG_state,factor_order = NULL){
  # re-orders factors in decreasing size of Lambda %*% F
  # based on current state
  # also re-orders Posterior

  current_state = BSFG_state$current_state

  # reorder factors based on var(Lambda) * var(F)
  Lambda = current_state$Lambda
  F = current_state$F

  # size is sum lambda_ij^2 * var(F_i)
  if(is.null(factor_order)) {
    sizes = colSums(Lambda^2) * colMeans(F^2)
    factor_order = order(sizes,decreasing=T)
  }

  reorder_params = c('Lambda','Lambda_prec','Plam',
                     'delta','tauh',
                     'F','B_F','U_F','F_h2','F_e_prec','tot_F_prec', 'B_F_prec'
  )

  # reorder currrent state
  for(param in reorder_params){
    if(! param %in% names(current_state)) next
    current_state[[param]] = current_state[[param]][,factor_order,drop=FALSE]
  }
  current_state$delta = matrix(c(current_state$tauh[1],current_state$tauh[-1]/current_state$tauh[-length(current_state$tauh)]),nrow=1)
  BSFG_state$current_state = current_state

  # reorder Posterior
  Posterior = BSFG_state$Posterior

  for(param in reorder_params){
    if(! param %in% names(Posterior)) next
    if(dim(Posterior[[param]])[1] == 0) next
    Posterior[[param]] = Posterior[[param]][,,factor_order,drop=FALSE]
  }

  BSFG_state$Posterior = Posterior

  return(BSFG_state)
}


rescale_factors_F = function(BSFG_state){
  # rescale factors based on F
  BSFG_state$current_state = within(BSFG_state$current_state,{
    F_sizes = colMeans(F^2)
    F = sweep(F,2,sqrt(F_sizes),'/')
    B_F = sweep(B_F,2,sqrt(F_sizes),'/')
    U_F = sweep(U_F,2,sqrt(F_sizes),'/')
    B_F_prec = sweep(B_F_prec,2,F_sizes,'*')
    Lambda = sweep(Lambda,2,sqrt(F_sizes),'*')
    delta_factor = c(F_sizes[1],exp(diff(log(F_sizes))))
    delta[] = delta / delta_factor
    tauh[] = cumprod(delta)
    Plam[] = sweep(Lambda_prec,2,tauh,'*')
  })
  return(BSFG_state)
}

#' Saves current state in Posterior
#'
#' Saves current state in Posterior
#' @seealso \code{\link{sample_BSFG}}, \code{\link{plot.BSFG_state}}
save_posterior_sample = function(BSFG_state) {
  # All parameters in current are matrices.
  # Posterior is a list of arrays.
  # All factor parameters are matrices with ncol == k
  # Posterior arrays are expanded / contracted as the number of factors changes (with update_k)

  current_state = BSFG_state$current_state
  Posterior = BSFG_state$Posterior

  total_samples = Posterior$total_samples + 1
  sp_num = Posterior$sp_num + 1
  Posterior$total_samples = total_samples
  Posterior$sp_num = sp_num

  current_state = within(current_state,{
    # re-transform random effects using RE_L (RE_L %*% diag(D) %*% t(RE_L) = bdiag(K_mats))
    U_R = as.matrix(BSFG_state$data_matrices$RE_L %*% U_R)
    U_F = as.matrix(BSFG_state$data_matrices$RE_L %*% U_F)
    # transform variables so that the variance of each column of F is 1.
    F_var = 1/tot_F_prec
    U_F = sweep(U_F,2,sqrt(F_var),'/')
    B_F = sweep(B_F,2,sqrt(F_var),'/')
    F = sweep(F,2,sqrt(F_var),'/')
    Lambda = sweep(Lambda,2,sqrt(F_var),'*')
    tauh[] = tauh * tot_F_prec
  })

  sp = dim(Posterior$Lambda)[1]

  for(param in Posterior$posteriorSample_params){
    ncol_current = ncol(current_state[[param]])
    ncol_Posterior = dim(Posterior[[param]])[3]
    if(ncol_current > ncol_Posterior){
      Posterior[[param]] = abind(Posterior[[param]],array(0,dim = c(sp,nrow(current_state[[param]]),ncol_current-ncol_Posterior)),along=3)
    }
    if(ncol_current < ncol_Posterior) {
      Posterior[[param]] = Posterior[[param]][,,1:ncol_current,drop=F]
    }
    Posterior[[param]][sp_num,,] = current_state[[param]]
  }

  for(param in Posterior$posteriorMean_params){
    Posterior[[param]] = (Posterior[[param]]*(total_samples - 1) + current_state[[param]])/total_samples
  }

  return(Posterior)
}

reset_Posterior = function(Posterior,BSFG_state){
  current_state = BSFG_state$current_state

  # re-transform random effects using RE_L
  current_state$U_F = BSFG_state$data_matrices$RE_L %*% current_state$U_F
  current_state$U_R = BSFG_state$data_matrices$RE_L %*% current_state$U_R

  for(param in Posterior$posteriorSample_params){
    Posterior[[param]] = array(0,dim = c(0,dim(current_state[[param]])))
    dimnames(Posterior[[param]])[2:3] = dimnames(current_state[[param]])
  }
  for(param in Posterior$posteriorMean_params) {
    if(Posterior$total_samples == 0) {
      Posterior[[param]] = array(0,dim = dim(current_state[[param]]))
      dimnames(Posterior[[param]]) = dimnames(current_state[[param]])
    }
  }
  Posterior$sp_num = 0
  Posterior
}

expand_Posterior = function(Posterior,size){
  for(param in Posterior$posteriorSample_params){
    Posterior[[param]] = abind(Posterior[[param]],array(NA,dim = c(size,dim(Posterior[[param]])[2:3])),along = 1)
  }
  Posterior
}

#' Resets Posterior samples
#'
#' Clears and resets the saved Posterior samples. Updates the burn parameter of
#'     \code{run_parameters} to reflect all previous samples in the chain now count as burnin
#' @seealso \code{\link{sample_BSFG}}, \code{\link{plot.BSFG_state}}
clear_Posterior = function(BSFG_state) {
  # resets Posterior samples if burnin was not sufficient
  Posterior = BSFG_state$Posterior
  run_parameters = BSFG_state$run_parameters

  run_parameters$burn = run_parameters$burn + run_parameters$thin*Posterior$total_samples

  Posterior$total_samples = 0
  Posterior = reset_Posterior(Posterior,BSFG_state)

  if(length(list.files(path = Posterior$folder))>0) system(sprintf('rm %s/*',Posterior$folder))
  Posterior$files = c()

  BSFG_state$Posterior = Posterior
  BSFG_state$run_parameters = run_parameters
  return(BSFG_state)
}


#' Saves a chunk of posterior samples
#'
#' Saves a chunk of posterior samples, each paramter in its own RData file.
#' The complete chain for a particular parameter can be re-loaded with \link{load_posterior_param}
#'
#' @param Posterior a Posterior list from a BSFG_state object
#' @param folder the folder to save the RData files
#' @return Posterior
save_posterior_chunk = function(BSFG_state){
  Posterior = BSFG_state$Posterior
  folder = Posterior$folder
  if(!dir.exists(folder)) dir.create(folder)
  file_suffix = sprintf('%d.rds',Posterior$total_samples)
  res = sapply(c(Posterior$posteriorSample_params,Posterior$posteriorMean_params),function(param) {
    file_name = sprintf('%s/%s_%s',folder,param,file_suffix)
    samples = Posterior[[param]]
    if(length(samples) > 0) {
      saveRDS(samples,file = file_name,compress = FALSE)
    }
  })
  print(sprintf('%d files',length(grep(file_suffix,list.files(path = folder)))))
  if(length(grep(file_suffix,list.files(path = folder)))>0) {
    Posterior$files = unique(c(Posterior$files,file_suffix))
  }
  Posterior = reset_Posterior(Posterior,BSFG_state)
  BSFG_state$Posterior = Posterior
  saveRDS(Posterior,file = sprintf('%s/Posterior_base.rds',folder))
  return(BSFG_state)
}

#' load the posterior samples of a single parameter from all saved chunks
#'
#' @param folder folder to find Posterior chunk files
#' @param param Name of parameter to load
#' @param samples vector of sample indices to load. If NULL, all samples loaded
#' @return array with all chuncks of Posterior samples appended together
load_posterior_param = function(BSFG_state,param,chunks = NULL){
  if(length(BSFG_state$Posterior$files) == 0) return(c())
  if(grepl('RData',BSFG_state$Posterior$files[1])) {
    return(load_posterior_param_old(BSFG_state,param,chunks))
  }
  folder = BSFG_state$Posterior$folder
  param_files = paste0(folder,'/',param,'_',BSFG_state$Posterior$files)
  all_files = list.files(path=folder,full.names = T)
  param_files = param_files[param_files %in% all_files]
  n_files = length(param_files)
  if(is.null(chunks)) chunks = 1:n_files
  param_files = na.omit(param_files[chunks])
  if(length(param_files) == 0) return(c())

  samples = readRDS(param_files[1])
  samples_dim = dim(samples)
  current_row = 0
  if(param %in% BSFG_state$Posterior$posteriorMean_params) {
    all_samples = samples / length(param_files)
  } else{
    all_samples = array(0,dim = c(BSFG_state$Posterior$total_samples*length(chunks)/n_files,samples_dim[2:3]))
    dimnames(all_samples) = dimnames(samples)
    all_samples[1:samples_dim[1],,] = samples
    current_row = samples_dim[1]
    if(current_row == 0) return(c())
  }

  # load other files
  if(length(param_files) > 1) {
    for(i in 2:length(param_files)){
      samples = readRDS(param_files[i])
      if(length(samples_dim) == 2) {
        all_samples = all_samples + samples / length(param_files)
      } else{
        samples_dim = dim(samples)
        if(current_row+samples_dim[1] > dim(all_samples)[1]){
          all_samples = abind(all_samples,array(0,dim = c(current_row+samples_dim[1] - dim(all_samples)[1],dim(all_samples)[2:3])),along = 1)
        }
        if(samples_dim[2] > dim(all_samples)[2]){
          all_samples = abind(all_samples,array(0,dim = c(samples_dim[1],samples_dim[2]-dim(all_samples)[2],dim(all_samples)[3])),along = 2)
        }
        if(samples_dim[3] > dim(all_samples)[3]){
          all_samples = abind(all_samples,array(0,dim = c(dim(all_samples)[1:2],samples_dim[3]-dim(all_samples)[3])),along = 3)
        }
        all_samples[current_row+1:dim(samples)[1],1:samples_dim[2],1:samples_dim[3]] = samples
        current_row = current_row+dim(samples)[1]
      }
    }
  }

  return(all_samples)
}

#' load the posterior samples of a single parameter from all saved chunks
#'
#' @param folder folder to find Posterior chunk files
#' @param param Name of parameter to load
#' @return array with all chuncks of Posterior samples appended together
load_posterior_param_old = function(BSFG_state,param,chunks=NULL){
  if(length(BSFG_state$Posterior$files) == 0) return(c())
  folder = BSFG_state$Posterior$folder
  n_files = length(BSFG_state$Posterior$files)
  param_files = paste0(folder,'/',param,'_',BSFG_state$Posterior$files)
  all_files = list.files(path=folder,full.names = T)
  param_files = param_files[param_files %in% all_files]
  n_files = length(param_files)
  if(is.null(chunks)) chunks = 1:n_files
  param_files = na.omit(param_files[chunks])
  if(length(param_files) == 0) return(c())

  load(param_files[1])
  samples_dim = dim(samples)
  current_row = 0
  if(param %in% BSFG_state$Posterior$posteriorMean_params) {
    all_samples = samples / length(param_files)
  } else{
    all_samples = array(0,dim = c(BSFG_state$Posterior$total_samples*length(chunks)/n_files,samples_dim[2:3]))
    all_samples[1:samples_dim[1],,] = samples
    current_row = samples_dim[1]
    if(current_row == 0) return(c())
  }

  # load other files
  if(length(param_files) > 1) {
    for(i in 2:length(param_files)){
      load(param_files[i])
      if(length(samples_dim) == 2) {
        all_samples = all_samples + samples / length(param_files)
      } else{
        samples_dim = dim(samples)
        if(current_row+samples_dim[1] > dim(all_samples)[1]){
          all_samples = abind(all_samples,array(0,dim = c(current_row+samples_dim[1] - dim(all_samples)[1],dim(all_samples)[2:3])),along = 1)
        }
        if(samples_dim[2] > dim(all_samples)[2]){
          all_samples = abind(all_samples,array(0,dim = c(samples_dim[1],samples_dim[2]-dim(all_samples)[2],dim(all_samples)[3])),along = 2)
        }
        if(samples_dim[3] > dim(all_samples)[3]){
          all_samples = abind(all_samples,array(0,dim = c(dim(all_samples)[1:2],samples_dim[3]-dim(all_samples)[3])),along = 3)
        }
        all_samples[current_row+1:dim(samples)[1],1:samples_dim[2],1:samples_dim[3]] = samples
        current_row = current_row+dim(samples)[1]
      }
    }
  }

  return(all_samples)
}

#' Re-loads a full Posterior list with all parameters
#'
#' @param BSFG_state a BSFG_state object (with empty Posterior)
#' @param params list of parameters to load. If NULL, all parameters will be loaded
#' @return Posterior list, as part of a BSFG_state object
reload_Posterior = function(BSFG_state,params = NULL){
  folder = BSFG_state$Posterior$folder
  Posterior = readRDS(paste(folder,'Posterior_base.rds',sep='/'))
  BSFG_state$Posterior = Posterior
  if(is.null(params)) params = c(Posterior$posteriorSample_params,Posterior$posteriorMean_params)
  for(param in params){
    Posterior[[param]] = load_posterior_param(BSFG_state,param)
    try(dimnames(Posterior[[param]]) <- dimnames(BSFG_state$Posterior[[param]]),silent=T)
  }
  Posterior
}


# pull out specific parameters from a specific sample from Posterior. Only returns
# parameters in listed in \code{terms}, and only if they have posterior samples (not posterior means)
make_current_state = function(Posterior,sample,terms){
  sample_terms = terms[terms %in% Posterior$posteriorSample_params]
  current_state = lapply(sample_terms,function(x) array(Posterior[[x]][sample,,],dim = dim(Posterior[[x]])[-1],dimnames = dimnames(Posterior[[x]])[-1]))
  names(current_state) = sample_terms
  current_state
}

#' Calculates the posterior mean of a function of parameters
#'
#' This function will apply the supplied function to each posterior sample of the chain. Variables
#'    referenced in FUN will be selected from the following search locations (in this order):
#'    1) Posterior$posteriorSample_params
#'    2) data_matrices
#'    3) priors
#'    4) Posterior$posteriorMean_params
#'    5) current_state
#'    6) calling environment (ex. sapply)
#'    7) global environment
#'
#' If mc.cores > 1, the operation will be parallelized. To reduce memory requirements,
#' if mc.cores > 1, a PSOCK cluster is created. This cluster only copyies in data when it
#' is actually used. The code is written so that only necessary objects are used by this cluster,
#' not the whole user's environment. Inside this cluster, mclapply is run, which forks the process inside the PSOCK
#' cluster. Creating the cluster takes some time, so often mc.cores=1 is faster. If the operation
#' is large and Posterior has many samples, this can still create memory issues, so limiting mc.cores
#' might still be necessary
#'
#' @param BSFG_state A BSFG_state object including a re-loaded Posterior list
#' @param FUN Operations to be applied to each posterior sample. Write as if this were operating
#'     within current_state. Can use priors, data_matrices, and other elements of current_state
#' @param samples (optional) vector of sample indexes to use in the computation
#' @param mc.cores (optional) number of cores to use for computations. See note about memory requirements.
#'
#' @return array of n_samples x dim1 x dim2 where dim1 and dim2 are the dimensions of the calculated
#'     parameter per posterior sample
get_posterior_FUN = function(BSFG_state,FUN,samples = NULL,mc.cores = 1) {
  FUN = match.call()[[3]]
  if(is(FUN,'character')){
    FUN = parse(text=FUN)
  }
  terms = all.vars(FUN)
  extra_terms = terms[terms %in% with(BSFG_state,c(names(current_state),names(data_matrices),names(priors),names(Posterior))) == F]
  extra_env = c()
  for(term in extra_terms){
    if(term %in% ls(parent.frame(2))) {
      extra_env[[term]] = parent.frame(2)[[term]]
    }
  }
  if(is.null(samples)) {  # count # available samples for the first term in terms (assuming it is in Posterior)
    term1 = terms[terms %in% BSFG_state$Posterior$posteriorSample_params][1]
    if(!term1 %in% names(BSFG_state$Posterior)) {
      print(sprintf('%s not in Posterior',term1))
      return(NULL)
    }
    samples = 1:dim(BSFG_state$Posterior[[term1]])[1]
  }
  base_env = with(BSFG_state,c(data_matrices,priors,Posterior[BSFG_state$Posterior$posteriorMean_params],current_state))
  base_env = c(base_env,extra_env)
  Posterior = BSFG_state$Posterior
  per_sample_fun = function(sample_index_i) {
    # get current sample of each of the terms in FUN
    current_sample = make_current_state(Posterior,sample_index_i,terms)
    # evaluate FUN in an environment constructed from current_sample, and BSFG_state, taking current_sample first
    env = c(current_sample,base_env)
    result = eval(FUN,envir = env)
    if(is(result,'Matrix')) result = as.matrix(result)
    result
  }
  sample_1_result <- tryCatch(per_sample_fun(1),
                              error = function(e) {
                                message(e)
                                return(NULL)
                              })
  if(is.null(sample_1_result)) return(NULL)
  dim_1 = dim(sample_1_result) # get the dimension of the returned value
  # calculate value for each sample

  if(mc.cores == 1) {
    res = do.call('c',lapply(samples,per_sample_fun))
  } else{
    cluster = parallel::makeCluster(1)
    doParallel::registerDoParallel(cluster)
    res = foreach::foreach(i=1,.export = c('make_current_state'),.combine = 'c') %dopar% {
      # This way, we make a new process, and only pass it the data needed for per_sample_fun
      do.call('c',parallel::mclapply(samples,per_sample_fun,mc.cores = mc.cores))
    }
    parallel::stopCluster(cluster)
    # res = do.call(c,mclapply(samples,per_sample_fun,mc.cores = mc.cores))
  }
  # re-formulate into an appropriate array with the first dimension as samples
  if(is.null(dim_1)) {
    dim_1 = length(sample_1_result)
    res = matrix(res,ncol = dim_1,byrow=T)
    colnames(res) = names(sample_1_result)
  } else {
    res = aperm(array(res,dim = c(dim_1,length(samples))),c(3,1,2))
    dimnames(res)[2:3] = dimnames(sample_1_result)
  }
  res
}


#' Calculates posterior mean of a function of parameters
#'
#' The synthetic parameter can be pre-calculated using \link{get_posterior_FUN}, or provided directly
#'     as an array with the first dimension the samples.
#'
#' @param X either an array of posterior samples (either a parameter from \code{Posterior}, or
#'     an object generated by \link{get_posterior_fun}), or the BSFG_state object with re-loaded Posterior
#' @param FUN (optional) if \code{X} is a BSFG_state object, the function to calculate the synthetic parameter
#' @param bychunk (optional) if \code{TRUE}, will re-load each chunk of the Posterior separately
#'     and calculate the posterior mean of \code{FUN} separately for each chunk, and then return
#'     the overall posterior mean. This saves memory because the whole posterior does not need to be
#'     loaded into memory at once. Only necessary terms from Posterior are loaded.
#'
#' @return posterior mean matrix
get_posterior_mean = function(X,FUN,bychunk = FALSE,mc.cores = 1,...){
  result = NULL
  if(!bychunk) {
    if(is(X,'BSFG_state')) {
      BSFG_state = X
      FUN = match.call()$FUN
      X = do.call(get_posterior_FUN,list(BSFG_state=BSFG_state,FUN=FUN,mc.cores=mc.cores))
    }
    if(length(dim(X)) == 3) {
      result = matrix(colMeans(matrix(X,nr = dim(X)[1])),nr = dim(X)[2])
      dimnames(result) = dimnames(X)[-1]
    }
    if(length(dim(X)) == 2) {
      result = colMeans(X)
      names(result) = colnames(X)
    }
  } else{
    if(!is(X,'BSFG_state')) stop('Provide a BSFG_state object as "X"')
    BSFG_state = X
    FUN = match.call()$FUN
    if(is(FUN,'character')){
      FUN = parse(text=FUN)
    }
    terms = all.vars(FUN)
    terms = terms[terms %in% BSFG_state$Posterior$posteriorSample_params]
    n_files = length(BSFG_state$Posterior$files)
    result = 0
    chunk=1
#
    pb = txtProgressBar(min=0,max = n_files,style=3)
    for(chunk in 1:n_files){
      for(term in terms){
        BSFG_state$Posterior[[term]] = load_posterior_param(BSFG_state,term,chunks = chunk)
      }
      if(length(BSFG_state$Posterior[[term]]) > 0) {
        samples = do.call(get_posterior_FUN,list(BSFG_state=BSFG_state,FUN=FUN,mc.cores=mc.cores))
        result = result + dim(samples)[1]*get_posterior_mean(samples,bychunk = FALSE)
        rm(samples)
        gc()
      }
      setTxtProgressBar(pb, chunk)
    }
    close(pb)
    result = result / BSFG_state$Posterior$total_samples
  }
  result
}



#' Calculates Highest Posteriod Density intervals of a function of parameters
#'
#' The synthetic parameter can be pre-calculated using \link{get_posterior_FUN}, or provided directly
#'     as an array with the first dimension the samples.#'
#' @param X either an array of posterior samples (either a parameter from \code{Posterior}, or
#'     an object generated by \link{get_posterior_FUN}), or the BSFG_state object with re-loaded Posterior
#' @param FUN (optional) if \code{X} is a BSFG_state object, the function to calculate the synthetic parameter
#' @param prob A numeric scalar in the interval (0,1) giving the target probability content of the intervals.
#'     The nominal probability content of the intervals is the multiple of 1/nrow(obj) nearest to prob
#' @param ... other parameters passed to \link{get_posterior_FUN}
#'
#' @return array with first dimension of length=2 giving the lower and upper limits of the parameters
#'     in the other two dimensions.
get_posterior_HPDinterval = function(X,FUN = NULL,prob = 0.95,...){
  if(is(X,'BSFG_state')) {
    BSFG_state = X
    FUN = match.call()[[3]]
    X = do.call(get_posterior_FUN,list(BSFG_state=BSFG_state,FUN=FUN))
  }
  dims = dim(X)
  if(length(dims) == 3) result = aperm(array(apply(X,3,function(x) HPDinterval(mcmc(x),prob=prob)),c(dims[2],2,dims[3])),c(2,1,3))
  if(length(dims) == 2) result = t(HPDinterval(mcmc(X),prob=prob))
  result
}


#' Converts dense Matrix to "matrix" faster!
#'
#' @param X a ddenseMatrix
#'
#' @return a matrix
#' @export
#'
#' @examples
toDense = function(X) {
  matrix(X@x,nrow(X))
}
