my_detectCores = function() {
  ncores = suppressWarnings(as.numeric(system('printenv SLURM_CPUS_PER_TASK',intern=T)))
  if(length(ncores) == 0 || is.na(ncores)) ncores = parallel::detectCores()
  ncores
}


make_model_setup = function(formula,data,relmat = NULL) {
  # ensure data has rownames
  if(is.null(rownames(data))) rownames(data) = 1:nrow(data)

  # check that RE's have approporiate levels in data
  RE_levels = list() # a list of levels for each of the random effects
  if(is.null(relmat)) relmat = list()
  for(re in names(relmat)) {
    # check that K is a matrix, then convert to Matrix
    if(is.list(relmat[[re]])){
      if(is.data.frame(relmat[[re]]$K)) relmat[[re]]$K = as.matrix(relmat[[re]]$K)
      if(is.matrix(relmat[[re]]$K)) relmat[[re]]$K = Matrix(relmat[[re]]$K,sparse=T)
      if(is.null(rownames(relmat[[re]]$K))) stop(sprintf('K %s must have rownames',re))
      RE_levels[[re]] = rownames(relmat[[re]]$K)
    } else{
      if(is.data.frame(relmat[[re]])) relmat[[re]] = as.matrix(relmat[[re]])
      if(is.matrix(relmat[[re]])) relmat[[re]] = Matrix(relmat[[re]],sparse=T)
      if(is.null(rownames(relmat[[re]]))) stop(sprintf('K %s must have rownames',re))
      RE_levels[[re]] = rownames(relmat[[re]])
    }
  }
  for(re in names(RE_levels)){
    if(!re %in% colnames(data)) stop(sprintf('Column "%s" required in data',re))
    data[[re]] = as.factor(data[[re]]) # ensure 'data[[re]]' is a factor
    if(!all(data[[re]] %in% RE_levels[[re]])) stop(sprintf('Levels of random effect %s missing.',re))
    data[[re]] = factor(data[[re]],levels = RE_levels[[re]]) # add levels to data[[re]]
  }

  # Use lme4 to evaluate formula in data
  # ensure there is a response
  response = 'y'
  while(response %in% all.vars(formula)){
    response = paste0(response,response)
  }
  formula = update(formula,sprintf('%s~.',response))
  data[[response]] = 1
  lmod <- lme4::lFormula(formula,data=data,#weights=weights,
                         control = lme4::lmerControl(check.nobs.vs.nlev = 'ignore',check.nobs.vs.nRE = 'ignore'))

  # compute RE_setup
  RE_terms = lmod$reTrms

  # construct the RE_setup list
  # contains:
  # Z: n x r design matrix
  # K: r x r PSD covariance matrix
  RE_setup = list()
  for(i in 1:length(RE_terms$cnms)){
    term = names(RE_terms$cnms)[i]
    n_factors = length(RE_terms$cnms[[i]])  # number of factors for this grouping factor

    # extract combined Z matrix
    combined_Zt = RE_terms$Ztlist[[i]]
    Zs_term = tapply(1:nrow(combined_Zt),gl(n_factors,1,nrow(combined_Zt),labels = RE_terms$cnms[[i]]),function(x) Matrix::t(combined_Zt[x,,drop=FALSE]))

    # extract K from relmat. If missing, assign to NULL
    K = NULL
    if(term %in% names(relmat)) {
      K = relmat[[term]]
      if(is(K,'dsCMatrix')) K = as(K,'dgCMatrix')
      if(nnzero(K)/length(K) > .5) K = as.matrix(K)  # if not really sparse
    }

    if(!is.null(K)) {
      if(!all(colnames(Zs_term[[1]]) %in% rownames(K))) stop('rownames of K not lining up with Z')
      K = K[colnames(Zs_term[[1]]),colnames(Zs_term[[1]])]
    } else {
      K = as(diag(1,ncol(Zs_term[[1]])),'dgCMatrix')
      rownames(K) = colnames(K) = colnames(Zs_term[[1]])
    }


    # make an entry in RE_setup for each random effect
    for(j in 1:n_factors){
      # name of variance component
      name = term
      if(n_factors > 1) name = paste(name,RE_terms$cnms[[i]][[j]],sep='.')
      while(name %in% names(RE_setup)) name = paste0(name,'.1') # hack for when same RE used multiple times

      # Z matrix
      Z = as(Zs_term[[j]],'dgCMatrix')


      RE_setup[[name]] = list(
        term = term,
        Z = Z,
        K = K
      )
    }

    # add names to RE_setup if needed
    n_RE = length(RE_setup)
    for(i in 1:n_RE){
      if(is.null(names(RE_setup)[i]) || names(RE_setup)[i] == ''){
        names(RE_setup)[i] = paste0('RE.',i)
      }
    }

  }
  return(list(lmod = lmod, RE_setup = RE_setup))
}
