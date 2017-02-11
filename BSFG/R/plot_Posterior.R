set_device = function(device) {
  # create new devices up to device
  devices = dev.list()
  devices = devices[names(devices) != 'RStudioGD']
  while(length(devices) < device){
    dev.new(noRStudioGD = TRUE)
    devices = dev.list()
    devices = devices[names(devices) != 'RStudioGD']
  }
  dev.set(devices[device])
}

trace_plot = function(data,main = NULL,ylim = NULL){
  if(is.null(ylim)) {
    range_y = range(data)
    max_range = max(abs(range_y))
    ylim = c(-max_range,max_range)
  }
  plot(NA,NA,xlim = c(0,nrow(data)),ylim = ylim,main= main,xlab = "iteration")
  for(i in 1:ncol(data)) lines(data[,i],col=i)
}

trace_plot_h2s = function(Posterior, n_factors = 8, device = NULL){
  if(!is.null(device)) {
    set_device(device)
  }

  F_h2_samples = Posterior$F_h2

  rows = ceiling(n_factors / 2)
  cols = 2*ceiling(n_factors / rows)
  par(mfrow=c(rows,cols))

  for(k in 1:min(n_factors,dim(F_h2_samples)[3])){
    trace_plot(as.matrix(F_h2_samples[,,k,drop=F]),main = sprintf('Factor %d h2s',k),ylim = c(0,1))
    hist(F_h2_samples[,1,k],breaks=seq(0,1,length=100),xlim = c(-0.1,1),main = sprintf('Factor %d',k))
    if(dim(F_h2_samples)[2] > 1) {
      for(i in 2:dim(F_h2_samples)[2]){
        hist(F_h2_samples[,1,k],breaks=seq(0,1,length=100),col = i)
      }
    }
  }
}

trace_plot_Lambda = function(Posterior, n_factors = 16, device = NULL){
  if(!is.null(device)) {
    set_device(device)
  }

  Lambda = Posterior$Lambda

  rows = ceiling(sqrt(n_factors))
  cols = ceiling(n_factors / rows)
  par(mfrow=c(rows,cols))

  for(k in 1:min(n_factors,dim(Lambda)[3])){
    o = order(-abs(colMeans(Lambda[,,k])))
    traces = Lambda[,o[1:5],k]
    trace_plot(traces,main = sprintf('Factor %d lambdas',k))
  }
}

plot_factor_correlations = function(Lambda, sim_Lambda){
  cors = abs(cor(Lambda,sim_Lambda))
  cors = rbind(cors,apply(cors,2,max))
  cors = cbind(cors,apply(cors,1,max))
  image(t(cors)[,nrow(cors):1],xlab = 'Actual factors', ylab = 'Estimated factors',xaxt='n',yaxt='n', main = 'Correlation of fitted and\nsimulated factor loadings'
  )
  axis(1,at=seq(0,1,length=ncol(cors)),labels = c(1:(ncol(cors)-1),'best'),las=2)
  axis(2,at=seq(0,1,length=nrow(cors)),labels = c('best',(nrow(cors)-1):1),las=2)
}

plot_element_wise_correlations = function(actual_G, estimated_G, main = NULL){
  xlim = ylim = range(c(actual_G,estimated_G))
  plot(c(actual_G),c(estimated_G),xlim = xlim,ylim=ylim,main = main,xlab = 'actual',ylab = 'estimated')
  abline(0,1)
}

plot_current_state_simulation = function(BSFG_state, device = NULL){
  require(cowplot)
  if(!is.null(device)) {
    set_device(device)
  }
  par(mfrow=c(3,3))

  setup = BSFG_state$setup
  run_parameters = BSFG_state$run_parameters
  run_variables = BSFG_state$run_variables

  current_state = within(BSFG_state$current_state,{
    # transform variables so that the variance of each column of F is 1.
    F_var = 1/tot_F_prec
    F_a = sweep(F_a,2,sqrt(F_var),'/')
    F = sweep(F,2,sqrt(F_var),'/')
    Lambda = sweep(Lambda,2,sqrt(F_var),'*')
  })

  Lambda = current_state$Lambda
  F_h2 = current_state$F_h2
  if(is.null(dim(F_h2))) F_h2 = matrix(F_h2,nrow=1)
  E_a_prec = current_state$tot_Y_prec / current_state$resid_h2
  resid_Y_prec = current_state$tot_Y_prec / (1-current_state$resid_h2)
  p = run_variables$p


  # correlation of factors
  sim_Lambda = setup$error_factor_Lambda
  plot_factor_correlations(Lambda, sim_Lambda)

  # element-wise correlations
  G_plots = lapply(1:nrow(F_h2),function(re) {
    G_est = tcrossprod(sweep(Lambda,2,sqrt(F_h2[re,]),'*')) + diag(c(current_state$resid_h2[re,]/current_state$tot_Y_prec))
    G_act = setup$G
    RE_name = rownames(F_h2)[re]
    if(is.null(RE_name)) RE_name = 'Va'
    plot_element_wise_correlations(G_act, G_est, main = sprintf('G: %s elements',RE_name))
  })
  E_est = tcrossprod(sweep(Lambda,2,sqrt(1-colSums(F_h2)),'*')) + diag(c((1-colSums(current_state$resid_h2))/current_state$tot_Y_prec))
  E_act = setup$R
  plot_element_wise_correlations(E_act, E_est,main = 'E elements')

  # factor h2s
  plot_factor_h2s(F_h2)

  # B's
  if(dim(setup$B)[1] > 1) {
    B = BSFG_state$current_state$B
    B_factor = BSFG_state$current_state$B_F   %*% t(BSFG_state$current_state$Lambda)
    plot(c(B),c(setup$B))
    abline(0,1)
    plot(c(B_factor),c(setup$B_F %*% t(setup$error_factor_Lambda)))
    abline(0,1)
    xlim = ylim = range(c(B[-1,],B_factor))
    plot(c(B_factor),c(B[-1,]),xlim = xlim,ylim=ylim);abline(0,1)
  }
}

plot_fixed_effects = function(B_act, B_resid,B_factor,B_total){
  xlim = ylim = range(unlist(c(B_act, B_resid,B_factor,B_total)))
  plot(c(B_act),c(B_resid),xlim=xlim,ylim=ylim,col=1,main = 'Fixed effects')
  points(c(B_act),c(B_factor),col=2)
  points(c(B_act),c(B_total),pch=21,col=0,bg=2)
  abline(0,1)
}

calc_posterior_mean_cov = function(Posterior,random_effect){
  with(Posterior,{
    p = dim(Lambda)[2]
    G = matrix(0,p,p)
    for(i in 1:sp_num){
      if(random_effect == 'Ve'){
        factor_h2s_i = 1-colSums(F_h2[i,,,drop=FALSE])
        resid_h2s_i = 1-colSums(resid_h2[i,,,drop=FALSE])
      } else{
        factor_h2s_i = F_h2[i,random_effect,]
        resid_h2s_i = resid_h2[i,random_effect,]
      }
      G_i = tcrossprod(sweep(Lambda[i,,],2,sqrt(factor_h2s_i),'*')) + diag(c(resid_h2s_i/tot_Y_prec[i,1,]))
      G = G + G_i/sp_num
    }
    G
  })
}

calc_posterior_mean_Lambda = function(Posterior){
  return(with(Posterior,apply(Lambda,c(2,3),mean)))
}

plot_factor_h2s = function(F_h2) {
  require(reshape2)
  if(is.null(rownames(F_h2))) rownames(F_h2) = 'Va'
  F_h2 = rbind(F_h2,1-colSums(F_h2))
  rownames(F_h2)[nrow(F_h2)] = 'Ve'
  colnames(F_h2) = 1:ncol(F_h2)
  F_h2_data = melt(F_h2)
  barplot(F_h2,main = 'Factor h2s',legend.text = rownames(F_h2))
}

plot_HPDIntervals_factor_h2 = function(Posterior){
  require(MCMCpack)
}

plot_posterior_simulation = function(BSFG_state, device = NULL){
  require(cowplot)
  if(!is.null(device)) {
    set_device(device)
  }
  par(mfrow=c(3,3))

  setup = BSFG_state$setup
  Posterior = BSFG_state$Posterior
  p = dim(Posterior$Lambda)[1]
  plot_element_wise_correlations(setup$R, calc_posterior_mean_cov(Posterior,'Ve'),main = 'E elements')
  for(RE_name in dimnames(Posterior$F_h2)[[2]]) {
    plot_element_wise_correlations(setup$G, calc_posterior_mean_cov(Posterior,RE_name),main = sprintf('G: %s elements',RE_name))
  }

  plot_factor_correlations(calc_posterior_mean_Lambda(Posterior),setup$error_factor_Lambda)

  plot_factor_h2s(apply(Posterior$F_h2,c(2,3),mean))


  # B's # These seem to be wrong.

  if(dim(setup$B)[1] > 1) {
    B_mean = apply(Posterior$B,c(2,3),mean)
    B_factor_mean = with(c(Posterior,BSFG_state$data_matrices), {
      if(ncol(X_F) == 0) return(rep(0,dim(Lambda)[1]))
      matrix(rowMeans(sapply(1:sp_num,function(i) B_F[i,,] %*% t(Lambda[i,,]))),nrow = ncol(X_F))
    })
    plot(c(B_mean),c(setup$B))
    abline(0,1)
    plot(c(B_factor_mean),c(setup$B_F %*% t(setup$error_factor_Lambda)))
    abline(0,1)
    xlim = ylim = range(c(B_mean[-1,],B_factor_mean))
    plot(c(B_factor_mean),c(B_mean[-1,]),xlim = xlim,ylim=ylim);abline(0,1)
    B_f_HPD = HPDinterval(mcmc(Posterior$B_F[,1,]))
    B_f_mean = colMeans(Posterior$B_F[,1,])

    plot(1:length(B_f_mean),B_f_mean,xlim = c(1,length(B_f_mean)),ylim = range(B_f_HPD),xlab = '',main = 'Posterior B_F')
    arrows(seq_along(B_f_mean),B_f_HPD[,1],seq_along(B_f_mean),B_f_HPD[,2],length=0)
    abline(h=0)
  }
}

plot_diagnostics_simulation = function(BSFG_state){
  plot_current_state_simulation(BSFG_state, device = 2)
  if(BSFG_state$Posterior$sp_num > 0) {
    plot_posterior_simulation(BSFG_state, device = 3)
    trace_plot_h2s(BSFG_state$Posterior,device = 4)
    trace_plot_Lambda(BSFG_state$Posterior,device = 5)
  }
}
