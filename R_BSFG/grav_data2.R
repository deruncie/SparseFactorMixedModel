library(microbenchmark)
library(MCMCpack)
library(BSFG)
library(ggplot2)
library(Rcpp)
library(qtl)
library(qtlcharts)
library(reshape)
library(splines)

# create a folder for holding the posterior samples of the current chain (multiple folders could be used for different chains)
rep = "4"
folder = sprintf('grav_spline_rep_%s',rep)
try(dir.create(folder))
setwd(folder)

data(grav)
grav <- calc.genoprob(grav, step=1)
grav <- reduce2grid(grav)
grav = sim.geno(grav,stepwidth = 'fixed',step=1)
phecol <- 1:nphe(grav)
out <- scanone(grav, phe=phecol, method="hk")

Y = grav$pheno
X = do.call(cbind,lapply(grav$geno,function(x) x$data))
X = apply(X,2,function(x) {x[is.na(x)] = mean(x,na.rm=T);x})
Chr = do.call(c,sapply(1:length(grav$geno),function(x) rep(sprintf('Chr%d',x),ncol(grav$geno[[x]]$data))))
X = do.call(cbind,lapply(grav$geno,function(x) apply(x$draws,c(1,2),mean)))
Chr = do.call(c,sapply(1:length(grav$geno),function(x) rep(sprintf('Chr%d',x),ncol(grav$geno[[x]]$draws))))

data = data.frame(ID = 1:nrow(Y))


times <- attr(grav, "time")
Y = Y - mean(unlist(Y))
Y = Y / sd(unlist(Y))
observations = data.frame(Y = c(t(Y)),Time = times, ID = rep(data$ID,each = ncol(Y)))
ggplot(observations,aes(x=Time,y=Y,group=ID)) + geom_line()

data_model_parameters = list(
  observations = observations,
  individual_model = Y~bs(Time,df=20,intercept=F), #poly(Time,2) +
  resid_Y_prec_shape = 2,
  resid_Y_prec_rate = 1/100
)

# initialize priors
run_parameters = BSFG_control(
  sampler = 'fast_BSFG',
  # sampler = 'general_BSFG',
  scale_Y = TRUE,
  simulation = FALSE,
  h2_divisions = 2,
  h2_step_size = .3,
  burn = 00,
  k_init = 10
)

priors = list(
  # fixed_resid_var = list(V = 5e5,   nu = 2.001),
  fixed_resid_var = list(V = 1,     nu = 3),
  QTL_resid_var = list(V = 1/10000,     nu = 3*1000),
  QTL_factors_var = list(V = 1/10000,     nu = 3*1000),
  tot_Y_var = list(V = 0.5,   nu = 3),
  tot_F_var = list(V = 18/20, nu = 20),
  # tot_F_var = list(V = 1, nu = 10000),
  delta_1   = list(shape = 2.1,  rate = 1/20),
  delta_2   = list(shape = 2, rate = 1),
  Lambda_df = 3,
  B_df      = 3,
  B_F_df    = 3
)
K = diag(1,nrow(data))
rownames(K) = data$ID
BSFG_state = BSFG_init(observations$Y, model=~1+(1|ID), data,
                       priors=priors,run_parameters=run_parameters,
                       data_model = regression_spline_model, data_model_parameters = data_model_parameters,
                       posteriorSample_params = c('Lambda','U_F','F','delta','tot_F_prec','F_h2','tot_Eta_prec','resid_h2', 'B', 'B_F', 'tau_B','tau_B_F','cis_effects','U_R','prec_B','prec_B_F'),
                       posteriorMean_params = c(),QTL_factors = X#, QTL_resid = X
)
BSFG_state$current_state$F_h2
# BSFG_state$current_state$B_F[] = 0
# BSFG_state$current_state$B[] = 0

h2_divisions = run_parameters$h2_divisions
BSFG_state$priors$h2_priors_resids = with(BSFG_state$data_matrices, sapply(1:ncol(h2s_matrix),function(x) {
  h2s = h2s_matrix[,x]
  pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10)
}))
BSFG_state$priors$h2_priors_resids[] = 0
BSFG_state$priors$h2_priors_resids[1] = 1
BSFG_state$priors$h2_priors_resids = BSFG_state$priors$h2_priors_resids/sum(BSFG_state$priors$h2_priors_resids)
BSFG_state$priors$h2_priors_factors = BSFG_state$priors$h2_priors_resids

rescale_model_matrices = function(var_Eta_factor,current_state){
  current_state = within(current_state,{
    var_Eta = var_Eta * var_Eta_factor
    Lambda = sweep(Lambda,1,sqrt(var_Eta_factor),'/')
    Plam = sweep(Plam,1,sqrt(var_Eta_factor),'/')
    B = sweep(B,2,sqrt(var_Eta_factor),'/')
    U_R = sweep(U_R,2,sqrt(var_Eta_factor),'/')
    tot_Eta_prec = tot_Eta_prec * var_Eta_factor
  })
  current_state$model_matrices = NULL
  current_state
}

n_samples = 200
scanone_results = list()
grav_eta = grav
for(i  in 1:100) {
  print(sprintf('Run %d',i))
  BSFG_state = sample_BSFG(BSFG_state,n_samples,grainSize=1)
  trace_plot(BSFG_state$Posterior$tot_F_prec[,1,])
  if(BSFG_state$current_state$nrun < BSFG_state$run_parameters$burn) {
    BSFG_state = reorder_factors(BSFG_state)
  }
  p1 = posterior_plot(BSFG_state$Posterior$B_F[,-1,1],colorGroup = Chr)
  p2 = trace_plot(BSFG_state$Posterior$B_F[,-1,1])
  BSFG_state = save_posterior_chunk(BSFG_state)
  print(BSFG_state)
  plot(BSFG_state)
  plot(with(BSFG_state$current_state,model_matrices[[1]]$X %*% Eta[1,]))
  lines(unlist(Y[1,]),lwd=2,col=2)
  Eta = BSFG_state$current_state$Eta
  colnames(Eta) = 1:ncol(Eta)
  a = melt(Eta)
  print(ggplot(a,aes(x=X2,y=value,group=X1)) + geom_line())
  a = melt(sweep(Eta,2,sqrt(BSFG_state$current_state$var_Eta),'*')[,-c(1:3)])
  print(ggplot(a,aes(x=X2,y=value,group=X1)) + geom_line())
  print(p1)
  print(p2)
  plot(BSFG_state$data_matrices$Y,BSFG_state$current_state$Y_fitted);abline(0,1)
  var_Eta = apply(BSFG_state$current_state$Eta,2,var)
  print(cbind(var_Eta,apply(BSFG_state$current_state$model_matrices[[1]]$X,2,var)))
  if(i < -5) {
    BSFG_state$current_state = rescale_model_matrices(var_Eta,BSFG_state$current_state)
  } else if(i < 10 || (i-1) %% 10 == 0) {
    # BSFG_state$current_state = update_k(BSFG_state)
    BSFG_state = reorder_factors(BSFG_state)
    BSFG_state = clear_Posterior(BSFG_state)
  }
  grav_eta$pheno = BSFG_state$current_state$F
  phecol = 1:ncol(grav_eta$pheno)
  scanone_results[[i]] = scanone(grav_eta, phe=phecol, method="hk")
}

spline_X = model.matrix(BSFG_state$current_state$Terms,data = data.frame(Time = seq(min(times),max(times),length=30)))
QTL_effects = 2:ncol(BSFG_state$data_matrices$X_F)
QTL_Chr = Chr

Posterior = reload_Posterior(BSFG_state)
total_samples = Posterior$total_samples
posterior_plot(Posterior$B_F[,-1,1],colorGroup = Chr)

total_samples = Posterior$total_samples
# total_samples = 500
B = aperm(array(sapply(1:total_samples,function(i) {
  with(Posterior,{
    0*B[i,,] + B_F[i,,] %*% t(Lambda[i,,])
  })
}),dim = c(dim(Posterior$B_F)[2],dim(Posterior$Lambda)[2],total_samples)),c(3,1,2))

image(t(Matrix(apply(B,c(2,3),mean)))[dim(B)[3]:1,],at=c(-Inf,seq(-1,1,length=11),Inf))
posterior_plot(B[,-1,1],colorGroup = Chr)

