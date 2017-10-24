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
times <- attr(grav, "time")
grav <- calc.genoprob(grav, step=1)
grav <- reduce2grid(grav)
# grav <- fill.geno(grav)
# grav = sim.geno(grav,stepwidth = 'fixed',step=1)
phecol <- 1:nphe(grav)
out <- scanone(grav, phe=phecol, method="hk")

Y = grav$pheno
X = do.call(cbind,lapply(grav$geno,function(x) x$prob[,,1]))
Chr = do.call(c,sapply(1:length(grav$geno),function(x) rep(sprintf('Chr%d',x),ncol(grav$geno[[x]]$data))))
# X = do.call(cbind,lapply(grav$geno,function(x) apply(x$draws,c(1,2),mean)))
# Chr = do.call(c,sapply(1:length(grav$geno),function(x) rep(sprintf('Chr%d',x),ncol(grav$geno[[x]]$draws))))

data = data.frame(ID = 1:nrow(Y))


Y = Y - mean(unlist(Y))
Y = Y / sd(unlist(Y))
observations = data.frame(Y = c(t(Y)),Time = times, ID = rep(data$ID,each = ncol(Y)))
ggplot(observations,aes(x=Time,y=Y,group=ID)) + geom_line()

observation_setup = list(
  observation_model = regression_model,
  observations = observations,
  individual_model = Y~bs(Time,df=20,intercept=F),
  resid_Y_prec_shape = 2,
  resid_Y_prec_rate = 1/100,
  do_not_penalize_bs = FALSE
)

# initialize priors
run_parameters = BSFG_control(
  scale_Y = FALSE,
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
  h2_priors_resids_fun = function(h2s,n) 1,#pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
  h2_priors_factors_fun = function(h2s,n) ifelse(h2s == 0 | h2s >= .99,n,n/(n-1)),
  Lambda_prior = list(
    sampler = sample_Lambda_prec_ARD,
    Lambda_df = 3,
    delta_1   = list(shape = 2.1,  rate = 1/20),
    delta_2   = list(shape = 6, rate = 2)
  ),
  B_prior = list(
    sampler = sample_B_prec_ARD,
    B_df      = 3,
    B_F_df    = 1
  )
)

K = diag(1,nrow(data))
rownames(K) = data$ID
BSFG_state = BSFG_init(observation_setup, model=~1+(1|ID), data,
                       priors=priors,run_parameters=run_parameters,
                       QTL_factors = X,
                       QTL_resid = X
                       # posteriorSample_params = c('Lambda','U_F','F','delta','tot_F_prec','F_h2','tot_Eta_prec','resid_h2', 'B', 'B_F', 'tau_B','tau_B_F','cis_effects','U_R','prec_B','prec_B_F'),
                       # posteriorMean_params = c(),QTL_factors = X#, QTL_resid = X
)
mm_rot = BSFG_state$run_parameters$observation_model_parameters$observation_setup$mm_rotation
Image(mm_rot)
BSFG_state$run_variables$resid_intercept = FALSE
# BSFG_state$current_state$F_h2
# BSFG_state$current_state$B_F[] = 0
# BSFG_state$current_state$B[] = 0
#
# h2_divisions = run_parameters$h2_divisions
# BSFG_state$priors$h2_priors_resids = with(BSFG_state$data_matrices, sapply(1:ncol(h2s_matrix),function(x) {
#   h2s = h2s_matrix[,x]
#   pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10)
# }))
# BSFG_state$priors$h2_priors_resids[] = 0
# BSFG_state$priors$h2_priors_resids[1] = 1
# BSFG_state$priors$h2_priors_resids = BSFG_state$priors$h2_priors_resids/sum(BSFG_state$priors$h2_priors_resids)
# BSFG_state$priors$h2_priors_factors = BSFG_state$priors$h2_priors_resids

rescale_Eta = function(BSFG_state){
  var_Eta_factor = apply(BSFG_state$current_state$Eta,2,var)
  current_state = within(BSFG_state$current_state,{
    var_Eta = var_Eta * var_Eta_factor
    Lambda = sweep(Lambda,1,sqrt(var_Eta_factor),'/')
    Plam = sweep(Plam,1,sqrt(var_Eta_factor),'/')
    B = sweep(B,2,sqrt(var_Eta_factor),'/')
    U_R = sweep(U_R,2,sqrt(var_Eta_factor),'/')
    tot_Eta_prec = tot_Eta_prec * var_Eta_factor
  })
  current_state
}

n_samples = 200
scanone_results = list()
scanone_results2 = list()
grav_eta = grav
n_samples = 100;
BSFG_state = clear_Posterior(BSFG_state)
for(i  in 1:110) {
  var_Eta_factor = apply(BSFG_state$current_state$Eta,2,var)
  print(var_Eta_factor)
  print(apply(BSFG_state$current_state$F,2,var))
  if(i %%5 == 0 && i < 20) {
    # BSFG_state$current_state$var_Eta = var_Eta_factor
    rescale_Eta = function(BSFG_state, var_Eta_factor){
      current_state = within(BSFG_state$current_state,{
        var_Eta = var_Eta * var_Eta_factor
        Eta = sweep(Eta,2,sqrt(var_Eta_factor),'/')
        Lambda = sweep(Lambda,1,sqrt(var_Eta_factor),'/')
        Plam = sweep(Plam,1,sqrt(var_Eta_factor),'*')
        B = sweep(B,2,sqrt(var_Eta_factor),'/')
        U_R = sweep(U_R,2,sqrt(var_Eta_factor),'/')
        tot_Eta_prec = tot_Eta_prec * var_Eta_factor
        delta = delta * sqrt(mean(var_Eta_factor))
      })
      current_state
    }
    var_Eta_factor = c(var_Eta_factor[1],rep(median(var_Eta_factor[-1]),length(var_Eta_factor)-1))
    BSFG_state$current_state = rescale_Eta(BSFG_state,var_Eta_factor)
  }
  if(i < 10) {
    # BSFG_state$data_matrices$X_F = 0*BSFG_state$data_matrices$X_F
    # BSFG_state$data_matrices$QTL_factors_X = 0*BSFG_state$data_matrices$QTL_factors_X
    # BSFG_state$priors$fixed_factors_prec_rate[-1] = 1e-10
    # BSFG_state$current_state$B_F_prec[] = 1e10 + 0*BSFG_state$current_state$B_F_prec
    # BSFG_state$current_state$B_F[] = 0*BSFG_state$current_state$B_F
    # BSFG_state$current_state = update_k(BSFG_state)

    # if(i>1) BSFG_state = rescale_factors_F(BSFG_state)
    BSFG_state = reorder_factors(BSFG_state)
    BSFG_state = clear_Posterior(BSFG_state)
  }
  print(sprintf('Run %d',i))
  BSFG_state = sample_BSFG(BSFG_state,n_samples,ncores=1)
  if(BSFG_state$Posterior$total_samples>0) trace_plot(BSFG_state$Posterior$tot_F_prec[,1,])
  if(BSFG_state$current_state$nrun < BSFG_state$run_parameters$burn) {
    BSFG_state = reorder_factors(BSFG_state)
  }
  BSFG_state = save_posterior_chunk(BSFG_state)
  print(BSFG_state)
  plot(BSFG_state)
  par(mfrow=c(1,2))
  plot(BSFG_state$run_parameters$observation_model_parameters$observation_setup$Y,BSFG_state$current_state$Y_fitted);abline(0,1)
  boxplot(BSFG_state$current_state$Eta)
  # dev.off()
  # BSFG_state$priors$fixed_factors_prec_rate[-1] = 3
  # BSFG_state$data_matrices$X_F = X_F
  # BSFG_state$data_matrices$QTL_factors_X = QTL_factors_X
  # if(i %%5 == 0) saveRDS(BSFG_state$current_state,file = sprintf('current_state_%d.rds',i))

  grav_eta$pheno = BSFG_state$current_state$F
  phecol = 1:ncol(grav_eta$pheno)
  scanone_results[[i]] = scanone(grav_eta, phe=phecol, method="hk")
  print(apply(scanone_results[[i]][,-c(1:2)],2,max))
  grav_eta$pheno = with(BSFG_state,data_matrices$X_F %*% current_state$B_F)
  phecol = 1:ncol(grav_eta$pheno)
  scanone_results2[[i]] = scanone(grav_eta, phe=phecol, method="hk")
  print(apply(scanone_results2[[i]][,-c(1:2)],2,max))
}

mm_rot = BSFG_state$run_parameters$observation_model_parameters$observation_setup$mm_rotation
vE = BSFG_state$current_state$var_Eta

boxplot(sweep(BSFG_state$current_state$Eta,2,sqrt(vE),'*') %*% t(mm_rot))

a = sweep(BSFG_state$current_state$Eta,2,sqrt(vE),'*') %*% t(mm_rot)
plot(NA,NA,xlim=c(0,ncol(a)),ylim = range(a))
b=apply(a[sample(1:nrow(a),10),],1,lines)

Image(sweep(BSFG_state$current_state$Lambda,1,sqrt(vE),'*'))
Image(mm_rot %*% sweep(BSFG_state$current_state$Lambda,1,sqrt(vE),'*'))
Image(tcrossprod(sweep(BSFG_state$current_state$Lambda,1,sqrt(vE),'*')))
Image(tcrossprod(mm_rot %*% sweep(BSFG_state$current_state$Lambda,1,sqrt(vE),'*')))


Terms = BSFG_state$run_parameters$observation_model_parameters$observation_setup$Terms
spline_X = model.matrix(Terms,data = data.frame(Time = seq(min(times),max(times),length=241))) %*% mm_rot
spline_X = sweep(spline_X,2,sqrt(vE),'*')
QTL_effects = 2:ncol(BSFG_state$data_matrices$X_F)
QTL_Chr = Chr

Posterior = reload_Posterior(BSFG_state)

pm_F = apply(Posterior$F,c(2,3),mean)
grav_eta$pheno = pm_F
phecol = 1:ncol(grav_eta$pheno)
out_F = scanone(grav_eta, phe=phecol, method="hk")
plot(out_F,lodcolumn = 1:3)
plot(out_F,lodcolumn = 1:3+3)
print(apply(out_F[,-c(1:2)],2,max))

pm_XBf = BSFG_state$data_matrices$X_F %*% apply(Posterior$B_F,c(2,3),mean)
grav_eta$pheno = pm_XBf
phecol = 1:ncol(grav_eta$pheno)
out_XBf = scanone(grav_eta, phe=phecol, method="hk")
print(apply(out_XBf[,-c(1:2)],2,max))

total_samples = Posterior$total_samples
posterior_plot(Posterior$B[,-1,21])

total_samples = Posterior$total_samples
# total_samples = 500
B = aperm(array(sapply(1:total_samples,function(i) {
  with(Posterior,{
    B[i,,] + B_F[i,,] %*% t(Lambda[i,,])
  })
}),dim = c(dim(Posterior$B_F)[2],dim(Posterior$Lambda)[2],total_samples)),c(3,1,2))

B_mean = apply(B,c(2,3),mean) %*% t(spline_X)
Image(t(B_mean)[dim(B_mean)[2]:1,-1],aspect=1/1.5)
image(t(Matrix(apply(B,c(2,3),mean)))[dim(B)[3]:1,],at=c(-Inf,seq(-.01,.01,length=11),Inf),aspect=1/1.5)
posterior_plot(B[,-1,1],colorGroup = Chr)

grav_eta$pheno = apply(Posterior$F,c(2,3),mean)
phecol = 1:ncol(grav_eta$pheno)
scanone_pm = scanone(grav_eta, phe=phecol, method="hk")
apply(scanone_pm[,-c(1:2)],2,max)
