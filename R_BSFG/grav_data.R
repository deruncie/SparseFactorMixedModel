library(microbenchmark)
library(MCMCpack)
library(BSFG)
library(ggplot2)
library(Rcpp)
library(qtl)
library(qtlcharts)
library(reshape)

# create a folder for holding the posterior samples of the current chain (multiple folders could be used for different chains)
rep = "4"
folder = sprintf('grav_rep_%s',rep)
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


# initialize priors
run_parameters = BSFG_control(
  sampler = 'fast_BSFG',
  # sampler = 'general_BSFG',
  scale_Y = TRUE,
  simulation = FALSE,
  h2_divisions = 10,
  h2_step_size = .3,
  burn = 00,
  k_init = 8
)

priors = list(
  # fixed_var = list(V = 5e5,   nu = 2.001),
  fixed_var = list(V = 1/10000,     nu = 3*1000),
  tot_Y_var = list(V = 0.5,   nu = 3),
  tot_F_var = list(V = 18/20, nu = 20),
  delta_1   = list(shape = 2.1,  rate = 1/20),
  delta_2   = list(shape = 3, rate = 1),
  Lambda_df = 10,
  B_df      = 3,
  B_F_df    = 3
)
K = diag(1,nrow(data))
rownames(K) = data$ID
BSFG_state = BSFG_init(Y, model=~1+(1|ID), data,
                       priors=priors,run_parameters=run_parameters,
                       posteriorSample_params = c('Lambda','U_F','F','delta','tot_F_prec','F_h2','tot_Eta_prec','resid_h2', 'B', 'B_F', 'tau_B','tau_B_F','cis_effects','U_R','prec_B','prec_B_F'),
                       posteriorMean_params = c()#,X_factor = X#, X_resid = X
)
BSFG_state$current_state$F_h2
# BSFG_state$current_state$B_F[] = 0

h2_divisions = run_parameters$h2_divisions
BSFG_state$priors$h2_priors_resids = with(BSFG_state$data_matrices, sapply(1:ncol(h2s_matrix),function(x) {
  h2s = h2s_matrix[,x]
  pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10)
}))
BSFG_state$priors$h2_priors_resids = BSFG_state$priors$h2_priors_resids/sum(BSFG_state$priors$h2_priors_resids)
BSFG_state$priors$h2_priors_factors = BSFG_state$priors$h2_priors_resids



n_samples = 200
for(i  in 1:100) {
  print(sprintf('Run %d',i))
  BSFG_state = sample_BSFG(BSFG_state,n_samples,grainSize=1)
  trace_plot(BSFG_state$Posterior$tot_F_prec[,1,])
  if(BSFG_state$current_state$nrun < BSFG_state$run_parameters$burn) {
    BSFG_state = reorder_factors(BSFG_state)
  }
  BSFG_state = save_posterior_chunk(BSFG_state)
  print(BSFG_state)
  plot(BSFG_state)
  if(i < 10 || (i-1) %% 20 == 0) {
    BSFG_state$current_state = update_k(BSFG_state)
    BSFG_state = reorder_factors(BSFG_state)
    BSFG_state = clear_Posterior(BSFG_state)
  }
}

Posterior = reload_Posterior(BSFG_state)
posterior_plot(Posterior$B_F[,-1,1],colorGroup = Chr)

total_samples = Posterior$total_samples
total_samples = 500
B = aperm(array(sapply(1:total_samples,function(i) {
  with(Posterior,{
    0*B[i,,] + B_F[i,,] %*% t(Lambda[i,,])
  })
}),dim = c(dim(Posterior$B_F)[2],dim(Posterior$Lambda)[2],total_samples)),c(3,1,2))

image(t(Matrix(apply(B,c(2,3),mean)))[dim(B)[3]:1,],at=c(-Inf,seq(-1,1,length=11),Inf))
posterior_plot(B[,-1,1],colorGroup = Chr)

