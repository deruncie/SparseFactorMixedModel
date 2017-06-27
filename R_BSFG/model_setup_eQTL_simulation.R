library(microbenchmark)
library(MCMCpack)
library(BSFG)

# set the directory to the location of the setup.RData or setup.mat file

# setwd("~/Box Sync/DER_projects/BSFG/R_BSFG/Sim_FE3_1")
seed = 1
set.seed(seed)

col_sha_geno_max_cor <- as.matrix(read.csv('~/Box Sync/DER_projects/NAM_CAM/BSFG_analysis/col_sha_markers_bsfg.csv', header = TRUE))
col_sha_geno_max_cor = col_sha_geno_max_cor[,1:(ncol(col_sha_geno_max_cor)/2)]
col_sha_geno_max_cor = col_sha_geno_max_cor[,sample(1:ncol(col_sha_geno_max_cor))]
new_halfSib_simulation_eQTL('Sim_eQTL_1', nSire=50,nRep=10,p=100, b=1, factor_h2s= c(rep(0,5),rep(0.3,5)),Va = 2, Ve = 2,Vb = 2,V_cis = 1,nSNP = 200,bSNP = 1,col_sha_geno_max_cor)


# create a folder for holding the posterior samples of the current chain (multiple folders could be used for different chains)
rep = "1"
folder = sprintf('R_eQTL_rep_%s',rep)
try(dir.create(folder))
setwd(folder)


# initialize priors
run_parameters = BSFG_control(
  sampler = 'fast_BSFG',
  # sampler = 'general_BSFG',
  scale_Y = FALSE,
  simulation = TRUE,
  h2_divisions = 20,
  h2_step_size = NULL,
  burn = 00,
  k_init = 10
)

priors = list(
  # fixed_var = list(V = 5e5,   nu = 2.001),
  fixed_var = list(V = 1,     nu = 3),
  QTL_resid_var = list(V = 1/1000,     nu = 3),
  QTL_factors_var = list(V = 1/1000,     nu = 3),
  tot_Y_var = list(V = 0.5,   nu = 3),
  tot_F_var = list(V = 18/20, nu = 20),
  delta_1   = list(shape = 2.1,  rate = 1/20),
  delta_2   = list(shape = 3, rate = 1),
  Lambda_df = 3,
  B_df      = 3,
  B_F_df    = 3,
  # h2_priors_resids_fun = function(h2s,n) pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
  # h2_priors_factors_fun = function(h2s,n) ifelse(h2s == 0,n,n/(n-1))
  h2_priors_resids = 1,
  h2_priors_factors = 1
)

print('Initializing')
# setup = load_simulation_data()
# setup = load_simulation_FE_data()
load('../setup_Sim_eQTL_1.RData')
# sim_data$data = with(sim_data,data.frame(TRT = X[,'TRT'],Line = colnames(Z1)[apply(Z1,1,function(x) which(x==1))]))
# sim_data$B = rbind(sim_data$mu,sim_data$B)
BSFG_state = with(setup,BSFG_init(Y, model=~1+(1|animal), data, #factor_model_fixed = ~1,
                                  QTL_factors = X_SNP,
                                  priors=priors,run_parameters=run_parameters,
                                  cis_genotypes = lapply(1:ncol(X_cis),function(x) matrix(X_cis[,x],ncol=1)),
                                  K_mats = list(animal = K),
                                  setup = setup))

# BSFG_state = with(setup,BSFG_init(list(observation_model = cis_eQTL_model,
#                                        Y = Y,
#                                        cis_genotypes = lapply(1:ncol(X_cis),function(x) matrix(X_cis[,x],ncol=1)),
#                                        ncores = parallel::detectCores()),
#                                   model=~1+(1|animal), data, #factor_model_fixed = ~1,
#                                   QTL_factors = X_SNP,
#                                   priors=priors,run_parameters=run_parameters,
#                                   K_mats = list(animal = K),
#                                   setup = setup))

BSFG_state$current_state$F_h2
BSFG_state$priors$h2_priors_resids
BSFG_state$priors$h2_priors_factors

save(BSFG_state,file="BSFG_state.RData")

BSFG_state = clear_Posterior(BSFG_state)


# # optional: To load from end of previous run, run above code, then run these lines:
# load('Posterior.RData')
# load('BSFG_state.RData')
# load('Priors.RData')
#  BSFG_state$current_state = current_state
#  BSFG_state$Posterior = Posterior
#  BSFG_state$priors = priors
# start_i = run_parameters$nrun;

# Run Gibbs sampler. Run in smallish chunks. Output can be used to re-start chain where it left off.
# burn in

n_samples = 20
for(i  in 1:100) {
  if(i < 10 || (i-1) %% 20 == 0) {
    BSFG_state$current_state = update_k(BSFG_state)
    BSFG_state = reorder_factors(BSFG_state)
    BSFG_state = clear_Posterior(BSFG_state)
  }
  print(sprintf('Run %d',i))
  BSFG_state = sample_BSFG(BSFG_state,n_samples,grainSize=1)
  trace_plot(BSFG_state$Posterior$tot_F_prec[,1,])
  # trace_plot(BSFG_state$Posterior$cis_effects[,1,1:10])
  # plot(setup$b_cis,BSFG_state$current_state$cis_effects);abline(0,1)
    if(BSFG_state$current_state$nrun < BSFG_state$run_parameters$burn) {
      BSFG_state = reorder_factors(BSFG_state)
    }
    BSFG_state = save_posterior_chunk(BSFG_state)
    print(BSFG_state)
    plot(BSFG_state)
}
BSFG_state$Posterior = reload_Posterior(BSFG_state)

trace_plot(BSFG_state$Posterior$tot_F_prec[,1,])
trace_plot(BSFG_state$Posterior$delta[,1,2:4])
trace_plot(BSFG_state$Posterior$B_F[,,1])
trace_plot(BSFG_state$Posterior$B_F[,1:3,2])
trace_plot(BSFG_state$Posterior$B_F[,,3])
trace_plot_h2s(BSFG_state$Posterior$F_h2)
trace_plot_Lambda(BSFG_state$Posterior$Lambda[1:100,,],n_factors = 10)
trace_plot(BSFG_state$Posterior$tau_B_F[1:100,1,1:10])

with(BSFG_state$Posterior,sapply(1:total_samples,function(i) colSums(1/prec_B_F[i,,])))
trace_plot(t(with(BSFG_state$Posterior,sapply(1:total_samples,function(i) colMeans(F[i,,]^2))))[,1:6])

plot(with(BSFG_state$Posterior,sapply(1:total_samples,function(i) colMeans(F[i,,]^2)))[1,],BSFG_state$Posterior$tot_F_prec[,1,1])

plot(BSFG_state$Posterior$tot_F_prec[,1,1],BSFG_state$Posterior$delta[,1,1])
plot(BSFG_state$Posterior$tot_F_prec[,1,3],BSFG_state$Posterior$B_F[,1,3])

posterior_plot(BSFG_state$Posterior$B_F[,,1])
posterior_plot(BSFG_state$Posterior$B_F[,,2])

order(BSFG_state$current_state$Lambda[,1])





posterior_plot = function(X,xlab='',ylab='',colorSig = T,ylim = NULL) {
  X = mcmc(X)
  Xdata = data.frame(ID = 1:ncol(X),Mean = colMeans(X),Median = apply(X,2,median))
  Xi = HPDinterval(X,prob = 0.95)
  Xdata$low_95 = Xi[,1]
  Xdata$high_95 = Xi[,2]
  Xi = HPDinterval(X,prob = 0.8)
  Xdata$low_80 = Xi[,1]
  Xdata$high_80 = Xi[,2]
  Xdata$color = sign(Xdata$low_95) == sign(Xdata$high_95)
  if(!colorSig) Xdata$color = 1
  if(is.null(ylim)) ylim = range(Xdata[,-1])

  p = ggplot(Xdata,aes(x=ID)) + geom_hline(yintercept = 0) +
    xlab(xlab) + ylab(ylab) + ylim(ylim)+
    geom_segment(aes(xend = ID,y = low_95,yend=high_95,color=color),size=.5) +
    geom_segment(aes(xend = ID,y = low_80,yend=high_80,color = color),size = .9) +
    geom_point(aes(y=Median,color = color)) +
    theme(legend.position = 'none')
  p
}


library(shinystan)
library(MCMCpack)
samples = BSFG_state$Posterior$F_h2
samples = list(A=BSFG_state$Posterior$Lambda[,,1])
samples[[1]][] = samples[[1]][,order(-abs(colMeans(samples[[1]])))]
my_sso <- as.shinystan(samples)
my_sso <- launch_shinystan(my_sso)
