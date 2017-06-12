library(microbenchmark)
library(MCMCpack)
library(BSFG)

# choose a seed for the random number generator. This can be a random seed (for analysis), or you can choose your seed so that
# you can repeat the MCMC exactly
seed = 1

# create a folder for holding the posterior samples of the current chain (multiple folders could be used for different chains)
rep = "1"
folder = sprintf('Rep_%s',rep)
try(dir.create(folder))
setwd(folder)

# initialize priors
run_parameters = BSFG_control(
  # sampler = 'fast_BSFG',
  sampler = 'general_BSFG',
  simulation   = FALSE,
  scale_Y = TRUE,
  h2_divisions = 20,
  h2_step_size = NULL,
  burn = 100
)

priors = list(
  fixed_var = list(V = 5e5,   nu = 2.001),
  # fixed_var = list(V = 1,     nu = 3),
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

BSFG_state = BSFG_init(Y, model=~1 + (1|Line),data,factor_model_fixed = NULL,priors=priors,run_parameters=run_parameters,K_mats = list(Line = K))

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
BSFG_state$run_parameters$simulation = FALSE
n_samples = 100;
for(i  in 1:70) {
  print(sprintf('Run %d',i))
  BSFG_state = sample_BSFG(BSFG_state,n_samples,ncores=1)
  if(BSFG_state$current_state$nrun < BSFG_state$run_parameters$burn) {
    BSFG_state = reorder_factors(BSFG_state)
  }
  BSFG_state = save_posterior_chunk(BSFG_state)
  print(BSFG_state)
  plot(BSFG_state)
}

# library(shinystan)
# library(MCMCpack)
# samples = BSFG_state$Posterior$F_h2
# samples = list(A=BSFG_state$Posterior$Lambda[,,1])
# samples[[1]][] = samples[[1]][,order(-abs(colMeans(samples[[1]])))]
# my_sso <- as.shinystan(samples)
# my_sso <- launch_shinystan(my_sso)

# F_h2 = get_posteriorMean(BSFG_state,'F_h2')
# Eta = get_posteriorMean(BSFG_state,'Eta')
# Lambda = get_posteriorMean(BSFG_state,'Lambda')
# F = get_posteriorMean(BSFG_state,'F')
# F_a = get_posteriorMean(BSFG_state,'F_a')
# B = get_posteriorMean(BSFG_state,'B')
# tot_Eta_prec = get_posteriorMean(BSFG_state,'tot_Eta_prec')
#
# LLt2 = matrix(rowMeans(apply(BSFG_state$Posterior$Lambda[sample(1:dim(BSFG_state$Posterior$Lambda)[1],50),,],1,tcrossprod)),dim(BSFG_state$Posterior$Lambda)[2])
