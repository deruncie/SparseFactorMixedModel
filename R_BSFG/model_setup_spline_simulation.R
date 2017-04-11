library(microbenchmark)
library(MCMCpack)
library(BSFG)

# set the directory to the location of the setup.RData or setup.mat file


# # choose a seed for the random number generator. This can be a random seed (for analysis), or you can choose your seed so that
# # you can repeat the MCMC exactly
seed = 1
new_halfSib_spline_simulation('Sim_FE_1', nSire=25,nRep=3,p=20, Time = 1:60, k=4, k_G=2, i_Va = 0.2, i_Ve = 0.2)
set.seed(seed)

# create a folder for holding the posterior samples of the current chain (multiple folders could be used for different chains)
rep = "2"
folder = sprintf('R_spline_rep_%s',rep)
try(dir.create(folder))
setwd(folder)


# initialize priors
run_parameters = BSFG_control(
  # sampler = 'fast_BSFG',
  sampler = 'general_BSFG',
  simulation = TRUE,
  h2_divisions = 20,
  h2_step_size = NULL,
  burn = 100
)

priors = list(
  # fixed_var = list(V = 5e5,   nu = 2.001),
  fixed_var = list(V = 1,     nu = 3),
  tot_Y_var = list(V = 0.5,   nu = 3),
  tot_F_var = list(V = 18/20, nu = 20),
  delta_1   = list(shape = 2.1,  rate = 1/20),
  delta_2   = list(shape = 3, rate = 1),
  Lambda_df = 3,
  B_df      = 3,
  B_F_df    = 3
)


print('Initializing')
# setup = load_simulation_data()
# setup = load_simulation_FE_data()
load('../setup.RData')
# options(error=recover)
BSFG_state = with(setup,BSFG_init(observations$Y, model=~X2+(1|animal), data, #factor_model_fixed = ~1,
                                  priors=priors,run_parameters=run_parameters,K_mats = list(animal = K),
                                  data_model = bs_model, data_model_parameters = list(observations = observations,df = 20,intercept = T,resid_Y_prec_shape = 2,resid_Y_prec_rate = 1),
                                  setup = setup))
BSFG_state$current_state$F_h2

h2_divisions = run_parameters$h2_divisions
BSFG_state$priors$h2_priors_resids = with(BSFG_state$data_matrices, sapply(1:ncol(h2s_matrix),function(x) {
  h2s = h2s_matrix[,x]
  pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10)
}))
BSFG_state$priors$h2_priors_resids = BSFG_state$priors$h2_priors_resids/sum(BSFG_state$priors$h2_priors_resids)
BSFG_state$priors$h2_priors_factors = c(h2_divisions-1,rep(1,h2_divisions-1))/(2*(h2_divisions-1))

save(BSFG_state,file="BSFG_state.RData")

BSFG_state = clear_Posterior(BSFG_state)


# # optional: To load from end of previous run, run above code, then run these lines:
# load('BSFG_state.RData')
# load('current_state.RData')
# load('Posterior/Posterior_base.RData')
# BSFG_state$current_state = current_state
# BSFG_state$Posterior = Posterior

# Run Gibbs sampler. Run in smallish chunks. Output can be used to re-start chain where it left off.
# burn in

n_samples = 100;
for(i  in 1:70) {
    print(sprintf('Run %d',i))
    # microbenchmark(
    #     BSFG_sampler(BSFG_state,2,detectCores()),
    #     BSFG_sampler(BSFG_state,2,1),
    #     times = 10
    #     )
    BSFG_state = sample_BSFG(BSFG_state,n_samples,ncores=1)
    if(BSFG_state$current_state$nrun < BSFG_state$run_parameters$burn) {
      BSFG_state = reorder_factors(BSFG_state)
    }
    BSFG_state = save_posterior_chunk(BSFG_state)
    print(BSFG_state)
    plot(BSFG_state)
}

BSFG_state$Posterior = reload_Posterior(BSFG_state)
plot(BSFG_state$Posterior$Eta,setup$Eta);abline(0,1)
setup$observations$Y_fitted = sapply(1:nrow(setup$observations),function(i) predict(BSFG_state$current_state$coefficients,newx = setup$observations$covariate[i]) %*% BSFG_state$Posterior$Eta[setup$observations$ID[i],])
with(setup$observations,plot(Y,Y_fitted));abline(0,1)
plot(c(with(BSFG_state$current_state,U_F %*% t(Lambda))),with(setup,c(as.matrix(U_F %*% t(error_factor_Lambda)))));abline(0,1)
plot(c(with(BSFG_state$current_state,U_R + U_F %*% t(Lambda))),with(setup,c(as.matrix(U_R + U_F %*% t(error_factor_Lambda)))));abline(0,1)
plot(BSFG_state$current_state$U_R,setup$U_R)

ggplot(setup$observations,aes(x=covariate,y=Y)) + geom_line(aes(group = ID,color = ID))
ggplot(setup$observations,aes(x=covariate,y=Y_fitted)) + geom_line(aes(group = ID,color = ID))

library(shinystan)
library(MCMCpack)
samples = BSFG_state$Posterior$F_h2
samples = list(A=BSFG_state$Posterior$Lambda[,,1])
samples[[1]][] = samples[[1]][,order(-abs(colMeans(samples[[1]])))]
my_sso <- as.shinystan(samples)
my_sso <- launch_shinystan(my_sso)
