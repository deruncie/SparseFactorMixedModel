library(microbenchmark)
library(MCMCpack)
library(BSFG)

# set the directory to the location of the setup.RData or setup.mat file


# # choose a seed for the random number generator. This can be a random seed (for analysis), or you can choose your seed so that
# # you can repeat the MCMC exactly
seed = 1
new_halfSib_simulation('Sim_FE_1', nSire=50,nRep=10,p=100, b=5, k=10, k_G=5, i_Va = 0.2, i_Ve = 0.2)
set.seed(seed)

# create a folder for holding the posterior samples of the current chain (multiple folders could be used for different chains)
rep = "2"
folder = sprintf('R_rep_%s',rep)
try(dir.create(folder))
setwd(folder)


# initialize priors
run_parameters = BSFG_control(
  sampler = 'fast_BSFG',
  # sampler = 'general_BSFG',
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
load('../setup.RData')
setup$data$Group = gl(3,1,length = nrow(setup$data))

# to test non-PD K matrix:
# X = matrix(sample(c(0,1),dim(setup$K)^2,replace=T),dim(setup$K)[1])
# Xs = sweep(X,2,colMeans(X),'-')
# K_new = tcrossprod(Xs)
# rownames(K_new) = rownames(setup$K)
# BSFG_state = with(setup,BSFG_init(Y, model=~Fixed1+Fixed2+Fixed3+Fixed4+(1+Fixed2 + Fixed3|Sire)+(Group|animal), #
                                  # data,priors=priors,run_parameters = run_parameters,K_mats = list(animal = K_new),
                                  # setup = setup))
# setup$Y[1:3] = NA
# setup$Y[sample(1:prod(dim(setup$Y)),5000)] = NA
BSFG_state = with(setup,BSFG_init(Y, model=~Fixed1+Fixed2+Fixed3+Fixed4+(1|animal), data, #factor_model_fixed = ~1,
                                  priors=priors,run_parameters=run_parameters,K_mats = list(animal = K),
                                  setup = setup))
BSFG_state$current_state$F_h2

h2_divisions = run_parameters$h2_divisions
BSFG_state$priors$h2_priors_resids = with(BSFG_state$data_matrices, sapply(1:ncol(h2s_matrix),function(x) {
    h2s = h2s_matrix[,x]
    pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10)
}))
BSFG_state$priors$h2_priors_resids = BSFG_state$priors$h2_priors_resids/sum(BSFG_state$priors$h2_priors_resids)
BSFG_state$priors$h2_priors_factors = c(h2_divisions-1,rep(1,h2_divisions-1))/(2*(h2_divisions-1))

# BSFG_state$priors$h2_priors_factors = with(BSFG_state$data_matrices, sapply(1:nrow(h2s_matrix),function(x) {
#     h2s = h2s_matrix[x,]
#     pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10)
# }))

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

n_samples = 100;
for(i  in 1:70) {
    print(sprintf('Run %d',i))
    BSFG_state = sample_BSFG(BSFG_state,n_samples,grainSize=1)
    if(BSFG_state$current_state$nrun < BSFG_state$run_parameters$burn) {
      BSFG_state = reorder_factors(BSFG_state)
    }
    BSFG_state = save_posterior_chunk(BSFG_state)
    print(BSFG_state)
    plot(BSFG_state)
}

library(shinystan)
library(MCMCpack)
samples = BSFG_state$Posterior$F_h2
samples = list(A=BSFG_state$Posterior$Lambda[,,1])
samples[[1]][] = samples[[1]][,order(-abs(colMeans(samples[[1]])))]
my_sso <- as.shinystan(samples)
my_sso <- launch_shinystan(my_sso)
