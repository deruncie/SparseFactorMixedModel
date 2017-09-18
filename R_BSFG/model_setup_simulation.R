library(microbenchmark)
library(MCMCpack)
library(BSFG)
library(Matrix)


# # choose a seed for the random number generator. This can be a random seed (for analysis), or you can choose your seed so that
# # you can repeat the MCMC exactly
seed = 2
new_halfSib_simulation('Sim_FE_1', nSire=50,nRep=10,p=100, b=5, factor_h2s= c(rep(0,5),rep(0.7,5)),Va = 2, Ve = 2,Vb = 2)
set.seed(seed)
load('setup.RData')


# create a folder for holding the posterior samples of the current chain (multiple folders could be used for different chains)
rep = "2"
folder = sprintf('R_rep_%s',rep)
try(dir.create(folder))
setwd(folder)

print('Initializing')
load('../setup.RData')

Y = setup$Y
K = setup$K
data = setup$data
X = setup$X

# setup$data$Group = gl(3,1,length = nrow(setup$data))

run_parameters = BSFG_control(
  sampler = 'fast_BSFG',
  # sampler = 'fast_missing_BSFG',
  # sampler = 'general_BSFG',
  scale_Y = FALSE,
  simulation = TRUE,
  h2_divisions = 200,
  h2_step_size = .3,
  burn = 100
)

priors = BSFG_priors(
  fixed_var = list(V = 1,     nu = 3),
  # tot_Y_var = list(V = 0.5,   nu = 3),
  tot_Y_var = list(V = 0.5,   nu = 10),
  tot_F_var = list(V = 18/20, nu = 20),
  h2_priors_resids_fun = function(h2s,n) 1,#pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
  h2_priors_factors_fun = function(h2s,n) 1,#ifelse(h2s == 0,n,n/(n-1)),
  Lambda_prior = list(
    sampler = sample_Lambda_prec_ARD,
    Lambda_df = 3,
    delta_1   = list(shape = 2.1,  rate = 1/20),
    delta_2   = list(shape = 3, rate = 1)
  ),
  # Lambda_prior = list(
  #   sampler = sample_Lambda_prec_ARD_v2,
  #   Lambda_df = 3,
  #   delta_1   = list(shape = 2,  rate = 40),
  #   delta_2   = list(shape = 2, rate = 3)
  # ),
  # Lambda_prior = list(
  #   sampler = sample_Lambda_prec_TPB,
  #   Lambda_A      = .5,
  #   Lambda_B      = .5,
  #   delta_1   = list(shape = 2.1,  rate = 1/20),
  #   delta_2   = list(shape = 3, rate = 1)
  # ),
  B_prior = list(
    sampler = sample_B_prec_ARD,
    B_df      = 3,
    B_F_df    = 3
  )
  # B_prior = list(
  #   sampler = sample_B_prec_TPB,
  #   B_A      = .5,
  #   B_B      = .5,
  #   B_omega  = 1/10,
  #   B_F_A      = .5,
  #   B_F_B      = .5,
  #   B_F_omega  = 1/10
  # )
)

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
data$ID = sample(1:nrow(data))
# diag(K) = diag(K) + 1e-6
# Y_full = Y
# nGroup = 5
# for(i in 1:nGroup) {
#   Y[(i+(1:nrow(Y)-1)) %% nGroup != 0,1:(100/nGroup)+(i-1)*(100/nGroup)] = NA
# }
# i = sample(1:500,250)
# Y[i,1:50] = NA
# Y[-i,-c(1:50)] = NA
# BSFG_state = BSFG_init(Y, model=~Fixed1+Fixed2+Fixed3+Fixed4+(1|Sire), data,# factor_model_fixed = ~0,
BSFG_state = BSFG_init(Y, model=~Fixed1+Fixed2+Fixed3+Fixed4+(1|animal), data,# factor_model_fixed = ~0,
# BSFG_state = BSFG_init(Y, model=~Fixed1+Fixed2+Fixed3+Fixed4+(1|ID), data, #factor_model_fixed = ~0,
# BSFG_state = BSFG_init(Y, model=~1+(1|animal), data, factor_model_fixed = ~0,
                                  K_mats = list(animal = K),
                                  run_parameters=run_parameters,
                                  priors=priors,
                                  setup = setup)
# BSFG_state$current_state$F_h2
# BSFG_state$priors$h2_priors_resids
# BSFG_state$priors$h2_priors_factors

save(BSFG_state,file="BSFG_state.RData")

BSFG_state = clear_Posterior(BSFG_state)

# class(BSFG_state) = c('fast_BSFG_noU',class(BSFG_state))

# load('current_state.RData')
# load('BSFG_state.RData')
# BSFG_state$current_state = current_state
# load('Posterior/Posterior_base.RData')
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
class(BSFG_state)[1] = 'general_BSFG'
BSFG_state = reorder_factors(BSFG_state)
BSFG_state = clear_Posterior(BSFG_state)
n_samples = 100;
for(i  in 1:110) {
    if(i %% 6 == 0){
      BSFG_state = reorder_factors(BSFG_state)
      BSFG_state = clear_Posterior(BSFG_state)
    }
    print(sprintf('Run %d',i))
    BSFG_state = sample_BSFG(BSFG_state,n_samples,grainSize=1,ncores = 1)
    if(BSFG_state$Posterior$total_samples>0) trace_plot(BSFG_state$Posterior$tot_F_prec[,1,])
    if(BSFG_state$Posterior$total_samples>0) trace_plot(log(BSFG_state$Posterior$delta[,1,]))
    if(BSFG_state$current_state$nrun < BSFG_state$run_parameters$burn) {
      BSFG_state = reorder_factors(BSFG_state)
      # BSFG_state$current_state = update_k(BSFG_state)
      BSFG_state = clear_Posterior(BSFG_state)
      BSFG_state$run_parameters$burn = max(c(BSFG_state$run_parameters$burn,BSFG_state$current_state$nrun+100))
      print(BSFG_state$run_parameters$burn)
    }
    try(print(apply(abs(cor(as.matrix(setup$F),apply(BSFG_state$Posterior$F,c(2,3),mean))),1,max)))
    BSFG_state = save_posterior_chunk(BSFG_state)
    print(BSFG_state)
    plot(BSFG_state)
}

BSFG_state$Posterior = reload_Posterior(BSFG_state)

pmG2 = get_posterior_mean(BSFG_state,tcrossprod(sweep(Lambda,2,sqrt(F_h2),'*')) + diag(resid_h2[1,]/tot_Eta_prec[1,]))
mean((pmG - setup$G)^2)

plot(BSFG_state$Posterior$Eta[is.na(Y)],Y_full[is.na(Y)]);abline(0,1)
cor(c(BSFG_state$Posterior$Eta[is.na(Y)]),c(Y_full[is.na(Y)]))
BSFG_state1 = BSFG_state

post_mean_Lambda = get_posterior_mean(BSFG_state,Lambda)
post_mean_F = get_posterior_mean(BSFG_state,U_F)
post_mean_Lambda = get_posterior_mean(BSFG_state, F %*% t(Lambda))
post_mean_P = get_posterior_FUN(BSFG_state,tcrossprod(Lambda))

XB = get_posterior_FUN(BSFG_state,'X %*% B')
XB2 = get_posterior_FUN(BSFG_state,X %*% B)
XB3 = get_posterior_FUN(BSFG_state,{a=X %*% B;a+Eta})
XB = get_posterior_FUN(BSFG_state,B)
B = get_posterior_mean(BSFG_state,'B')
B2 = get_posteriorMean(BSFG_state,'B')

U = get_posterior_mean(BSFG_state,U_F %*% t(Lambda) + U_R)
plot(U,with(setup,U_F %*% t(error_factor_Lambda) + U_R));abline(0,1)


U = get_posterior_mean(BSFG_state,U_F %*% t(Lambda) + U_R)
plot(U,with(setup,U_F %*% t(error_factor_Lambda) + U_R))


G = get_posterior_FUN(BSFG_state,tcrossprod(sweep(Lambda,2,sqrt(F_h2),'*')) + diag(resid_h2[1,]/tot_Eta_prec[1,]))
i = 1
G_samples = get_posterior_FUN(BSFG_state,Lambda %*% diag(F_h2[i,]) %*% t(Lambda) + diag(resid_h2[i,]/tot_Eta_prec[i,]))
U_samples = get_posterior_FUN(BSFG_state, Z %*% U_F %*% t(Lambda) + Z %*% U_R)
image(cov2cor(get_posterior_mean(G)))
plot(get_posterior_mean(G)~setup$G,ylim=c(-4,4));abline(0,1)

library(shinystan)
library(MCMCpack)
samples = BSFG_state$Posterior$F_h2
samples = list(A=BSFG_state$Posterior$Lambda[,,1])
samples[[1]][] = samples[[1]][,order(-abs(colMeans(samples[[1]])))]
my_sso <- as.shinystan(samples)
my_sso <- launch_shinystan(my_sso)
