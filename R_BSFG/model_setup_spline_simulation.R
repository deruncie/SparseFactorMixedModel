library(microbenchmark)
library(MCMCpack)
library(BSFG)

# set the directory to the location of the setup.RData or setup.mat file


# # choose a seed for the random number generator. This can be a random seed (for analysis), or you can choose your seed so that
# # you can repeat the MCMC exactly
seed = 1
set.seed(seed)
new_halfSib_spline_simulation('Sim_FE_1', nSire=2,nRep=20,p=4, Time = 1:60, k=4, k_G=4, factor_h2s = rep(0.9,4),resid_h2 = c(.3,.3,0,0))
load('setup.RData')

ggplot(setup$observations,aes(x=covariate,y=Y)) + geom_line(aes(group = ID))

# create a folder for holding the posterior samples of the current chain (multiple folders could be used for different chains)
rep = "2"
folder = sprintf('R_spline_rep_%s',rep)
try(dir.create(folder))
setwd(folder)


# initialize priors
run_parameters = BSFG_control(
  # sampler = 'fast_BSFG',
  sampler = 'general_BSFG',
  scale_Y = FALSE,
  simulation = FALSE,
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


data_model_parameters = list(
  observations = setup$observations,
  df = 40,
  # degree=6,
  intercept = TRUE,
  individual_model = ~SPLINE,
  resid_Y_prec_shape = 2,
  resid_Y_prec_rate = 1
)

# options(error=recover)
BSFG_state = with(setup,BSFG_init(observations$Y, model=~X2+(1|animal), data, #factor_model_fixed = ~1,
                                  priors=priors,run_parameters=run_parameters,K_mats = list(animal = K),
                                  data_model = bs_model, data_model_parameters = data_model_parameters,
                                  posteriorSample_params = c('Lambda','U_F','F','delta','tot_F_prec','F_h2','tot_Eta_prec','resid_h2', 'B', 'B_F', 'tau_B','tau_B_F','cis_effects','U_R'),
                                  posteriorMean_params = c()
                                  ))
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
    BSFG_state = sample_BSFG(BSFG_state,n_samples,ncores=1)
    if(BSFG_state$current_state$nrun < BSFG_state$run_parameters$burn) {
      BSFG_state = reorder_factors(BSFG_state)
    }
    BSFG_state = save_posterior_chunk(BSFG_state)
    print(BSFG_state)
    plot(BSFG_state)
    plot(BSFG_state$data_matrices$Y,BSFG_state$current_state$Y_fitted);abline(0,1)
}

Posterior = reload_Posterior(BSFG_state)
data = setup$observations
data$Y_fitted = apply(Posterior$Y_fitted,c(2,3),mean)
data$Y_fitted_low = apply(Posterior$Y_fitted,c(2,3),function(x) HPDinterval(mcmc(x))[1])
data$Y_fitted_high = apply(Posterior$Y_fitted,c(2,3),function(x) HPDinterval(mcmc(x))[2])
# data$Y_fitted = BSFG_state$current_state$Y_fitted
ggplot(data,aes(x=covariate,y=Y)) + geom_line(aes(group = ID)) + geom_line(aes(y=Y_fitted,group=ID),color ='red')
ggplot(subset(data,ID %in% c(1,2,21,22)),aes(x=covariate,y=Y)) + geom_ribbon(aes(ymin = Y_fitted_low,ymax = Y_fitted_high,group = ID),alpha = 0.2) + geom_line(aes(group = ID)) + geom_line(aes(y=Y_fitted,group=ID),color ='red')





BSFG_state$Posterior = reload_Posterior(BSFG_state)
plot(BSFG_state$Posterior$Eta,setup$Eta);abline(0,1)
setup$observations$Y_fitted = sapply(1:nrow(setup$observations),function(i) predict(BSFG_state$current_state$coefficients,newx = setup$observations$covariate[i]) %*% BSFG_state$Posterior$Eta[setup$observations$ID[i],])
with(setup$observations,plot(Y,Y_fitted));abline(0,1)
plot(c(with(BSFG_state$current_state,U_F %*% t(Lambda))),with(setup,c(as.matrix(U_F %*% t(error_factor_Lambda)))));abline(0,1)
plot(c(with(BSFG_state$current_state,U_R + U_F %*% t(Lambda))),with(setup,c(as.matrix(U_R + U_F %*% t(error_factor_Lambda)))));abline(0,1)
plot(BSFG_state$current_state$U_R,setup$U_R)

ggplot(setup$observations,aes(x=covariate,y=Y)) + geom_line(aes(group = ID,color = ID))
ggplot(setup$observations,aes(x=covariate,y=Y_fitted)) + geom_line(aes(group = ID,color = ID))


Posterior_P = aperm(array(sapply(1:Posterior$total_samples,function(i) {
  with(Posterior,{
    Lambda[i,,] %*% t(Lambda[i,,]) + diag(1/tot_Eta_prec[i,,])
  })
}),dim = c(rep(dim(Posterior$Lambda)[2],2),Posterior$total_samples)),c(3,1,2))
dimnames(Posterior_P)[2:3] = dimnames(Posterior$Eta)[2]

library(heatmap3)
i = 1:dim(Posterior_P)[2]
heatmap3(cov2cor(apply(Posterior_P,c(2,3),mean))[i,i]^2,Rowv = NA,Colv=NA)
trace_plot(Posterior_P[,i,1])

Posterior_G = aperm(array(sapply(1:Posterior$total_samples,function(i) {
  with(Posterior,{
    Lambda[i,,] %*% diag(F_h2[i,,]) %*% t(Lambda[i,,]) + diag(resid_h2[i,,]/tot_Eta_prec[i,,])
  })
}),dim = c(rep(dim(Posterior$Lambda)[2],2),Posterior$total_samples)),c(3,1,2))
dimnames(Posterior_G)[2:3] = dimnames(Posterior$Eta)[2]
heatmap3(cov2cor(apply(Posterior_G,c(2,3),mean))[i,i]^2,Rowv = NA,Colv=NA)
trace_plot(Posterior_G[,i,6])
boxplot(Posterior_G[,i[2:10],6]);abline(h=0)


E_B = with(c(BSFG_state$data_matrices,Posterior),sapply(1:total_samples,function(i) X[1,] %*% B[i,,]))
E_B_F = with(c(BSFG_state$data_matrices,Posterior),sapply(1:total_samples,function(i) X_F[1,] %*% B_F[i,,] %*% t(Lambda[i,,])))
E_U_F = with(c(BSFG_state$data_matrices,Posterior),sapply(1:total_samples,function(i) Z[1,] %*% U_F[i,,] %*% t(Lambda[i,,])))
E_E_F = with(c(BSFG_state$data_matrices,Posterior),sapply(1:total_samples,function(i) (F[i,1,] - Z[1,] %*% U_F[i,,]-X_F[1,] %*% B_F[i,,]) %*% t(Lambda[i,,])))
E_F = with(c(BSFG_state$data_matrices,Posterior),sapply(1:total_samples,function(i) (F[i,1,]) %*% t(Lambda[i,,])))
E_U_R = with(c(BSFG_state$data_matrices,Posterior),sapply(1:total_samples,function(i) Z[1,] %*% U_R[i,,]))

i=2;plot(E_B[,i] + E_B_F[,i] + E_U_F[,i] + E_E_F[,i] + E_U_R[,i],Posterior$Eta[i,1,]);abline(0,1)

summary(apply(E_B,2,var))
summary(apply(E_F,2,var))
summary(apply(E_U_R,2,var))

Eta_mean = apply(Posterior$Eta,c(2,3),mean)
Eta_mean_R = sweep(Eta_mean,2,colMeans(Eta_mean),'-')
