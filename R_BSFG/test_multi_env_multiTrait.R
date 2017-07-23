# the goal of this is to use MCMCglmm and BSFG to fit a model with two correlated random effects for one trait
# The first example is the case of the same trait in two environments - a genetic correlation, but no environmental correlation

library(reshape2)
nG = 100
nE = 2
nR = 4
data = expand.grid(Geno = factor(1:nG), Env = LETTERS[1:nE], Rep = 1:nR)

geno_base = rnorm(nG)
geno_1 = geno_base + 0.5*rnorm(nG)
geno_2 = geno_base + 0.5*rnorm(nG)
data$y = rnorm(nrow(data)) + (data$Env == unique(data$Env)[1])*geno_1[data$Geno] + (data$Env == unique(data$Env)[2])*geno_2[data$Geno]
data$ID = paste(data$Geno,data$Rep,data$Env)
d = c()
for(i in 1:2){
  di = data
  di$Gene = i
  di$y = di$y + rnorm(nrow(di))
  d = rbind(d,di)
}
data = d
data$Gene = factor(data$Gene)

data2 = dcast(data,Geno+Rep+ID~Env+Gene,value.var = 'y')



# library(lme4)
# lme1 = lmer(y~Env:Gene + (0+Env:Gene|Geno),data)
# as.data.frame(VarCorr(lme1))

library(MCMCglmm)
n_trait = length(unique(data$Env:data$Gene))

# mcmc1 = MCMCglmm(y ~ Env,data=data,random = ~us(Env):Geno,rcov = ~us(Env):units,prior = list(G = list(G1=list(V=diag(1,2),nu=3)),R = list(V=diag(1,2),nu=3)))
# summary(mcmc1)
#
mcmc2 = MCMCglmm(y ~ Gene:Env,data=data,random = ~us(Env:Gene):Geno,rcov = ~idh(Env:Gene):units,
                 prior = list(G = list(G1=list(V=diag(1,n_trait),nu=n_trait+2)),R = list(V=diag(1,n_trait),nu=n_trait+1)))
summary(mcmc2)
i = c(matrix(1:n_trait,nc=nE,byrow=T))
Gcor = cov2cor(matrix(colMeans(mcmc2$VCV[,1:n_trait^2]),n_trait))[i,i]
image(Matrix(Gcor))

mcmc3 = MCMCglmm(cbind(A,B)~trait,data=data2,random = ~us(trait):Geno,rcov = ~us(trait):units,prior = list(G = list(G1=list(V=diag(1,2),nu=3)),R = list(V=diag(1,2),nu=3)),family = rep('gaussian',2))
summary(mcmc3)

mcmc4 = MCMCglmm(cbind(A,B)~trait,data=data2,random = ~us(trait):Geno,rcov = ~idh(trait):units,prior = list(G = list(G1=list(V=diag(1,2),nu=3)),R = list(V=diag(1,2),nu=3)),family = rep('gaussian',2))
summary(mcmc4)

library(BSFG)


run_parameters = BSFG_control(
  sampler = 'fast_BSFG',
  # sampler = 'general_BSFG',
  scale_Y = FALSE,
  simulation = FALSE,
  h2_divisions = 100,
  h2_step_size = 0.3,
  burn = 100,
  k_init = 2
)

priors = BSFG_priors(
  fixed_var = list(V = 1,     nu = 3),
  # tot_Y_var = list(V = 0.5,   nu = 3),
  tot_Y_var = list(V = 0.5,   nu = 10),
  tot_F_var = list(V = 18/20, nu = 20),
  h2_priors_resids_fun = function(h2s,n) 1,#pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
  h2_priors_factors_fun = function(h2s,n) ifelse(h2s == 0 | h2s >= .99,n,n/(n-1)),
  Lambda_prior = list(
    sampler = sample_Lambda_prec_ARD,
    Lambda_df = 3,
    delta_1   = list(shape = 2.1,  rate = 1/20),
    delta_2   = list(shape = 3, rate = 1)
  ),
  B_prior = list(
    sampler = sample_B_prec_ARD,
    B_df      = 3,
    B_F_df    = 3
  )
)

Y = data2[,-c(1:3)]
BSFG_state = BSFG_init(Y, model=~1+(1|Geno), data2,
                       run_parameters=run_parameters,
                       priors=priors)

BSFG_state = reorder_factors(BSFG_state)
BSFG_state = clear_Posterior(BSFG_state)
n_samples = 100;
for(i  in 1:20) {
  if(i == 10){
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
  BSFG_state = save_posterior_chunk(BSFG_state)
  print(BSFG_state)
  plot(BSFG_state)
}


BSFG_state$Posterior = reload_Posterior(BSFG_state)

Gs = get_posterior_FUN(BSFG_state,Lambda %*% diag(F_h2[1,]) %*% t(Lambda) + diag(resid_h2[1,]/tot_Eta_prec[1,]))
get_posterior_mean(Gs)
image(Matrix(cov2cor(get_posterior_mean(Gs))))
effectiveSize(Gs[,,1])

Rs = get_posterior_FUN(BSFG_state,Lambda %*% diag(1-F_h2[1,]) %*% t(Lambda) + diag((1-resid_h2[1,])/tot_Eta_prec[1,]))
get_posterior_mean(Rs)
trace_plot(cbind(Gs[,1,2],Rs[,1,2]))
trace_plot(cbind(Gs[,1,1],Rs[,1,1]))
trace_plot(mcmc3$VCV[,c(1,5),drop=F])
trace_plot(mcmc3$VCV[,c(1,5)+2,drop=F])

Gcor = get_posterior_FUN(BSFG_state,cov2cor(Lambda %*% diag(F_h2[1,]) %*% t(Lambda) + diag(resid_h2[1,]/tot_Eta_prec[1,])))
get_posterior_mean(Gcor)
plot(Gcor[,1,2])


Gs1 = get_posterior_FUN(BSFG_state,rowSums(Lambda^2 * F_h2[rep(1,nrow(Lambda)),]))
Gs2 = get_posterior_FUN(BSFG_state,resid_h2[1,]/tot_Eta_prec[1,])



run_parameters = BSFG_control(
  sampler = 'fast_BSFG',
  # sampler = 'general_BSFG',
  scale_Y = FALSE,
  simulation = FALSE,
  h2_divisions = 10,
  h2_step_size = 0.3,
  burn = 100,
  k_init = 15
)

####
observation_setup = list(
  observation_model = regression_model,
  observations = data,
  individual_model = y~0+Gene:Env,
  resid_Y_prec_shape = 2,
  resid_Y_prec_rate = 1
)



BSFG_state = BSFG_init(observation_setup, model=~1+(1|Geno), data2,
                       run_parameters=run_parameters,
                       priors=priors)

BSFG_state = reorder_factors(BSFG_state)
BSFG_state = clear_Posterior(BSFG_state)
n_samples = 100;
for(i  in 1:20) {
  if(i == 10){
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
  BSFG_state = save_posterior_chunk(BSFG_state)
  print(BSFG_state)
  plot(BSFG_state)
}
