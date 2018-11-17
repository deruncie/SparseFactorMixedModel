library(microbenchmark)
library(MCMCpack)
library(BSFG)
library(cowplot)
library(splines)
library(Matrix)

# set the directory to the location of the setup.RData or setup.mat file

#
# # # choose a seed for the random number generator. This can be a random seed (for analysis), or you can choose your seed so that
# # # you can repeat the MCMC exactly
seed = 1
set.seed(seed)
p=4
new_halfSib_spline_simulation('Sim_FE_1', nSire=2,nRep=200,p=p, Time = seq(1,60,length=10), k=4, k_G=4)#, factor_h2s = rep(0.9,4),resid_h2 = c(.3,.3,0,0))
load('setup.RData')
#
ggplot(setup$observations,aes(x=covariate,y=Y)) + geom_line(aes(group = ID))
#
# create a folder for holding the posterior samples of the current chain (multiple folders could be used for different chains)
rep = "2"
folder = sprintf('R_spline_rep_%s',rep)
try(dir.create(folder))
setwd(folder)

# spline simulation
n = 500
degree=4
range = 10
Times = seq(1-10,60+10,length=60)
data_ID = data.frame(ID=1:(n*2),animal = 1:(n*2),group = rep(1:2,each=n))
data_tall = list()
Terms = delete.response(terms(model.frame(~poly(Times,degree=degree))))
b = rnorm(degree)*3
for(id in 1:n){
  shift = 0*runif(1,0,30)
  my_Times = sample(1:60,10)
  my_Times = seq(1,60,length=60) #+ runif(1,-5,5)/2
  X = model.matrix(Terms,data=data.frame(Times = my_Times-shift))
  # y = X %*% c(0,b) + rnorm(length(my_Times),0,.1)
  y = X %*% c(0,b+.1*rnorm(degree)*3) + rnorm(length(my_Times),0,.1)
  data_tall[[id]] = data.frame(ID = id,covariate = my_Times,Y=y,group = 1)
}
b = b+.5*rnorm(degree,0,c(0,0,rep(1,degree-2)))*3
for(id in n+1:(n)){
  shift = 0*runif(1,0,30)
  my_Times = sample(1:60,10)
  my_Times = seq(1,60,length=60) #+ runif(1,-5,5)/2
  X = model.matrix(Terms,data=data.frame(Times = my_Times-shift))
  # y = X %*% c(0,b) + rnorm(length(my_Times),0,.1)
  y = X %*% c(0,b+.1*rnorm(degree)*3) + rnorm(length(my_Times),0,.1)
  data_tall[[id]] = data.frame(ID = id,covariate = my_Times,Y=y,group = 2)
}
data_tall = do.call(rbind,data_tall)
data_tall = subset(data_tall,covariate > 10 & covariate < 55)
ggplot(data_tall,aes(x=covariate,y=Y,group=ID)) + geom_line(aes(color=group))

K = diag(1,2*n)
rownames(K) = 1:(2*n)
setup = list(data = data_ID,observations = data_tall,K = K)


library(rrBLUP)
mm= model.matrix(~0+bs(covariate,df=30,intercept=F),data_tall)
D = matrix(0,ncol(mm),ncol(mm))
D[lower.tri(D)] = 1
diag(D) = 1
coefs = tapply(1:nrow(data_tall),data_tall$ID,function(x) {
  sub_data = data_tall[x,]
  X = mm[x,]
  XD = X %*% D
  res = mixed.solve(y = sub_data$Y,Z = XD,X = matrix(1,nrow = nrow(X),nc=1))
  beta = c(res$beta[1],res$u)
  return(beta)
})
coefs = do.call(rbind,coefs)



# initialize priors
run_parameters = BSFG_control(
  # sampler = 'fast_BSFG',
  # sampler = 'general_BSFG',
  scale_Y = FALSE,
  simulation = FALSE,
  h2_divisions = 20,
  h2_step_size = NULL,
  burn = 1000,
  thin=10,
  K = 10
)

priors = BSFG_priors(
  tot_Y_var = list(V = 1,   nu = 3),
  tot_F_var = list(V = 18/20, nu = 20e6),
  h2_priors_resids_fun = function(h2s,n) 1,#pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
  h2_priors_factors_fun = function(h2s,n) 1#ifelse(h2s == 0 | h2s >= .99,n,n/(n-1))
  ,Lambda_prior = list(
    sampler = sample_Lambda_prec_horseshoe,
    prop_0 = 0.5,
    delta_l = list(shape = 3, rate = 1),
    delta_iterations_factor = 100
  ),
  B2_prior = list(
    sampler = sample_B2_prec_horseshoe,
    prop_0 = 0.1
  ),
  cis_effects_prior = list(
    prec = 1e-6
  )
)


print('Initializing')
# setup = load_simulation_data()
# setup = load_simulation_FE_data()
# load('../setup.RData')


# data_model_parameters = list(
#   observations = setup$observations,
#   df = 40,
#   # degree=6,
#   intercept = TRUE,
#   individual_model = ~SPLINE,
#   resid_Y_prec_shape = 2,
#   resid_Y_prec_rate = 1
# )
setup$observations$A = setup$observations$Y
setup$observations$B = setup$observations$Y
setup$observations$C = setup$observations$Y
observation_setup = list(
  observation_model = regression_model,
  observations = setup$observations,
  # individual_model = Y~ poly(covariate,4)+bs(covariate,df=5,intercept=T),
  # individual_model = Y~ poly(covariate,4)+bs(covariate,df=5,intercept=T)+bs(covariate,df=20,intercept=T)+bs(covariate,df=40,intercept=F),
  # individual_model = cbind(A,B,C)~poly(covariate,4)+bs(covariate,df=5,intercept=T)+bs(covariate,df=20,intercept=T)+bs(covariate,df=40,intercept=F),
  individual_model = A~1+bs(covariate,df=30,intercept=F),
  resid_Y_prec_shape = 2,
  resid_Y_prec_rate = 1,
  do_not_penalize_bs = TRUE
)

setup$data$ID = unique(setup$observations$ID)
# observation_setup$observation_setup = regression_model(observation_setup,list(data_matrices=list(data=setup$data)))
#
# a = regression_model(observation_setup,list(data_matrices=list(data=setup$data)))


# options(error=recover)
data = setup$data
BSFG_state = setup_model_BSFG(observation_setup,formula = ~group+(1|animal),
                              data = data,
                              run_parameters = run_parameters,
                              run_ID = 'spline_1')
BSFG_state = set_priors_BSFG(BSFG_state,priors)
BSFG_state = initialize_variables_BSFG(BSFG_state)
BSFG_state = initialize_BSFG(BSFG_state)
BSFG_state$Posterior$posteriorSample_params = c(BSFG_state$Posterior$posteriorSample_params,'Lambda_c2')
BSFG_state = clear_Posterior(BSFG_state)

# BSFG_state = with(setup,BSFG_init(observation_setup, model=~group+(1|animal), data, #factor_model_fixed = ~1,
# # BSFG_state = with(setup,BSFG_init(Eta, model=~1+(1|animal), data, #factor_model_fixed = ~1,
#                                   priors=priors,run_parameters=run_parameters,K_mats = list(animal = K),
#                                   posteriorSample_params = c('Lambda','U_F','F','delta','tot_F_prec','F_h2','tot_Eta_prec','resid_h2', 'B', 'B_F', 'B_tau','B_F_tau','cis_effects','U_R'),
#                                   posteriorMean_params = c()
#                                   ))
# BSFG_state$current_state$F_h2
# BSFG_state$priors$h2_priors_resids
# BSFG_state$priors$h2_priors_factors

# save(BSFG_state,file="BSFG_state.RData")
#
#
# Eta = BSFG_state$current_state$Eta
# BSFG_state = initialize_variables(BSFG_state)
# BSFG_state$current_state$Eta = Eta

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
  var_Eta_factor = apply(BSFG_state$current_state$Eta,2,var)
  # var_Eta_factor[] = mean(var_Eta_factor)
  print(var_Eta_factor)
  print(apply(BSFG_state$current_state$F,2,var))
  if(i %%5 == 0 && i < 20) {
    # BSFG_state$current_state$var_Eta = var_Eta_factor
    rescale_Eta = function(BSFG_state, var_Eta_factor){
      current_state = within(BSFG_state$current_state,{
        var_Eta = var_Eta * var_Eta_factor
        Eta = sweep(Eta,2,sqrt(var_Eta_factor),'/')
        Lambda = sweep(Lambda,1,sqrt(var_Eta_factor),'/')
        # Plam = sweep(Plam,1,sqrt(var_Eta_factor),'*')
        # B = sweep(B,2,sqrt(var_Eta_factor),'/')
        U_R = sweep(U_R,2,sqrt(var_Eta_factor),'/')
        tot_Eta_prec = tot_Eta_prec * var_Eta_factor
        delta = delta * sqrt(mean(var_Eta_factor))
      })
      current_state
    }
    var_Eta_factor = c(1,rep(median(var_Eta_factor[-1]),length(var_Eta_factor)-1))
    BSFG_state$current_state = rescale_Eta(BSFG_state,var_Eta_factor)
  }
    print(sprintf('Run %d',i))
    BSFG_state = sample_BSFG(BSFG_state,n_samples,ncores=1)
    if(BSFG_state$current_state$nrun < BSFG_state$run_parameters$burn) {
      BSFG_state = reorder_factors(BSFG_state)
      BSFG_state = clear_Posterior(BSFG_state)
    }
    BSFG_state = save_posterior_chunk(BSFG_state)
    print(BSFG_state)
    plot(BSFG_state)
    par(mfrow=c(1,2))
    plot(BSFG_state$run_parameters$observation_model_parameters$observation_setup$Y,BSFG_state$current_state$Y_fitted);abline(0,1)
    boxplot(BSFG_state$current_state$Eta)
}

observation_model_state = with(BSFG_state,run_parameters$observation_model(run_parameters$observation_model_parameters,BSFG_state)$state)
BSFG_state$current_state[names(observation_model_state)] = observation_model_state

Posterior = reload_Posterior(BSFG_state)
data = setup$observations
trait = 1
data$Y_fitted = colMeans(Posterior$Y_fitted[,,trait])
data$Y_fitted_low = apply(Posterior$Y_fitted[,,trait],2,function(x) HPDinterval(mcmc(x))[1])
data$Y_fitted_high = apply(Posterior$Y_fitted[,,trait],2,function(x) HPDinterval(mcmc(x))[2])
# data$Y_fitted = BSFG_state$current_state$Y_fitted
ggplot(data,aes(x=covariate,y=Y)) + geom_line(aes(group = ID)) + geom_line(aes(y=Y_fitted,group=ID),color ='red')
ggplot(subset(data,ID %in% sample(unique(data$ID),1)),aes(x=covariate,y=Y)) + geom_ribbon(aes(ymin = Y_fitted_low,ymax = Y_fitted_high,group = ID),alpha = 0.2) + geom_line(aes(group = ID)) + geom_line(aes(y=Y_fitted,group=ID),color ='red')

Lambda = get_posterior_mean(BSFG_state,Lambda)

newx = seq(15,50,length=100)
new_MM = model.matrix(BSFG_state$run_parameters$observation_model_parameters$observation_setup$Terms,data.frame(covariate=newx)) #%*% BSFG_state$run_parameters$observation_model_parameters$observation_setup$mm_rotation
# new_MM = model.matrix(BSFG_state$run_parameters$observation_model_parameters$observation_setup$Terms,data.frame(Time=newx))
new_MM = sweep(new_MM,2,sqrt(BSFG_state$current_state$var_Eta),'*')
a=new_MM %*% t(BSFG_state$current_state$Eta)
ids = sample(1:ncol(a),1)
plot(NA,NA,xlim = range(newx),ylim = range(a))
aa=sapply(ids,function(i) lines(newx,a[,i]))
with(subset(data,ID %in% ids),points(covariate,Y))



Eta = BSFG_state$current_state$Eta[1:10,]
plot(NA,NA,xlim = c(1,ncol(Eta)),ylim = range(Eta))
aa=sapply(1:nrow(Eta), function(i) lines(Eta[i,]))

Image(cor(Eta))
Image(cov2cor(tcrossprod(BSFG_state$current_state$Lambda)+diag(1/BSFG_state$current_state$tot_Eta_prec[1,])))

BSFG_state$Posterior = reload_Posterior(BSFG_state)
plot(apply(BSFG_state$Posterior$Eta,c(2,3),mean),setup$Eta);abline(0,1)
# setup$observations$Y_fitted = sapply(1:nrow(setup$observations),function(i) predict(BSFG_state$current_state$coefficients,newx = setup$observations$covariate[i]) %*% BSFG_state$Posterior$Eta[setup$observations$ID[i],])
setup$observations$Y_fitted = sapply(1:nrow(setup$observations),function(i) model.matrix(BSFG_state$run_parameters$observation_model_parameters$observation_setup$Terms,data.frame(covariate=setup$observations$covariate[i])) %*% colMeans(BSFG_state$Posterior$Eta[,setup$observations$ID[i],]))
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
dimnames(Posterior_P)[2:3] = dimnames(Posterior$Eta)[3]

library(heatmap3)
i = 1:dim(Posterior_P)[2]
heatmap3(cov2cor(apply(Posterior_P,c(2,3),mean))[i,i]^2,Rowv = NA,Colv=NA)
trace_plot(Posterior_P[,i,1])

Posterior_G = aperm(array(sapply(1:Posterior$total_samples,function(i) {
  with(Posterior,{
    Lambda[i,,] %*% diag(F_h2[i,,]) %*% t(Lambda[i,,]) + diag(resid_h2[i,,]/tot_Eta_prec[i,,])
  })
}),dim = c(rep(dim(Posterior$Lambda)[2],2),Posterior$total_samples)),c(3,1,2))
dimnames(Posterior_G)[2:3] = dimnames(Posterior$Eta)[3]
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
