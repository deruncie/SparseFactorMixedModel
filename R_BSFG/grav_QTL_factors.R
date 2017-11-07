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
# X = do.call(cbind,lapply(grav$geno,function(x) x$prob[,,1]))
X = do.call(cbind,lapply(grav$geno,function(x) x$data-1))
Chr = do.call(c,sapply(1:length(grav$geno),function(x) rep(sprintf('Chr%d',x),ncol(grav$geno[[x]]$data))))
pos = do.call(c,lapply(grav$geno,function(x) x$map))
# X = do.call(cbind,lapply(grav$geno,function(x) apply(x$draws,c(1,2),mean)))
# Chr = do.call(c,sapply(1:length(grav$geno),function(x) rep(sprintf('Chr%d',x),ncol(grav$geno[[x]]$draws))))

y = unlist(Y[,ncol(Y)])
data = data.frame(ID = 1:nrow(Y),Y = (y-mean(y))/sd(y))

tall_data = melt(X)
colnames(tall_data)[1:2] = c('ID','Pos')
tall_data$Pos = rep(pos,each = nrow(X))
tall_data$Chr = rep(Chr,each = nrow(X))
tall_data = na.omit(tall_data)


# initialize priors
run_parameters = BSFG_control(
  scale_Y = TRUE,
  simulation = FALSE,
  h2_divisions = 100,
  h2_step_size = NULL,
  burn = 10000,
  thin = 10,
  k_init = 12
)

priors = BSFG_priors(
  fixed_var = list(V = 1,     nu = 3),
  # QTL_resid_var = list(V = 1/10000,     nu = 3*1000),
  # QTL_factors_var = list(V = 1/10000,     nu = 3),
  # tot_Y_var = list(V = 0.5,   nu = 3),
  tot_Y_var = list(V = 0.5,   nu = 3),
  tot_F_var = list(V = 18/20, nu = 20e5),
  h2_priors_resids_fun = function(h2s,n) 1,#pmax(pmin(ddirichlet(c(h2s,1-sum(h2s)),rep(2,length(h2s)+1)),10),1e-10),
  h2_priors_factors_fun = function(h2s,n) 1,#ifelse(h2s == 0 | h2s >= .99,n,n/(n-1)),
  Lambda_prior = list(
    sampler = sample_Lambda_prec_ARD,
    Lambda_df = 3,
    delta_1   = list(shape = 2.1,  rate = 1/20),
    delta_2   = list(shape = 1.5, rate = 1)
  ),
  B_prior = list(
    sampler = sample_B_prec_ARD,
    B_df      = 3,
    B_F_df    = 2
  )
)

observation_setup = list(
  observation_model = regression_model,
  observations = tall_data,
  individual_model = value~0+Chr:bs(Pos,df=20),
  # individual_model = value~0+DaylTemp + DaylTemp:(Vern*oscillation),
  resid_Y_prec_shape = 2,
  resid_Y_prec_rate = 1
)
a=regression_model(observation_setup,list(data_matrices = list(data = data)))
# Y = a$state$Eta
# Y[a$state$Y_missing] = NA
#
# colnames(Y) = sub('value::variable','',colnames(Y))
# Y = Y[,order(colnames(Y))]

# Y[] = NA
# Y_temp = Y
# Y[is.na(Y)] = rnorm(sum(is.na(Y)))
BSFG_state = BSFG_init(observation_setup, model=~Y+(1|ID), data, factor_model_fixed = ~Y,#+(1|Treatment.chamber)
                       # BSFG_state = BSFG_init(Y, model=~1+(1|genotype), wide_data,
                       priors=priors,run_parameters=run_parameters)

n_samples = 100;
for(i  in 1:110) {
  var_Eta_factor = apply(BSFG_state$current_state$Eta,2,var)
  print(var_Eta_factor)
  if(i == 8) {
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
    BSFG_state$current_state = rescale_Eta(BSFG_state,var_Eta_factor)
    BSFG_state$run_parameters$lambda_propto_Vp = FALSE
  }


  # if(i < 10 || (i-1) %% 20 == 0) {
  if(i < 10) {
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
  plot(BSFG_state$current_state$Y,BSFG_state$current_state$Y_fitted);abline(0,1)
}
BSFG_state$Posterior = reload_Posterior(BSFG_state)

mm_rotation = a$mm_rotation
vE = BSFG_state$current_state$var_Eta

newdata = data.frame(Chr = Chr,Pos = pos)
MM = with(BSFG_state$run_parameters$observation_model_parameters$observation_setup,model.matrix(Terms,newdata)[,a$nonzero_cols] %*% mm_rotation)
MM = sweep(MM,2,sqrt(vE),'*')

plot(MM %*% BSFG_state$current_state$Lambda[,1])
plot(MM %*% t(BSFG_state$current_state$Eta[1,,drop=F]))

L = get_posterior_mean(BSFG_state,Lambda)
plot(MM %*% L[,1])

y_eff = get_posterior_mean(BSFG_state,(B_F %*% t(Lambda) + B[2,,drop=F]) %*% t(MM))
plot(y_eff[2,])
scanone_pm = scanone(grav, pheno.col = ncol(Y),method="hk")
plot(scanone_pm)
