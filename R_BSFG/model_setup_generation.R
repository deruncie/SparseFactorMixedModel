setwd("~/Runcie Lab/pedigree")

data=read.csv('IDDammSirelnUpdate(1).csv',header = TRUE)
source("~/Runcie Lab/pedigree/setup_pedigree_GeneLine.R")
setup=setup_pedigree_GeneLine(data,GenerationCode=0,LineCode=c(1,2,3))


setwd("~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG")
model_path = "~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG"
source(paste(model_path,'fast_BSFG_sampler_init.R',sep='/'))
source(paste(model_path,'fast_BSFG_sampler_current.R',sep='/'))
source(paste(model_path,'BSFG_functions.R',sep='/'))
source(paste(model_path,'plotting_diagnostics.R',sep='/'))
source(paste(model_path,'setup_pedigree.R',sep='/'))
sourceCpp(paste(model_path,'BSFG_functions_c.cpp',sep='/'))
source(paste(model_path,'fast_BSFG_sampler_init_fixedlambda.R',sep='/'))
source(paste(model_path,'fast_BSFG_sampler_fixedlambda.R',sep='/'))
source(paste(model_path,'fast_BSFG_sampler_generation.R',sep='/'))

library(Rcpp)
library(RcppArmadillo)

setwd('Sim_2')
folder = "generation0"
try(dir.create(folder))
setwd(folder)

# initialize priors
run_parameters = list(
  b0           = 1,
  b1           = 0.0005,
  epsilon      = 1e-2,
  prop         = 1.00,
  h2_divisions = 50,
  save_freq    = 100,
  burn = 2000,       #100
  thin = 400,         #2
  draw_iter = 200
)

#load("./priors.RData")
priors = list(
  k_init              =   20,
  resid_Y_prec_shape  =   2,
  resid_Y_prec_rate   =   1/10,
  E_a_prec_shape      =   2,
  E_a_prec_rate       =   1/10,
  W_prec_shape        =   2,
  W_prec_rate         =   1/10,
  Lambda_df           =   1.5,         # 3/2, 3, 1
  delta_1_shape       =   2.1,
  delta_1_rate        =   1/20,
  delta_2_shape       =   3,           # 3, 3/2, 9/2
  delta_2_rate        =   1,
  h2_priors_factors   =   c(run_parameters$h2_divisions-1,rep(1,run_parameters$h2_divisions-1))/(2*(run_parameters$h2_divisions-1))
)

print('Initializing')
save(priors,file = 'Priors.RData')

# Initialize Chain, prep runs
BSFG_state = fast_BSFG_sampler_init(priors,run_parameters)
save(BSFG_state,file="BSFG_state.RData")

n_samples = 4000;
for(i  in 1:10) {
  print(sprintf('Run %d',i))
  BSFG_state = fast_BSFG_sampler_generation(BSFG_state,n_samples)
  print(i)
} 
