#model_path = "~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG"
model_path="."
library(Rcpp)
library(RcppArmadillo)
source(paste(model_path,'BSFG_functions.R',sep='/'))
sourceCpp(paste(model_path,'BSFG_functions_c.cpp',sep='/'))
source(paste(model_path,'fast_BSFG_sampler_init_fixedlambda.R',sep='/'))
source(paste(model_path,'fast_BSFG_sampler_fixedlambda.R',sep='/'))

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

# For the 1st comparation: choose 5 as old BSFG_state. Pairs(5,6; 5,7)
# For the 2nd comparation: choose 5 as old BSFG_state. Pairs(5,14; 5,17)
setwd('~/../Desktop/5/Lambda1.5_delta2shape3')
# fixed lambda
# load 'BSFG_state' as the old BSFG_state
load("BSFG_state.RData")
priors = BSFG_state$priors

# YNew and YOld dataset should be in the upper directory of current dir
# 'setup datasets' will have the same file names by using setup_pedigree.R function :'setup.RData'. 
# Need to change the names of 'setup dataset' files as the example below
YNew="setup6.RData"
YOld="setup5.RData"
BSFG_state = fast_BSFG_sampler_init_fixedlambda(priors,run_parameters,YNew,YOld)

# load 'BSFG_fixedlambda' as the new BSFG_state
BSFG_state = fast_BSFG_sampler_fixedlambda(BSFG_state,YNew,YOld)
load("BSFG_fixedlambda.RData")
