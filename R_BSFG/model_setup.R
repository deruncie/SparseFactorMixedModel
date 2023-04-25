
setwd("~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG")
#model_path = '~/Box Sync/DER_projects/BSFG/R_BSFG'
#model_path = "../SparseFactorMixedModel/R_BSFG"
#model_path = "~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG"

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



# set the directory to the location of the setup.RData or setup.mat file
#setwd('Example_simulation')
setwd('1')
# choose a seed for the random number generator. This can be a random seed (for analysis), or you can choose your seed so that 
# you can repeat the MCMC exactly
#seed = sample(1:1e3,1)
#set.seed(seed)

# create a folder for holding the posterior samples of the current chain (multiple folders could be used for different chains)
folder = "Lambda1.5_delta2shape3"   #original
#rep = "fixedlambda1"
rep = "test666"
folder = "generation0"
folder = sprintf('R_rep_%s',rep)
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

#BSFG_state = clear_Posterior(BSFG_state)


# # optional: To load from end of previous run, run above code, then run these lines:
# load('Posterior.RData')
# load('BSFG_state.RData')
# load('Priors.RData')
#  BSFG_state$current_state = current_state
#  BSFG_state$Posterior = Posterior
#  BSFG_state$priors = priors
# start_i = run_parameters$nrun;

# Run Gibbs sampler. Run in smallish chunks. Output can be used to re-start chain where it left off.
# c1=fix(clock);
n_samples = 4000;
for(i  in 1:1) {
    print(sprintf('Run %d',i))
    BSFG_state = fast_BSFG_sampler_current(BSFG_state,n_samples)
    BSFG_state = fast_BSFG_sampler_generation(BSFG_state,n_samples)
    print(i)
} 

## function only for plotting
setwd('~/Runcie Lab/SparseFactorMixedModel/R_BSFG/R_rep/R_tune_5')
load("BSFG_state.RData")
plot_diagnostics(BSFG_state)

#fixed lambda
#load("BSFG_state.RData")
#priors = BSFG_state$priors
BSFG_state = fast_BSFG_sampler_init_fixedlambda(priors,run_parameters,YNew="setup_LineCode2.RData",YOld="setup_LineCode1.RData")

load("BSFG_fixedlambda.RData")
BSFG_state = fast_BSFG_sampler_fixedlambda(BSFG_state)



setwd('~/Runcie Lab/SparseFactorMixedModel/R_BSFG/Sim_1/R_rep_test')
ComparingGMatrix_plot("F_a")

setwd('~/../Desktop')
load('BSFG_state.RData')

data = read.csv('ScentDataAll (2).csv')
setup = setup_pedigree(data,3)
