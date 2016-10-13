library(Rcpp)
library(RcppArmadillo)
model_path = "~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG"
data = read.csv('~/Runcie Lab/pedigree/IDDammSirelnUpdate2.csv',header = TRUE)
source(paste(model_path,'plotting_diagnostics.R',sep='/'))
source(paste(model_path,'setup_pedigree.R',sep='/'))

# First get the setup dataset
# 1
setwd('~/Runcie Lab/SparseFactorMixedModel_v2/MCMC')
GenerationCode=c(0,1,2,3)
LineCode = c(0,1,2,3)
fixedfactor="Generation,Line"

folder = "1"
try(dir.create(folder))
setwd(folder)
setup=setup_pedigree(data=data,GenerationCode,LineCode,fixedfactor,randomfactor=FALSE)

# 2
setwd('~/Runcie Lab/SparseFactorMixedModel_v2/MCMC')
GenerationCode=c(0,1,2,3)
LineCode = c(0,1)
fixedfactor="Generation"

folder = "2"
try(dir.create(folder))
setwd(folder)
setup=setup_pedigree(data=data,GenerationCode,LineCode,fixedfactor,randomfactor=FALSE)

# 3
setwd('~/Runcie Lab/SparseFactorMixedModel_v2/MCMC')
GenerationCode=c(0,1,2,3)
LineCode = c(0,2)
fixedfactor="Generation"

folder = "3"
try(dir.create(folder))
setwd(folder)
setup=setup_pedigree(data=data,GenerationCode,LineCode,fixedfactor,randomfactor=FALSE)

# 4
setwd('~/Runcie Lab/SparseFactorMixedModel_v2/MCMC')
GenerationCode=c(0,1,2,3)
LineCode = c(0,3)
fixedfactor= "Generation"

folder = "4"
try(dir.create(folder))
setwd(folder)
setup=setup_pedigree(data=data,GenerationCode,LineCode,fixedfactor,randomfactor=FALSE)

# 5
setwd('~/Runcie Lab/SparseFactorMixedModel_v2/MCMC')
GenerationCode=3
LineCode = c(0,1)
fixedfactor= NA

folder = "5"
try(dir.create(folder))
setwd(folder)
setup=setup_pedigree(data=data,GenerationCode,LineCode,fixedfactor,randomfactor=FALSE)

# 6
setwd('~/Runcie Lab/SparseFactorMixedModel_v2/MCMC')
GenerationCode=3
LineCode = c(0,2)
fixedfactor= NA

folder = "6"
try(dir.create(folder))
setwd(folder)
setup=setup_pedigree(data=data,GenerationCode,LineCode,fixedfactor,randomfactor=FALSE)

# 7
setwd('~/Runcie Lab/SparseFactorMixedModel_v2/MCMC')
GenerationCode=3
LineCode = c(0,3)
fixedfactor= NA

folder = "7"
try(dir.create(folder))
setwd(folder)
setup=setup_pedigree(data=data,GenerationCode,LineCode,fixedfactor,randomfactor=FALSE)

# 8
setwd('~/Runcie Lab/SparseFactorMixedModel_v2/MCMC')
GenerationCode=c(0,1)
LineCode = c(0,1,2,3)
fixedfactor= "Generation"

folder = "8"
try(dir.create(folder))
setwd(folder)
setup=setup_pedigree(data=data,GenerationCode,LineCode,fixedfactor,randomfactor=FALSE)

# 9
setwd('~/Runcie Lab/SparseFactorMixedModel_v2/MCMC')
GenerationCode=c(0,1)
LineCode = c(0,1,2,3)
fixedfactor= "Generation,Line"

folder = "9"
try(dir.create(folder))
setwd(folder)
setup=setup_pedigree(data=data,GenerationCode,LineCode,fixedfactor,randomfactor=FALSE)

# 14
setwd('~/Runcie Lab/SparseFactorMixedModel_v2/MCMC')
GenerationCode=2
LineCode = 1
fixedfactor= NA

folder = "14"
try(dir.create(folder))
setwd(folder)
setup=setup_pedigree(data=data,GenerationCode,LineCode,fixedfactor,randomfactor=FALSE)

# 15
setwd('~/Runcie Lab/SparseFactorMixedModel_v2/MCMC')
GenerationCode=2
LineCode = 2
fixedfactor= NA

folder = "15"
try(dir.create(folder))
setwd(folder)
setup=setup_pedigree(data=data,GenerationCode,LineCode,fixedfactor,randomfactor=FALSE)

# 16
setwd('~/Runcie Lab/SparseFactorMixedModel_v2/MCMC')
GenerationCode=2
LineCode = 3
fixedfactor= NA

folder = "16"
try(dir.create(folder))
setwd(folder)
setup=setup_pedigree(data=data,GenerationCode,LineCode,fixedfactor,randomfactor=FALSE)

# 17
setwd('~/Runcie Lab/SparseFactorMixedModel_v2/MCMC')
GenerationCode=1
LineCode = 1
fixedfactor= NA

folder = "17"
try(dir.create(folder))
setwd(folder)
setup=setup_pedigree(data=data,GenerationCode,LineCode,fixedfactor,randomfactor=FALSE)

# 18
setwd('~/Runcie Lab/SparseFactorMixedModel_v2/MCMC')
GenerationCode=1
LineCode = 2
fixedfactor= NA

folder = "18"
try(dir.create(folder))
setwd(folder)
setup=setup_pedigree(data=data,GenerationCode,LineCode,fixedfactor,randomfactor=FALSE)

# 19
setwd('~/Runcie Lab/SparseFactorMixedModel_v2/MCMC')
GenerationCode=1
LineCode = 3
fixedfactor= NA

folder = "19"
try(dir.create(folder))
setwd(folder)
setup=setup_pedigree(data=data,GenerationCode,LineCode,fixedfactor,randomfactor=FALSE)

#-----------------------------------------------------
# Second for loop
setwd('~/Runcie Lab/SparseFactorMixedModel_v2/MCMC')
G_MCMC = readRDS("G_MCMC.Rdata")

G_BSFG_post = list()
G_BSFG_sc = list()
G_BSFG = list()

for (i in 1:15){
setwd('~/Runcie Lab/SparseFactorMixedModel_v2/MCMC')
G_est_all = G_MCMC[i][[1]]
traitnames = rownames(G_est_all[[1]])
G_BSFG_sc[[i]] = G_est_all
sp_num = length(G_est_all) 
p = ncol(G_est_all[[1]])
setwd(paste0(i))
#save(G_est_all, file = "G_est_all.RData")
load("setup.RData")
Y = setup$Y
VY = apply(Y,2,var,na.rm=T)

#posterior mean
G_est = sqrt(diag(VY))%*%(Reduce("+",G_est_all)/sp_num)%*%t(sqrt(diag(VY)))
#save(G_est,file = "G_est_posterior.RData")
colnames(G_est) = traitnames
rownames(G_est) = traitnames
G_BSFG_post[[i]] = G_est

#unscale G_est_all
G_est_un_all = list()
for(j in 1:sp_num) {
  G_est_un_all[[j]] = sqrt(diag(VY))%*%G_est_all[[j]]%*%t(sqrt(diag(VY)))
  colnames(G_est_un_all[[j]]) = traitnames
  rownames(G_est_un_all[[j]]) = traitnames
}
G_BSFG[[i]] = G_est_un_all
}
names(G_BSFG) = names(G_MCMC)
names(G_BSFG_sc) = names(G_MCMC)
names(G_BSFG_post) = names(G_MCMC)
save(G_BSFG,file = "G_BSFG.RData")
save(G_BSFG_sc,file = "G_BSFG_sc.RData")
save(G_BSFG_post,file = "G_BSFG_post.RData")
