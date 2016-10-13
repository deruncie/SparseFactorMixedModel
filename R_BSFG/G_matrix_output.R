library(Rcpp)
library(RcppArmadillo)
model_path = "~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG"
data = read.csv('~/Runcie Lab/pedigree/IDDammSirelnUpdate2.csv',header = TRUE)
source(paste(model_path,'plotting_diagnostics.R',sep='/'))
source(paste(model_path,'setup_pedigree.R',sep='/'))

#-----------------------------------------------------
# Second for loop
setwd('~/Runcie Lab/SparseFactorMixedModel_v2/MCMC')
G_BSFG_post = list()
G_BSFG_sc = list()
G_BSFG = list()

for (i in c(1:9,14:19)) {
setwd(sprintf('~/Runcie Lab/SparseFactorMixedModel_v2/MCMC/%d/%d',i,i))
  load("setup.RData")
  Y = setup$Y
  VY = apply(Y,2,var,na.rm=T)
  
setwd("Lambda1.5_delta2shape3")
load("BSFG_state.RData")
# Scaled G Matrix
G_BSFG_sc[[i]] = G_Matrix_Comp(BSFG_state)
G_est_all = G_BSFG_sc[[i]]
sp_num = length(G_est_all)
p = BSFG_state$run_variables$P

traitnames = BSFG_state$traitnames


#setwd(paste0(i))
#save(G_est_all, file = "G_est_all.RData")

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
#names(G_BSFG) = names(G_MCMC)
#names(G_BSFG_sc) = names(G_MCMC)
#names(G_BSFG_post) = names(G_MCMC)
save(G_BSFG,file = "G_BSFG.RData")
save(G_BSFG_sc,file = "G_BSFG_sc.RData")
save(G_BSFG_post,file = "G_BSFG_post.RData")
