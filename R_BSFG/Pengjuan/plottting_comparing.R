setwd("~/Runcie Lab/SparseFactorMixedModel_v2/MCMC/5/5/Lambda1.5_delta2shape3")
BSFGOld = "BSFG_state.RData"
BSFGNew = "BSFG_fixedlambda65.RData"
BSFGStar = "BSFG_fixedlambda75.RData"
target = "F_h2"
#it looks like we do not need so much latent factors. choose the first l factors.
#choose the first l latent traits and make them a vector for all iterations 
l=5  # 2:6
pos_fa = function(spn,pos.fa,n){
  pos.fa.rb = vector()
  for (i in 1:spn){
    pos.fa.rb = c(pos.fa.rb,pos.fa[1:(l*n),i])
  }
  return(pos.fa.rb)}
# old
load(BSFGOld)
spn = dim(BSFG_state$Posterior[[target]])[2]
old.n   = dim(BSFG_state$data_matrices$Y)[1]
old.k   = nrow(BSFG_state$Posterior[[target]])
old.pos.fa = BSFG_state$Posterior$F_a

old.pos.fa.rb = pos_fa(spn,old.pos.fa,old.n)
# new
load(BSFGNew)
spn = dim(BSFG_state$Posterior[[target]])[2]
new.n   = dim(BSFG_state$data_matrices$Y)[1]
new.k   = nrow(BSFG_state$Posterior[[target]])
new.pos.fa = BSFG_state$Posterior$F_a
new.pos.fa.rb = pos_fa(spn,new.pos.fa,new.n)
# star : represent the third comparison
load(BSFGStar)
spn = dim(BSFG_state$Posterior[[target]])[2]
star.n   = dim(BSFG_state$data_matrices$Y)[1]
star.k   = nrow(BSFG_state$Posterior[[target]])
star.pos.fa = BSFG_state$Posterior$F_a
star.pos.fa.rb = pos_fa(spn,star.pos.fa,star.n)
# combine
# column "old", "new", "star"
pos.fa.b = cbind(c(old.pos.fa.rb,new.pos.fa.rb,star.pos.fa.rb),rep(c("old","new","star"),c(old.n*l*spn,new.n*l*spn,star.n*l*spn)))
# column "factor"
pos.fa.b = cbind(pos.fa.b,c(rep(rep(sprintf("f%d",1:l),each=old.n),spn),rep(rep(sprintf("f%d",1:l),each=new.n),spn),rep(rep(sprintf("f%d",1:l),each=star.n),spn)))
colnames(pos.fa.b) = c("Fa","model", "factor")
pos.fa.b=as.data.frame(pos.fa.b)
#save(pos.fa.b,file="pos_fa_bind.RData")
load("pos_fa_bind.RData")
# plot
# 1.boxplot of distribution of Fa
library(ggplot2)
ggplot(pos.fa.b, aes(x = factor, y = as.numeric(Fa), color=model,fill=model))+
  geom_boxplot()
  #geom_tile("The distribution of Latent traits")
  #geom_abline(intercept=0,slope=0,lty=2)

ggplot(pos.fa.b)+
  geom_density(aes(x = factor, y = as.numeric(Fa), color=model))
  #geom_tile("The distribution of Latent traits")
# 2.For each trait, compare the correlation of G1 vs G2
load(BSFGOld)
old.traces.G = G_Traces_Comp(BSFG_state)
old.G = G_Matrix_Comp(BSFG_state)
p = BSFG_state$run_variables$p
# do this in case the value of k differ for this two models and there in no lambda in BSFG_state_fixed case
Lambda = BSFG_state$Posterior$Lambda
load(BSFGNew)
BSFG_state$Posterior$Lambda = Lambda
new.traces.G = G_Traces_Comp(BSFG_state)
new.G = G_Matrix_Comp(BSFG_state)
#For all the elements in G
#pairs(t(old.traces.G)~t(new.traces.G))
#For each trait
old.G = array(unlist(old.G),dim = c(nrow(old.G[[1]]),ncol(old.G[[1]]),length(old.G)))
new.G = array(unlist(new.G),dim = c(nrow(new.G[[1]]),ncol(new.G[[1]]),length(new.G)))
old.G.b = apply(old.G,2,rbind)
new.G.b = apply(new.G,2,rbind)
#plot(old.G.b[,2],new.G.b[,2])
#combine old.G and new.G as a dataframe
#G.b = c(old.G.b[,2],new.G.b[,2])
#G.b = cbind(G.b, rep(c("old","new"),each=p*spn))
#G.b = as.data.frame(G.b)
#colnames(G.b) = c("covariance","model")
# 1st triat
each_trait_G_pair = function(i){
G.b = cbind(old.G.b[,i],new.G.b[,i])
G.b = as.data.frame(G.b)
colnames(G.b) = c("Old","New")
library(ggplot2)
ggplot(G.b,aes(y = as.numeric(Old),x=as.numeric(New)))+
  #geom_dotplot(binwidth = 15)
  geom_point()+
  geom_abline(intercept = 0,slope = 1)
}
each_trait_G_pair(1)
each_trait_G_pair(2)
each_trait_G_pair(3)
# 3.delta z
# The distance of two matrix(in predicted phenotype value/level)

