setwd("~/Runcie Lab/SparseFactorMixedModel_v2/MCMC/5/5/Lambda1.5_delta2shape3")
BSFGOld = "BSFG_state.RData"
BSFGNew = "BSFG_fixedlambda65.RData"
BSFGStar = "BSFG_fixedlambda75.RData"
target = "F_h2"
#it looks like we do not need so much latent factors. choose the first l factors.
#choose the first l latent traits and make them a vector for all iterations 
l=5 # 2:6
pos_fa = function(spn,pos.fa,n){
  pos.fa.rb = vector()
  for (i in 1:spn){
    #pos.fa.rb = c(pos.fa.rb,pos.fa[1:(l*n),i])
    pos.fa.rb = c(pos.fa.rb,pos.fa[1:l,i]) #F_h2
  }
  return(pos.fa.rb)}
# old
load(BSFGOld)
spn = dim(BSFG_state$Posterior[[target]])[2]
old.n   = dim(BSFG_state$data_matrices$Y)[1]
old.k   = nrow(BSFG_state$Posterior[[target]])
#old.pos.fa = BSFG_state$Posterior$F_a
old.pos.fh2 = BSFG_state$Posterior$F_h2

#old.pos.fa.rb = pos_fa(spn,old.pos.fa,old.n)
old.pos.fh2.rb = pos_fa(spn,old.pos.fh2,old.n)
# new
load(BSFGNew)
spn = dim(BSFG_state$Posterior[[target]])[2]
new.n   = dim(BSFG_state$data_matrices$Y)[1]
new.k   = nrow(BSFG_state$Posterior[[target]])
#new.pos.fa = BSFG_state$Posterior$F_a
new.pos.fh2 = BSFG_state$Posterior$F_h2
#new.pos.fa.rb = pos_fa(spn,new.pos.fa,new.n)
new.pos.fh2.rb = pos_fa(spn,new.pos.fh2,new.n)
# star : represent the third comparison
load(BSFGStar)
spn = dim(BSFG_state$Posterior[[target]])[2]
star.n   = dim(BSFG_state$data_matrices$Y)[1]
star.k   = nrow(BSFG_state$Posterior[[target]])
#star.pos.fa = BSFG_state$Posterior$F_a
star.pos.fh2 = BSFG_state$Posterior$F_h2
#star.pos.fa.rb = pos_fa(spn,star.pos.fa,star.n)
star.pos.fh2.rb = pos_fa(spn,star.pos.fh2,star.n)
# combine
# column "old", "new", "star"
#pos.fa.b = cbind(c(old.pos.fa.rb,new.pos.fa.rb,star.pos.fa.rb),rep(c("old","new","star"),c(old.n*l*spn,new.n*l*spn,star.n*l*spn)))
pos.fh2.b = data.frame(c(old.pos.fh2.rb,new.pos.fh2.rb,star.pos.fh2.rb),rep(c("old","new","star"),c(l*spn,l*spn,l*spn)))
# column "factor"
#pos.fa.b = data.frame(pos.fa.b,c(rep(rep(sprintf("f%d",1:l),each=old.n),spn),rep(rep(sprintf("f%d",1:l),each=new.n),spn),rep(rep(sprintf("f%d",1:l),each=star.n),spn)))
pos.fh2.b = data.frame(pos.fh2.b,c(rep(sprintf("f%d",1:l),spn),rep(sprintf("f%d",1:l),spn),rep(sprintf("f%d",1:l),spn)))
#colnames(pos.fa.b) = c("Fa","model", "factor")
colnames(pos.fh2.b) = c("Fh2","model", "factor")
#pos.fa.b=as.data.frame(pos.fa.b)
#save(pos.fa.b,file="pos_fa_bind.RData")
save(pos.fh2.b,file="pos_fh2_bind.RData")
#load("pos_fa_bind.RData")
load("pos_fh2_bind.RData")

# plot
# 1.boxplot of distribution of Fa
library(ggplot2)
ggplot(pos.fh2.b, aes(x = factor, y = Fh2, color=model,fill=model))+
  geom_boxplot()+
  ggtitle("The distribution of Latent traits")
  #geom_abline(intercept=0,slope=0,lty=2)

ggplot(pos.fh2.b,aes( x = Fh2, color=model))+
  geom_density()+
  facet_wrap(~factor,scales = "free_y")+
  ggtitle("The distribution of Latent traits")

# 2.For each trait, compare the correlation of G1 vs G2
source('~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG/BSFG_functions.R')
load(BSFGOld)
#old.traces.G = G_Traces_Comp(BSFG_state)
old.G = G_Matrix_Comp(BSFG_state)
p = BSFG_state$run_variables$p
# do this in case the value of k differ for this two models and there in no lambda in BSFG_state_fixed case
Lambda = BSFG_state$Posterior$Lambda
load(BSFGNew)
load(BSFGStar)
BSFG_state$Posterior$Lambda = Lambda
#new.traces.G = G_Traces_Comp(BSFG_state)
new.G = G_Matrix_Comp(BSFG_state)
#For all the elements in G
#pairs(t(old.traces.G)~t(new.traces.G))
#For each trait
old.G = array(unlist(old.G),dim = c(nrow(old.G[[1]]),ncol(old.G[[1]]),length(old.G)))
new.G = array(unlist(new.G),dim = c(nrow(new.G[[1]]),ncol(new.G[[1]]),length(new.G)))
old.G.b = apply(old.G,2,rbind)
new.G.b = apply(new.G,2,rbind)
save(old.G.b,file = "old_G_b.RData")
load("old_G_b.RData")
save(new.G.b,file = "new_G_b.RData")
load("new_G_b.RData")
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
  geom_abline(intercept = 0,slope = 1)+
  ggtitle(paste0(BSFG_state$traitnames[i]))
}
pdf("each_trait_G_pair.pdf")
#par(mfrow = c(3,3))
for (i in 1:p){
  each_trait_G_pair(i)
}
dev.off()
each_trait_G_pair(1)
each_trait_G_pair(2)
each_trait_G_pair(3)
each_trait_G_pair(4)
each_trait_G_pair(5)
# 3.delta z
# The distance of two matrix(in predicted phenotype value/level)
#selection gradient
beta = 0.3609
beta_v = c(0,beta,rep(0,16))
traces_traits = matrix(,p,sp_num)
for(j in 1:sp_num) {
  Lj = matrix(Posterior$Lambda[,j],p,k)
  h2j = Posterior$F_h2[,j]
  G_Lj = Lj %*%  diag(sqrt(h2j))
  Gj = G_Lj %*%  t(G_Lj) + diag(1/Posterior$E_a_prec[,j])
  Gj = sqrt(diag(VY))%*%Gj%*%t(sqrt(diag(VY)))
  traces_traits[,j] = Gj%*%beta_v
}
boxplot(t(traces_traits))
save(traces_traits,file = 'traces_traits.RData')
# For each G matrix
setwd("~/Runcie Lab/SparseFactorMixedModel_v2/MCMC/19/19/Lambda1.5_delta2shape3")
load("G_BSFG_sc.RData")
