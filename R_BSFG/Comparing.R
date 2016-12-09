# First method
# Fixed lambda
# all the BSFG G matrix
load("~/../Desktop/G_BSFG_sc.RData")
# 5 G matrix
BSFG.sc5 = G_BSFG_sc[5]

model_path = "~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG"
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
# 5 BSFG_state
setwd("~/Runcie Lab/SparseFactorMixedModel_v2/MCMC/5/5/Lambda1.5_delta2shape3")
load("BSFG_state.RData")
priors = BSFG_state$priors

YNew="setup6.RData"
YOld="setup5.RData"
BSFG_state = fast_BSFG_sampler_init_fixedlambda(priors,run_parameters,YNew,YOld)
BSFG_state = fast_BSFG_sampler_fixedlambda(BSFG_state,YNew,YOld)

#load("BSFG_fixedlambda65.RData")
BSFGOld = "BSFG_state.RData"
BSFGNew = "BSFG_fixedlambda65.RData"
BSFG3 = 
# Comparisons:

  load(BSFGOld)
  spn = dim(BSFG_state$Posterior[[target]])[2]
  old.n   = dim(BSFG_state$data_matrices$Y)[1]
  old.k   = nrow(BSFG_state$Posterior[[target]])
  #pos = BSFG_state$Posterior[[target]][,spn]
  #pos = BSFG_state$Posterior[[target]]
  old.pos.fa = BSFG_state$Posterior$F_a
  old.pos.fa.rb = vector()
  l=5
  for (i in 1:spn){
  old.pos.fa.rb = c(old.pos.fa.rb,old.pos.fa[1:(l*old.n),i])
  }
  #combine posterior of F_a into a big matrix

  #old.pos.fa.matrix = matrix(0,nr=0,nc=l)
  #for (i in 1:spn){
  #  old.fa.matrix = matrix(old.pos.fa[,i],nr=old.n,nc=k)
  #  old.pos.fa.matrix = rbind(old.pos.fa.matrix,old.fa.matrix[,1:l])
  #}
  #boxplot(old.pos.fa.matrix)
  # combine a column which specify the latent factor
  
  
  load(BSFGNew)
  spn = dim(BSFG_state$Posterior[[target]])[2]
  new.n   = dim(BSFG_state$data_matrices$Y)[1]
  new.k   = nrow(BSFG_state$Posterior[[target]])
  #pos = BSFG_state$Posterior[[target]][,spn]
  #pos = BSFG_state$Posterior[[target]]
  new.pos.fa = BSFG_state$Posterior$F_a
  #combine posterior of F_a into a big matrix
  l=5
  new.pos.fa.rb = vector()
  for (i in 1:spn){
    new.pos.fa.rb = c(new.pos.fa.rb,new.pos.fa[1:(l*new.n),i])
  }
  #new.pos.fa.matrix = matrix(0,nr=0,nc=l)
  #for (i in 1:spn){
  #  new.fa.matrix = matrix(new.pos.fa[,i],nr=new.n,nc=k)
  #  new.pos.fa.matrix = rbind(new.pos.fa.matrix,new.fa.matrix[,1:l])
  #}
  #boxplot(new.pos.fa.matrix)
  
 
  
#combine old pos.fa.matrix and new pos.fa.matrix into one dataframe 
# use ggplot to plot
#it looks like we do not need so much latent factors. choose the first8 factors.
  #pos.fa.matrix = cbind(rbind(old.pos.fa.matrix,new.pos.fa.matrix),rep(c("old","new"),c(old.n*spn,new.n*spn)))
  # column "old", "new"
  pos.fa.b = cbind(c(old.pos.fa.rb,new.pos.fa.rb),rep(c("old","new"),c(old.n*l*spn,new.n*l*spn)))
  # column "factor"
  pos.fa.b = cbind(pos.fa.b,c(rep(rep(sprintf("f%d",1:l),each=old.n),spn),rep(rep(sprintf("f%d",1:l),each=new.n),spn)))
  colnames(pos.fa.b) = c("Fa","model", "factor")
  pos.fa.b=as.data.frame(pos.fa.b)
  #pos.fa.sample = sample(1:nrow(pos.fa.b),200*(old.n+new.n))
  #pos.fa.sample = pos.fa.b[pos.fa.sample,]
  #plot boxplot of distribution of Fa
  library(ggplot2)
  ggplot(pos.fa.b,aes(x = factor, y = as.numeric(Fa), color=model,fill=model))+
    geom_boxplot()+
    geom_abline(a=0,b=0,lty=2)
  # For each trait, compare the correlation of G1 vs G2
  
  # 
  
  
  if (target!="F_h2"){
    pos = matrix(pos,nr=n)
    k   = nrow(BSFG_state$Posterior[[target]])/n
  }  
  #load data from new population
  load(BSFGNew)
  n   = dim(BSFG_state$data_matrices$Y)[1]
  star = BSFG_state$Posterior[[target]][,spn]    
  if (target!="F_h2"){
    star = matrix(star,nr=n)
    pdf(sprintf("comparing_%s_densityplot.pdf",target))
    for(i in 1:k){
      plot(density(pos[,i]),main = sprintf("%s %d",target,i),col = "blue",type = "l",xlab = "#obs")
      # plot(density(F_a_pos[,i]),main = sprintf("%d",i),col = "blue",type = "l",ylim = c(min(F_a_pos)-5,max(F_a_pos)+5),xlab = "#obs")
      lines(density(star[,i]),col="red",type = "l")
      abline
      legend("topright",legend = c("original","new"),col = c("blue","red"),text.col = c("blue","red"),bty = "n",pch = 1)
    }
    dev.off()
  }else{
    pdf(sprintf("comparing_%s_densityplot.pdf",target))
    plot(density(pos),main = sprintf("%s",target),col = "blue",type = "l",xlab = "#obs")
    # plot(density(F_a_pos[,i]),main = sprintf("%d",i),col = "blue",type = "l",ylim = c(min(F_a_pos)-5,max(F_a_pos)+5),xlab = "#obs")
    lines(density(star),col="red",type = "l")
    abline
    legend("topright",legend = c("original","new"),col = c("blue","red"),text.col = c("blue","red"),bty = "n",pch = 1)
    dev.off()
  }
}

ComparingGMatrix_plot("F_h2",BSFGNew,BSFGOld)



