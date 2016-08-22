setwd('~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG/Sim_1/R_rep_test')
setwd('~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG/Sim_1/R_rep_fixedlambda')

#function for plotting 3D plot
factor_3Dplot = function(x1,x2,x3){
  load('BSFG_state.RData')

  n_old = BSFG_state$run_variables$n
  F_a_old = BSFG_state$Posterior$F_a
  F_a_old = matrix(F_a_old[,ncol(F_a_old)],nr=n_old)
  #Y_mean = apply(na.omit(BSFG_state$data_matrices$Y),2,mean)
  #Y_sd = apply(na.omit(BSFG_state$data_matrices$Y),2,sd)
  #F_a_old = apply(F_a_old,2,function(x) x*Y_sd+Y_mean)
  
  load('BSFG_fixedlambda.RData')
  F_a_new = BSFG_state$Posterior$F_a
  n_new = BSFG_state$run_variables$n
  F_a_new = matrix(F_a_new[,ncol(F_a_new)],nr=n_new)
  #F_a_new = F_a_new*Y_sd+Y_mean
  df = data.frame(rbind(cbind(F_a_old[,x1],F_a_old[,x2],F_a_old[,x3]),
                        cbind(F_a_new[,x1],F_a_new[,x2],F_a_new[,x3])))
  
  
  df$factor = factor(c(rep("green",n_old),rep("yellow",n_new)))
  
  require(rgl)
  plot3d(df$X1,df$X2,df$X3,col=df$factor,xlab="x1",ylab='x2',zlab='x3')
  title3d(main = "3D plot of three factors")
  rgl.bbox(color=c("#333377","black"), emission="#333377",
           specular="#3333FF", shininess=5, alpha=0.8 )
#aspect3d(1,1,1)
}

#-------------------------------------------
#plot~~~~~
#latent factor we choose
x1=1
x2=2
x3=5
factor_3Dplot(x1,x2,x3)

#How much contribution for each of these factors
load('BSFG_state.RData')
p = BSFG_state$run_variables$p

Lambda = BSFG_state$Posterior$Lambda
Lambda = matrix(Lambda[,ncol(Lambda)],nr=p)

#the size of x1, x2, x3 columns
