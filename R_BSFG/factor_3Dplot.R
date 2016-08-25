setwd('~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG/Sim_1/R_rep_test')
setwd('~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG/Sim_1/R_rep_fixedlambda')

#function for plotting 3D plot
rgl_add_axes <- function(x, y, z, axis.col = "grey",
                         xlab = "", ylab="", zlab="", show.plane = TRUE, 
                         show.bbox = FALSE, bbox.col = c("#333377","black"))
{ 
  
  lim <- function(x){c(-max(abs(x)), max(abs(x))) * 1.1}
  # Add axes
  xlim <- lim(x); ylim <- lim(y); zlim <- lim(z)
  rgl.lines(xlim, c(0, 0), c(0, 0), color = axis.col)
  rgl.lines(c(0, 0), ylim, c(0, 0), color = axis.col)
  rgl.lines(c(0, 0), c(0, 0), zlim, color = axis.col)
  
  # Add a point at the end of each axes to specify the direction
  axes <- rbind(c(xlim[2], 0, 0), c(0, ylim[2], 0), 
                c(0, 0, zlim[2]))
  rgl.points(axes, color = axis.col, size = 3)
  
  # Add axis labels
  rgl.texts(axes, text = c(xlab, ylab, zlab), color = axis.col,
            adj = c(0.5, -0.8), size = 2)
  
  # Add plane
  if(show.plane) 
    xlim <- xlim/1.1; zlim <- zlim /1.1
  rgl.quads( x = rep(xlim, each = 2), y = c(0, 0, 0, 0),
             z = c(zlim[1], zlim[2], zlim[2], zlim[1]))
  
  # Add bounding box decoration
  if(show.bbox){
    rgl.bbox(color=c(bbox.col[1],bbox.col[2]), alpha = 0.5, 
             emission=bbox.col[1], specular=bbox.col[1], shininess=5, 
             xlen = 3, ylen = 3, zlen = 3) 
  }
}


load('BSFG_state.RData')
n_old = BSFG_state$run_variables$n
F_a_old = BSFG_state$Posterior$F_a
F_a_old = matrix(F_a_old[,ncol(F_a_old)],nr=n_old)
#Y_mean = apply(na.omit(BSFG_state$data_matrices$Y),2,mean)
#Y_sd = apply(na.omit(BSFG_state$data_matrices$Y),2,sd)
#F_a_old = apply(F_a_old,2,function(x) x*Y_sd+Y_mean)

load('BSFG_fixedlambda2.RData')
F_a_new2 = BSFG_state$Posterior$F_a
n_new2 = BSFG_state$run_variables$n
F_a_new2 = matrix(F_a_new2[,ncol(F_a_new2)],nr=n_new2)
#F_a_new = F_a_new*Y_sd+Y_mean

load('BSFG_fixedlambda3.RData')
F_a_new3 = BSFG_state$Posterior$F_a
n_new3 = BSFG_state$run_variables$n
F_a_new3 = matrix(F_a_new3[,ncol(F_a_new3)],nr=n_new3)





#=======================================================
factor_3Dplot = function(x1,x2,x3){
  
  if(exists('n_new3')){
    df = data.frame(rbind(cbind(F_a_old[,x1],F_a_old[,x2],F_a_old[,x3]),
                          cbind(F_a_new2[,x1],F_a_new2[,x2],F_a_new2[,x3]),
                          cbind(F_a_new3[,x1],F_a_new3[,x2],F_a_new3[,x3])))
    df$factor = factor(c(rep("black",n_old),rep("red",n_new2),
                         rep("blue",n_new3)))
  }else{
    df = data.frame(rbind(cbind(F_a_old[,x1],F_a_old[,x2],F_a_old[,x3]),
                          cbind(F_a_new2[,x1],F_a_new2[,x2],F_a_new2[,x3])))
    df$factor = factor(c(rep("green",n_old),rep("yellow",n_new2)))
  }
  
  require(rgl)
  plot3d(df$X1,df$X2,df$X3,col=df$factor,type = "p",  xlab="x1",ylab='x2',zlab='x3')
  title3d(main = "3D plot of three factors")
  rgl.bbox(color=c("#333377","black"), emission="#333377",
           specular="#3333FF", shininess=5, alpha=0.8 )
  lim <- function(x){c(-max(abs(x)), max(abs(x))) * 1.1}
  # Add axes
  axis.col = "grey"
  xlim <- lim(df$X1); ylim <- lim(df$X2); zlim <- lim(df$X3)
  rgl.lines(xlim, c(0, 0), c(0, 0), color = axis.col)
  rgl.lines(c(0, 0), ylim, c(0, 0), color = axis.col)
  rgl.lines(c(0, 0), c(0, 0), zlim, color = axis.col)
  
  aspect3d(1,1,1)
 #create a movie 3d
  movie3d(spin3d(axis = c(0, 0, 1)), duration = 10,dir = getwd())
}

#-------------------------------------------
#plot~~~~~
#latent factor we choose
x1=1
x2=2
x3=10
factor_3Dplot(x1,x2,x3)

#How much contribution for each of these factors
load('BSFG_state.RData')
p = BSFG_state$run_variables$p

Lambda = BSFG_state$Posterior$Lambda
Lambda = matrix(Lambda[,ncol(Lambda)],nr=p)

#the size of x1, x2, x3 columns


