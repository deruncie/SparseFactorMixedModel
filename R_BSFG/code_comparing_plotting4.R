source('~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG/BSFG_functions.R')
setwd("~/Runcie Lab/SparseFactorMixedModel_v2/MCMC")
MCMC = readRDS("G_MCMC.Rdata")

traitnames = rownames(MCMC[[1]][[1]])

old.G = MCMC[[5]]

p = 18
spn = 2000

new.G = MCMC[[1]]

star.G = MCMC[[8]]

four.G = MCMC[[12]]

# The distance of two matrix(in predicted phenotype value/level)
#selection gradient
beta = 0.3609
beta.v = c(0,beta,rep(0,16))
traces.traits = matrix(,p,spn)

#old.G in list form
old.traces.traits = sapply(old.G,function(x)x%*%beta.v)
#old.traits.posmean = apply(old.traces.traits,1,mean)
new.traces.traits = sapply(new.G,function(x)x%*%beta.v)
#new.traits.posmean = apply(new.traces.traits,1,mean)
star.traces.traits = sapply(star.G,function(x)x%*%beta.v)
four.traces.traits = sapply(four.G,function(x)x%*%beta.v)
#star.traits.posmean = apply(star.traces.traits,1,mean)
#boxplot(old.traces_traits[2,])
#boxplot(t(old.traces.traits))
#boxplot(new.traces_traits[2,])
#boxplot(t(new.traces.traits))
# use ggplot 
traces.traits.df = data.frame(c(c(old.traces.traits),c(new.traces.traits),c(star.traces.traits),c(four.traces.traits)))
traces.traits.df$model = rep(c("old","new","star","four"),c(p*spn,p*spn,p*spn,p*spn))

traces.traits.df$model = factor(traces.traits.df$model, level=c("old","new","star","four"))
#traces.traits.df$model = factor(traces.traits.df$model, level=c("star","new","old","four"))

traces.traits.df$traits = rep(traitnames,4*spn)
names(traces.traits.df) = c("response","model","traits")

summarySE <- function(data=NULL, measurevar, groupvars=NULL, na.rm=FALSE,
                      conf.interval=.95, .drop=TRUE) {
  library(plyr)
  
  # New version of length which can handle NA's: if na.rm==T, don't count them
  length2 <- function (x, na.rm=FALSE) {
    if (na.rm) sum(!is.na(x))
    else       length(x)
  }
  
  # This does the summary. For each group's data frame, return a vector with
  # N, mean, and sd
  datac <- ddply(data, groupvars, .drop=.drop,
                 .fun = function(xx, col) {
                   c(
                     q1 = quantile   (xx[[col]],0.025, na.rm=na.rm)- mean   (xx[[col]], na.rm=na.rm),
                     q2   = quantile     (xx[[col]],0.975, na.rm=na.rm)- mean   (xx[[col]], na.rm=na.rm),
                     mean = mean   (xx[[col]], na.rm=na.rm)
                     
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
  datac <- rename(datac, c("mean" = measurevar))
  
  #datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  #ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  #datac$ci <- datac$se * ciMult
  
  return(datac)
}

traces.traits.df.sum = summarySE(traces.traits.df,measurevar = "response",groupvars = c("model","traits"))
names(traces.traits.df.sum) = c("model","traits","q1","q2","response")


#save(traces.traits.df.sum,file = "plot6.RData")


library(ggplot2)
#ggplot(traces.traits.df,aes(x =traits,y=response,col=model,fill=model))+
# geom_boxplot()
#plot(x=rep(1:18,2),y=c(old.traits.posmean,new.traits.posmean),col=rep(c(1,2),eac#h=18),type = "o")
pd <- position_dodge(width = 0.6)
#load("plot5.RData")
ggplot(traces.traits.df.sum,aes(x = traits, y = response, color=model))+
  ylab("Response")+
  geom_errorbar(aes(ymin=q1+response, ymax=q2+response), position=pd)+
  geom_point(position=pd,size=1)


#5. PCA decomposition

old.traces.eigen = sapply(old.G, function(x)eigen(x)$values)
new.traces.eigen = sapply(new.G, function(x)eigen(x)$values)
star.traces.eigen = sapply(star.G, function(x)eigen(x)$values)
four.traces.eigen = sapply(four.G, function(x)eigen(x)$values)

traces.eigen.df = data.frame(c(c(old.traces.eigen),c(new.traces.eigen),c(star.traces.eigen),c(four.traces.eigen)))
traces.eigen.df$model = rep(c("old","new","star","four"),c(p*spn,p*spn,p*spn,p*spn))

traces.eigen.df$model = factor(traces.eigen.df$model, level=c("old","new","star","four"))
#traces.eigen.df$model = factor(traces.eigen.df$model, level=c("star","new","old","four"))

traces.eigen.df$lambda = rep(paste0("PC",1:18),4*spn)
names(traces.eigen.df) = c("response","model","PC")

traces.eigen.df.sum = summarySE(traces.eigen.df,measurevar = "response",groupvars = c("model","PC"))
names(traces.eigen.df.sum) = c("model","PC","q1","q2","response")

# reorder the PC value in the data
traces.eigen.df.sum$PC = factor(traces.eigen.df.sum$PC,level=paste0("PC",1:18))

# size dataframe
old.traces.size = apply(old.traces.eigen,2,sum)
new.traces.size = apply(new.traces.eigen,2,sum)
star.traces.size = apply(star.traces.eigen,2,sum)
four.traces.size = apply(four.traces.eigen,2,sum)
size.df = data.frame(old.traces.size,new.traces.size,star.traces.size,four.traces.size)
library(reshape2)
size.df = melt(size.df,value.name = "Size")
size.df$variable = factor(size.df$variable,levels = c("star.traces.size","new.traces.size","old.traces.size","four.traces.size"))


old.traces.pc1c = apply(old.traces.eigen,2,function(x)sum(x[2])/sum(x))
new.traces.pc1c = apply(new.traces.eigen,2,function(x)sum(x[2])/sum(x))
star.traces.pc1c = apply(star.traces.eigen,2,function(x)sum(x[2])/sum(x))
four.traces.pc1c = apply(four.traces.eigen,2,function(x)sum(x[2])/sum(x))

# PC1
old.traces.pc1c = apply(old.traces.eigen,2,function(x)x[1]/sum(x))
new.traces.pc1c = apply(new.traces.eigen,2,function(x)x[1]/sum(x))
star.traces.pc1c = apply(star.traces.eigen,2,function(x)x[1]/sum(x))
four.traces.pc1c = apply(four.traces.eigen,2,function(x)x[1]/sum(x))
pc1c.df = data.frame(old.traces.pc1c,new.traces.pc1c,star.traces.pc1c,four.traces.pc1c)
library(reshape2)
pc1c.df = melt(pc1c.df,value.name = "PC1.Contribution")
pc1c.df$variable = factor(pc1c.df$variable,levels = c("star.traces.pc1c","new.traces.pc1c","old.traces.pc1c","four.traces.pc1c"))


# size of each model 
library(dplyr)
traces.eigen.df.sum %>%
  group_by(model) %>%
  summarize(size = sum(response))

library(ggplot2)
ggplot(size.df, aes(x="",y=Size, fill=variable))+
  geom_boxplot(alpha=0.85)+
    theme(axis.title.x = element_blank())


ggplot(pc1c.df, aes(x="",y=PC1.Contribution, fill=variable))+
  geom_boxplot(alpha=0.85)+
  theme(axis.title.x = element_blank())


pd <- position_dodge(width = 0.6)

ggplot(traces.eigen.df.sum,aes(x = PC, y = response, color=model))+
   ylab("Variance")+
  geom_errorbar(aes(ymin=q1+response, ymax=q2+response), position=pd)+
  geom_point(position=pd,size=1)+
  xlab("Principle Components")


# 6. angle and length of deltaZ and first eigenvector

# distribution within population
dist_within_theta = function(traces.traits){
  theta = vector(,spn)
  for (i in 1:spn){
    pairs = randcomb(1:spn, 2)
    a = traces.traits[,pairs[1]]
    b = traces.traits[,pairs[2]]
    a.length = sqrt(sum(a * a)) 
    b.length = sqrt(sum(b * b)) 
    theta[i] = acos(sum(a*b) / (a.length * b.length))
  }
  return(theta)
}

# first eigen
old.traces.vector = sapply(old.G, function(x)eigen(x)$vector[,1])
new.traces.vector = sapply(new.G, function(x)eigen(x)$vector[,1])
star.traces.vector = sapply(star.G, function(x)eigen(x)$vector[,1])
four.traces.vector = sapply(four.G, function(x)eigen(x)$vector[,1])

old.theta = dist_within_theta(old.traces.vector)
new.theta = dist_within_theta(new.traces.vector)
star.theta = dist_within_theta(star.traces.vector)
four.theta = dist_within_theta(four.traces.vector)
# deltaZ
old.theta = dist_within_theta(old.traces.traits) #old.traces.traits is from deltaZ_posterior_mean
new.theta = dist_within_theta(new.traces.traits)
star.theta = dist_within_theta(star.traces.traits)
four.theta = dist_within_theta(four.traces.traits)

# distribution between population
dist_between_theta = function(Atraces.traits,Btraces.traits){
  theta = vector(,spn)
  for (i in 1:spn){
    pairs = randcomb(1:spn, 2)
    a = Atraces.traits[,pairs[1]]
    b = Btraces.traits[,pairs[2]]
    a.length = sqrt(sum(a * a)) 
    b.length = sqrt(sum(b * b)) 
    theta[i] = acos(sum(a*b) / (a.length * b.length))
  }
  return(theta)
}

# first eigen vector
old.new.theta = dist_between_theta(old.traces.vector,new.traces.vector)
# K-S test
#ks.test(old.theta,old.new.theta) #different
old.star.theta = dist_between_theta(old.traces.vector,star.traces.vector)
four.star.theta = dist_between_theta(four.traces.vector,star.traces.vector)
new.star.theta = dist_between_theta(new.traces.vector,star.traces.vector)
new.four.theta = dist_between_theta(new.traces.vector,four.traces.vector)
old.four.theta = dist_between_theta(old.traces.vector,four.traces.vector)

# deltaZ
old.new.theta = dist_between_theta(old.traces.traits,new.traces.traits)
# K-S test
#ks.test(old.theta,old.new.theta) #different
old.star.theta = dist_between_theta(old.traces.traits,star.traces.traits)
new.star.theta = dist_between_theta(new.traces.traits,star.traces.traits)
#--------------------------------
# plot
theta.df = data.frame(old.theta,new.theta,star.theta,four.theta, old.new.theta, old.star.theta, new.star.theta,four.star.theta,new.four.theta,old.four.theta)

library(reshape2)
theta.df = melt(theta.df, value.name = "Degree")
names(theta.df) = c("labels", "theta")
ggplot(theta.df, aes(x = theta, col=labels))+
  geom_density(aes(fill=labels),alpha=0.1)
#+scale_color_manual(breaks = c("old.theta","new.theta","star.theta","old.new.theta","old.star.theta","new.star.theta"),
#                     values=rainbow(6))

#-------------------------------
# 7.length of deltaZ

length_deltaZ = function(traces.traits){
  #posmean.deltaZ = apply(traces.traits,1,mean)
  #deltaZ.length = sqrt(sum(posmean.deltaZ * posmean.deltaZ)) 
  deltaZ.length = apply(traces.traits,2,function(x)sqrt(sum(x * x)))
  return(deltaZ.length)
}
old.deltaZ.length = length_deltaZ(old.traces.traits)
new.deltaZ.length = length_deltaZ(new.traces.traits)
star.deltaZ.length = length_deltaZ(star.traces.traits)
four.deltaZ.length = length_deltaZ(four.traces.traits)

deltaZ.length.df = data.frame(old.deltaZ.length,new.deltaZ.length,star.deltaZ.length,four.deltaZ.length)
library(reshape2)
deltaZ.length.df = melt(deltaZ.length.df,value.name = "Length")
deltaZ.length.df$variable = factor(deltaZ.length.df$variable,levels = c("star.deltaZ.length","new.deltaZ.length","old.deltaZ.length","four.deltaZ.length"))



ggplot(deltaZ.length.df, aes(x="",y=Length, fill=variable))+
  geom_boxplot(alpha=0.85)+
   theme(axis.title.x = element_blank())
 

