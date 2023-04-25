source('~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG/BSFG_functions.R')
setwd("~/Runcie Lab/SparseFactorMixedModel_v2/MCMC")
MCMC = readRDS("G_MCMC.Rdata")
#setwd("~/Runcie Lab/SparseFactorMixedModel_v2/MCMC/5/5/Lambda1.5_delta2shape3")
#BSFGOld = "BSFG_state.RData"
#BSFGNew = "BSFG_state14.RData"
#BSFGStar = "BSFG_state17.RData"
#target = "F_h2"
#load(BSFGOld)
#traitnames = BSFG_state$traitnames
traitnames = rownames(MCMC[[1]][[1]])

#old.G = G_Matrix_Comp(BSFG_state)
old.G = MCMC[[4]]
#p = BSFG_state$run_variables$p
#spn = dim(BSFG_state$Posterior[[target]])[2]
p = 18
spn = 2000
#load(BSFGNew)
#new.G = G_Matrix_Comp(BSFG_state)
new.G = MCMC[[3]]

#load(BSFGStar)
#star.G = G_Matrix_Comp(BSFG_state)
star.G = MCMC[[2]]


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
#star.traits.posmean = apply(star.traces.traits,1,mean)
#boxplot(old.traces_traits[2,])
#boxplot(t(old.traces.traits))
#boxplot(new.traces_traits[2,])
#boxplot(t(new.traces.traits))
# use ggplot 
traces.traits.df = data.frame(c(c(old.traces.traits),c(new.traces.traits),c(star.traces.traits)))
traces.traits.df$model = rep(c("old","new","star"),c(p*spn,p*spn,p*spn))

traces.traits.df$model = factor(traces.traits.df$model, level=c("old","new","star"))
traces.traits.df$model = factor(traces.traits.df$model, level=c("star","new","old"))

traces.traits.df$traits = rep(traitnames,3*spn)
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


save(traces.traits.df.sum,file = "plot6.RData")


library(ggplot2)
#ggplot(traces.traits.df,aes(x =traits,y=response,col=model,fill=model))+
# geom_boxplot()
#plot(x=rep(1:18,2),y=c(old.traits.posmean,new.traits.posmean),col=rep(c(1,2),eac#h=18),type = "o")
pd <- position_dodge(width = 0.6)
load("plot5.RData")
plot5 <- ggplot(traces.traits.df.sum,aes(x = traits, y = response, color=model))+
  scale_color_manual(breaks =rev(c("star","new","old")),
                     values=rev(c( "blue","red2","grey0")),labels=rev(c("Low Line","High Line","Control Line")))+
  ylab("Response")+
  geom_errorbar(aes(ymin=q1+response, ymax=q2+response), position=pd)+
  geom_point(position=pd,size=1)+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=15),legend.text=element_text(size=12),
        legend.title=element_text(size=12))+
  labs(fill = "G (F3)",color = "G (F3)")


library(ggplot2)
#ggplot(traces.traits.df,aes(x =traits,y=response,col=model,fill=model))+
# geom_boxplot()
#plot(x=rep(1:18,2),y=c(old.traits.posmean,new.traits.posmean),col=rep(c(1,2),eac#h=18),type = "o")
#pd <- position_dodge(width = 0.3)
load("plot6.RData")
plot6<- ggplot(traces.traits.df.sum,aes(x = traits, y = response, color=model))+
  scale_color_manual(breaks = c("star","new","old"),
                     values=c("grey60","grey30","grey0"),labels=c("F1","F2","F3"))+
  ylab("Response")+
  geom_errorbar(aes(ymin=q1+response, ymax=q2+response), position=pd)+
  geom_point(position=pd,size=1)+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=15),legend.text=element_text(size=12),
        legend.title=element_text(size=12))+
  labs(fill = "G (Control)",color = "G (Control)")


plots1=list()
plots1=lapply(c(1,3),function(i) get(paste0("plot",i)))
plots2=list()
plots2=lapply(c(2,4),function(i) get(paste0("plot",i)))
plots3=list()
plots3=lapply(c(5,6),function(i) get(paste0("plot",i)))
plots4=list(g11,g22) 
library(gridExtra)
ggsave("BSFG_comparison_plots1.pdf", marrangeGrob(grobs = plots1, nrow=2, ncol=1,top = NULL),width = 297, height = 210, units = "mm")
ggsave("BSFG_comparison_plots2.pdf", marrangeGrob(grobs = plots2, nrow=2, ncol=1,top = NULL),width = 297, height = 210, units = "mm")
ggsave("BSFG_comparison_plots3.pdf", marrangeGrob(grobs = plots3, nrow=2, ncol=1,top = NULL),width = 297, height = 210, units = "mm")
ggsave("EigenvaluesComparison.pdf", marrangeGrob(grobs = plots4, nrow=2, ncol=1,top = NULL),width = 297, height = 210, units = "mm")


#5. PCA decomposition

source('~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG/BSFG_functions.R')
#load(BSFGOld)
#old.G = G_Matrix_Comp(BSFG_state)
old.G = MCMC[[4]]
new.G = MCMC[[3]]
star.G = MCMC[[2]]
#p = BSFG_state$run_variables$p
#traitnames = BSFG_state$traitnames
#load("BSFG_state6.RData")
#new.G = G_Matrix_Comp(BSFG_state)
#load("BSFG_state7.RData")
#star.G = G_Matrix_Comp(BSFG_state)

old.traces.eigen = sapply(old.G, function(x)eigen(x)$values)
new.traces.eigen = sapply(new.G, function(x)eigen(x)$values)
star.traces.eigen = sapply(star.G, function(x)eigen(x)$values)

traces.eigen.df = data.frame(c(c(old.traces.eigen),c(new.traces.eigen),c(star.traces.eigen)))
traces.eigen.df$model = rep(c("old","new","star"),c(p*spn,p*spn,p*spn))

traces.eigen.df$model = factor(traces.eigen.df$model, level=c("old","new","star"))
traces.eigen.df$model = factor(traces.eigen.df$model, level=c("star","new","old"))

traces.eigen.df$lambda = rep(paste0("PC",1:18),3*spn)
names(traces.eigen.df) = c("response","model","PC")

traces.eigen.df.sum = summarySE(traces.eigen.df,measurevar = "response",groupvars = c("model","PC"))
names(traces.eigen.df.sum) = c("model","PC","q1","q2","response")

# reorder the PC value in the data
traces.eigen.df.sum$PC = factor(traces.eigen.df.sum$PC,level=paste0("PC",1:18))

# size dataframe
old.traces.size = apply(old.traces.eigen,2,sum)
new.traces.size = apply(new.traces.eigen,2,sum)
star.traces.size = apply(star.traces.eigen,2,sum)
size.df = data.frame(old.traces.size,new.traces.size,star.traces.size)
library(reshape2)
size.df = melt(size.df,value.name = "Size")
size.df$variable = factor(size.df$variable,levels = c("star.traces.size","new.traces.size","old.traces.size"))

# PC1
old.traces.pc1c = apply(old.traces.eigen,2,function(x)x[1]/sum(x))
new.traces.pc1c = apply(new.traces.eigen,2,function(x)x[1]/sum(x))
star.traces.pc1c = apply(star.traces.eigen,2,function(x)x[1]/sum(x))
pc1c.df = data.frame(old.traces.pc1c,new.traces.pc1c,star.traces.pc1c)
library(reshape2)
pc1c.df = melt(pc1c.df,value.name = "PC1.Contribution")
pc1c.df$variable = factor(pc1c.df$variable,levels = c("star.traces.pc1c","new.traces.pc1c","old.traces.pc1c"))



save(traces.eigen.df.sum , file = "eigen567.RData")
save(traces.eigen.df.sum, file = "eigen51417.RData")
save(size.df, file = "size567.RData")
save(size.df, file = "size51417.RData")
save(pc1c.df, file = "pc1c567.RData")
save(pc1c.df, file = "pc1c51417.RData")

load("eigen567.RData")
load("size567.RData")
load("pc1c567.RData")
# size of each model 
library(dplyr)
traces.eigen.df.sum %>%
  group_by(model) %>%
  summarize(size = sum(response))

library(ggplot2)
size567.boxplot <- ggplot(size.df, aes(x="",y=Size, fill=variable))+
  geom_boxplot(alpha=0.85)+
  scale_fill_manual(breaks =rev(c("star.traces.size","new.traces.size","old.traces.size")),
                     values=rev(c( "blue","red2","grey0")),labels=rev(c("Low Line","High Line","Control Line")))+
  theme(axis.title.x = element_blank())+
  labs(fill = "G (F3)")

pc1c567.boxplot <- ggplot(pc1c.df, aes(x="",y=PC1.Contribution, fill=variable))+
  geom_boxplot(alpha=0.85)+
  scale_fill_manual(breaks =rev(c("star.traces.pc1c","new.traces.pc1c","old.traces.pc1c")),
                    values=rev(c( "blue","red2","grey0")),labels=rev(c("Low Line","High Line","Control Line")))+
  theme(axis.title.x = element_blank())+
  labs(fill = "G (F3)")

pd <- position_dodge(width = 0.6)

g11<-ggplot(traces.eigen.df.sum,aes(x = PC, y = response, color=model))+
  scale_color_manual(breaks =rev(c("star","new","old")),
                     values=rev(c( "blue","red2","grey0")),labels=rev(c("Low Line","High Line","Control Line")))+
  ylab("Variance")+
  geom_errorbar(aes(ymin=q1+response, ymax=q2+response), position=pd)+
  geom_point(position=pd,size=1)+
  xlab("Principle Components")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=15),legend.text=element_text(size=12),
        legend.title=element_text(size=12))+
  labs(color = "G (F3)")
#-------------------------------
# combine two boxplot
library(gridExtra)
library("cowplot")
plot_grid(size567.boxplot, size51417.boxplot, ncol = 2, nrow = 1)
plot_grid(pc1c567.boxplot, pc1c51417.boxplot, ncol = 2, nrow = 1)
#-------------------------------
load("eigen51417.RData")
load("size51417.RData")
load("pc1c51417.RData")
# size of each model
library(dplyr)
traces.eigen.df.sum %>%
  group_by(model) %>%
  summarize(size = sum(response))

# boxplot of size of each model
library(ggplot2)
size51417.boxplot <- ggplot(size.df, aes(x="",y=Size, fill=variable))+
  geom_boxplot()+
  scale_fill_manual(breaks = c("star.traces.size","new.traces.size","old.traces.size"),
                     values=c("grey60","grey30","grey0"),labels=c("F1","F2","F3"))+
  theme(axis.title.x = element_blank())+
  labs(fill = "G (Control)")

pc1c51417.boxplot = ggplot(pc1c.df, aes(x="",y=PC1.Contribution, fill=variable))+
  geom_boxplot()+
  scale_fill_manual(breaks = c("star.traces.pc1c","new.traces.pc1c","old.traces.pc1c"),
                    values=c("grey60","grey30","grey0"),labels=c("F1","F2","F3"))+
  theme(axis.title.x = element_blank())+
  labs(fill = "G (Control)")

g22<-ggplot(traces.eigen.df.sum,aes(x = reorder(PC,PC), y = response, color=model))+
  scale_color_manual(breaks = c("star","new","old"),
                     values=c("grey60","grey30","grey0"),labels=c("F1","F2","F3"))+
  ylab("Response")+
  geom_errorbar(aes(ymin=q1+response, ymax=q2+response), position=pd)+
  geom_point(position=pd,size=1)+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=15),legend.text=element_text(size=12),
        legend.title=element_text(size=12),legend.margin=unit(0.455,"cm"))+
  labs(color = "G (Control)")+
  xlab("Principle Components")+
  ylab("Variance")
# statistical test for the difference
# a.ANOVA test
class(traces.eigen.df)
anova(lm(response~model, data = traces.eigen.df))
anova(lm(response~model:PC, data = traces.eigen.df)) #which one to use?
anova(lm(response~model*PC, data = traces.eigen.df))

# b.Bartlett test

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

old.theta = dist_within_theta(old.traces.vector)
new.theta = dist_within_theta(new.traces.vector)
star.theta = dist_within_theta(star.traces.vector)
# deltaZ
old.theta = dist_within_theta(old.traces.traits) #old.traces.traits is from deltaZ_posterior_mean
new.theta = dist_within_theta(new.traces.traits)
star.theta = dist_within_theta(star.traces.traits)


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
new.star.theta = dist_between_theta(new.traces.vector,star.traces.vector)


# deltaZ
old.new.theta = dist_between_theta(old.traces.traits,new.traces.traits)
# K-S test
#ks.test(old.theta,old.new.theta) #different
old.star.theta = dist_between_theta(old.traces.traits,star.traces.traits)
new.star.theta = dist_between_theta(new.traces.traits,star.traces.traits)
#--------------------------------
# plot
theta.df = data.frame(old.theta,new.theta,star.theta, old.new.theta, old.star.theta, new.star.theta)

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

deltaZ.length.df = data.frame(old.deltaZ.length,new.deltaZ.length,star.deltaZ.length)
library(reshape2)
deltaZ.length.df = melt(deltaZ.length.df,value.name = "Length")
deltaZ.length.df$variable = factor(deltaZ.length.df$variable,levels = c("star.deltaZ.length","new.deltaZ.length","old.deltaZ.length"))
save(deltaZ.length.df, file="traitslength567.RData")
save(deltaZ.length.df, file="traitslength51417.RData")

load("traitslength567.RData")
traits.length567.boxplot <- ggplot(deltaZ.length.df, aes(x="",y=Length, fill=variable))+
  geom_boxplot(alpha=0.85)+
  scale_fill_manual(breaks =rev(c("star.deltaZ.length","new.deltaZ.length","old.deltaZ.length")),
                    values=rev(c( "blue","red2","grey0")),labels=rev(c("Low Line","High Line","Control Line")))+
  theme(axis.title.x = element_blank())+
  labs(fill = "G (F3)")

load("traitslength51417.RData")
traits.length51417.boxplot <- ggplot(deltaZ.length.df, aes(x="",y=Length, fill=variable))+
  geom_boxplot()+
  scale_fill_manual(breaks = c("star.deltaZ.length","new.deltaZ.length","old.deltaZ.length"),
                    values=c("grey60","grey30","grey0"),labels=c("F1","F2","F3"))+
  theme(axis.title.x = element_blank())+
  labs(fill = "G (Control)")

library(gridExtra)
library("cowplot")
plot_grid(traits.length567.boxplot, traits.length51417.boxplot, ncol = 2, nrow = 1)
