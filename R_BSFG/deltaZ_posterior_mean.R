#source('~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG/BSFG_functions.R')
#setwd("~/Runcie Lab/SparseFactorMixedModel_v2/MCMC/5/5/Lambda1.5_delta2shape3")
#BSFGOld = "BSFG_state.RData"
#BSFGNew = "BSFG_fixedlambda65.RData"
#BSFGStar = "BSFG_fixedlambda75.RData"
#BSFGNew = "BSFG_state6.RData"
#BSFGStar = "BSFG_state7.RData"
#target = "F_h2"
#load(BSFGOld)
#traitnames = BSFG_state$traitnames
#old.traces.G = G_Traces_Comp(BSFG_state)
#old.G = G_Matrix_Comp(BSFG_state)
#p = BSFG_state$run_variables$p
#spn = dim(BSFG_state$Posterior[[target]])[2]
# posterior mean
# do this in case the value of k differ for this two models and there in no lambda in BSFG_state_fixed case
#Lambda = BSFG_state$Posterior$Lambda
#load(BSFGNew)
#BSFG_state$Posterior$Lambda = Lambda
#new.G = G_Matrix_Comp(BSFG_state)

#load(BSFGStar)
#BSFG_state$Posterior$Lambda = Lambda
#star.G = G_Matrix_Comp(BSFG_state)
#new.traces.G = G_Traces_Comp(BSFG_state)


# The distance of two matrix(in predicted phenotype value/level)
#selection gradient
#beta = 0.3609
#beta.v = c(0,beta,rep(0,16))
#traces.traits = matrix(,p,spn)

#old.G in list form
#old.traces.traits = sapply(old.G,function(x)x%*%beta.v)
#old.traits.posmean = apply(old.traces.traits,1,mean)
#new.traces.traits = sapply(new.G,function(x)x%*%beta.v)
#new.traits.posmean = apply(new.traces.traits,1,mean)
#star.traces.traits = sapply(star.G,function(x)x%*%beta.v)
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
