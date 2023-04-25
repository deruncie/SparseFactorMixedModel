setwd("~/Runcie Lab/SparseFactorMixedModel_v2/MCMC/5/5/Lambda1.5_delta2shape3")
BSFGOld = "BSFG_state1.RData"
BSFGNew = "BSFG_state2.RData"
BSFGStar = "BSFG_state3.RData"
target = "F_h2"

#function change list into a array
L_to_A = function(l){
  l.to.a = array(unlist(l),dim = c(nrow(l[[1]]),ncol(l[[1]]),length(l)))
  return(l.to.a)
}

#To interpret lambda/latent factors
#Lambda.array = L_to_A(BSFG_state$Posterior$Lambda)
#take posterior mean
Lambda.pos = apply(Lambda,1,mean)
Lambda.pos = matrix(Lambda.pos,nr=p)

Lambda.pos[,2]=(abs(Lambda.pos[,1])+abs(Lambda.pos[,2]))/2
Lambda.pos=Lambda.pos[,-1]

rownames(Lambda.pos) = BSFG_state$traitnames
pdf("interpret_lambda1.pdf")
par(mfrow=c(3,3))
for (i in 1:17){
  names(Lambda.pos[,i]) = rownames(Lambda.pos)
plot(Lambda.pos[,i], xlab = paste0(i),col=1:18,)
  
abline(0,0)
}
dev.off()
save(Lambda.pos,file="LambdaMatrix.RData")
#it looks like we do not need so much latent factors. choose the first l factors.
#choose the first l latent traits and make them a vector for all iterations 
l=8 # 2:6
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
old.k   = nrow(BSFG_state$Posterior[[target]])-1
#old.pos.fa = BSFG_state$Posterior$F_a
old.pos.fh2 = BSFG_state$Posterior$F_h2
old.pos.fh2[2,] = (old.pos.fh2[1,]+old.pos.fh2[2,])/2
old.pos.fh2 = old.pos.fh2[-1,]
#old.pos.fa.rb = pos_fa(spn,old.pos.fa,old.n)
old.pos.fh2.rb = pos_fa(spn,old.pos.fh2,old.n)
# new
load(BSFGNew)
spn = dim(BSFG_state$Posterior[[target]])[2]
new.n   = dim(BSFG_state$data_matrices$Y)[1]
new.k   = nrow(BSFG_state$Posterior[[target]])-1
#new.pos.fa = BSFG_state$Posterior$F_a
new.pos.fh2 = BSFG_state$Posterior$F_h2
new.pos.fh2[2,] = (new.pos.fh2[1,]+new.pos.fh2[2,])/2
new.pos.fh2 = new.pos.fh2[-1,]
#new.pos.fa.rb = pos_fa(spn,new.pos.fa,new.n)
new.pos.fh2.rb = pos_fa(spn,new.pos.fh2,new.n)
# star : represent the third comparison
load(BSFGStar)
spn = dim(BSFG_state$Posterior[[target]])[2]
star.n   = dim(BSFG_state$data_matrices$Y)[1]
star.k   = nrow(BSFG_state$Posterior[[target]])-1
#star.pos.fa = BSFG_state$Posterior$F_a
star.pos.fh2 = BSFG_state$Posterior$F_h2
star.pos.fh2[2,] = (star.pos.fh2[1,]+star.pos.fh2[2,])/2
star.pos.fh2 = star.pos.fh2[-1,]
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
#save(pos.fh2.b,file="pos_fh2_bind.RData")
#load("pos_fa_bind.RData")
#load("pos_fh2_bind.RData")
pos.fh2.b$model = factor(pos.fh2.b$model, level=c("old","new","star"))
pos.fh2.b$model = factor(pos.fh2.b$model, level=c("star","new","old"))
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
                     q0 = min(xx[[col]], na.rm=na.rm),
                     q1 = quantile   (xx[[col]],0.025, na.rm=na.rm),
                     q2 = median(xx[[col]], na.rm=na.rm),
                      q3   = quantile     (xx[[col]],0.975, na.rm=na.rm),
                     q4 = max  (xx[[col]], na.rm=na.rm)
                     
                   )
                 },
                 measurevar
  )
  
  # Rename the "mean" column    
 # datac <- rename(datac, c("mmm" = measurevar))
  
  #datac$se <- datac$sd / sqrt(datac$N)  # Calculate standard error of the mean
  
  # Confidence interval multiplier for standard error
  # Calculate t-statistic for confidence interval: 
  # e.g., if conf.interval is .95, use .975 (above/below), and use df=N-1
  #ciMult <- qt(conf.interval/2 + .5, datac$N-1)
  #datac$ci <- datac$se * ciMult
  
  return(datac)
}

pos.fh2.b.sum = summarySE(pos.fh2.b,measurevar = "Fh2",groupvars = c("model","factor"))

names(pos.fh2.b.sum) = c("model","factor","q0","q1" ,"q2","q3","q4" )

save(pos.fh2.b,file = "plot3.RData")

save(pos.fh2.b.sum,file="plot4.RData")

# plot
# 1.boxplot of distribution of Fh2
library(ggplot2)
load("plot1.RData")
plot1 <- ggplot(pos.fh2.b, aes(x = factor, y = Fh2, fill=model,color=model))+
  geom_boxplot(alpha = 0.6,position="dodge")+
  scale_color_manual(breaks =rev(c("star","new","old")),
                     values=rev(c( "blue","red2","grey0")),labels=rev(c("Low Line","High Line","Control Line")))+
  scale_fill_manual(breaks =rev(c("star","new","old")),
                    values=rev(c( "blue","red2","grey0")),labels=rev(c("Low Line","High Line","Control Line")))+
  ylab("Variance")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=15),legend.text=element_text(size=12),
        legend.title=element_text(size=12))+
  labs(fill = "G (F3)",color = "G (F3)")


# 95% boxplot
load("plot2.RData")
plot2<- ggplot(pos.fh2.b.sum, aes(x = factor, color=model,fill=model))+
  geom_boxplot(alpha = 0.5,aes(ymin=q0,lower=q1,middle=q2,upper=q3,ymax=q4),stat = "identity")+
  scale_color_manual(breaks =rev(c("star","new","old")),
                     values=rev(c( "blue","red2","grey0")),labels=rev(c("Low Line","High Line","Control Line")))+
  scale_fill_manual(breaks =rev(c("star","new","old")),
                    values=rev(c( "blue","red2","grey0")),labels=rev(c("Low Line","High Line","Control Line")))+
  #geom_boxplot(alpha = 0.5)+
  ylab("Variance")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=15),legend.text=element_text(size=12),
        legend.title=element_text(size=12))+
  labs(fill = "G (F3)",color = "G (F3)")

  #geom_abline(intercept=0,slope=0,lty=2)
load("plot3.RData")
plot3<- ggplot(pos.fh2.b, aes(x = factor,y = Fh2, color=model,fill=model))+
  geom_boxplot(alpha = 0.5)+
  #geom_boxplot(alpha = 0.5,aes(ymin=q0,lower=q1,middle=q2,upper=q3,ymax=q4),stat = "identity")+
  scale_color_manual(breaks = c("star","new","old"),
                     values=c("grey60","grey30","grey0"),labels=c("F1","F2","F3"))+
  scale_fill_manual(breaks = c("star","new","old"),
                    values=c("grey60","grey30","grey0"),labels=c("F1","F2","F3"))+
  ylab("Variance")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=15),legend.text=element_text(size=12),
        legend.title=element_text(size=12),legend.margin=unit(0.455,"cm"))+
  labs(fill = "G (Control)",color = "G (Control)")

load("plot4.RData")
plot4 <- ggplot(pos.fh2.b.sum, aes(x = factor, color=model,fill=model))+
  #geom_boxplot(alpha = 0.5)+
  geom_boxplot(alpha = 0.5,aes(ymin=q0,lower=q1,middle=q2,upper=q3,ymax=q4),stat = "identity")+
  scale_color_manual(breaks = c("star","new","old"),
                     values=c("grey60","grey30","grey0"),labels=c("F1","F2","F3"))+
  scale_fill_manual(breaks = c("star","new","old"),
                    values=c("grey60","grey30","grey0"),labels=c("F1","F2","F3"))+
  ylab("Variance")+
  theme(axis.text=element_text(size=10),
        axis.title=element_text(size=15),legend.text=element_text(size=12),
        legend.title=element_text(size=12),legend.margin=unit(0.455,"cm"))+
  labs(fill = "G (Control)",color = "G (Control)")
# percentage og heribility
load('BSFG_state.RData')
p = BSFG_state$run_variables$p

Lambda = BSFG_state$Posterior$Lambda
Lambda = matrix(Lambda[,ncol(Lambda)],nr=p)
old.G = G_Matrix_Comp(BSFG_state)

G_pos = function(G){
G = array(unlist(G),dim = c(nrow(G[[1]]),ncol(G[[1]]),length(G)))
G.pos = apply(G,c(1,2),mean)
return(G.pos)
}
old.G.pos = G_pos(old.G)
#the size of x1, x2, x3 column

#Contribution of ith Factor
con_factor = function(l,G.pos){
Lambda_l = t(Lambda[,l])%*%Lambda[,l]
S = sum(diag(G.pos))
return(Lambda_l/S)
}

con_factor(3,old.G.pos)
sapply(1:8,function(x)con_factor(x,old.G.pos))
# 2.For each trait, compare the correlation of G1 vs G2
source('~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG/BSFG_functions.R')
load(BSFGOld)
#old.traces.G = G_Traces_Comp(BSFG_state)
old.G = G_Matrix_Comp(BSFG_state)
p = BSFG_state$run_variables$p
# posterior mean
# do this in case the value of k differ for this two models and there in no lambda in BSFG_state_fixed case
Lambda = BSFG_state$Posterior$Lambda
load(BSFGNew)
load("BSFG_state7.RData")
load(BSFGStar)
BSFG_state$Posterior$Lambda = Lambda
#new.traces.G = G_Traces_Comp(BSFG_state)
new.G = G_Matrix_Comp(BSFG_state)
#For all the elements in G
#pairs(t(old.traces.G)~t(new.traces.G))
#For each trait

old.G = array(unlist(old.G),dim = c(nrow(old.G[[1]]),ncol(old.G[[1]]),length(old.G)))
new.G = array(unlist(new.G),dim = c(nrow(new.G[[1]]),ncol(new.G[[1]]),length(new.G)))
old.G.pos = apply(old.G,c(1,2),mean)
new.G.pos = apply(new.G,c(1,2),mean)
#new.G.b = apply(new.G,2,rbind)
#save(old.G.pos,file = "old_G_b.RData")
#load("old_G_b.RData")
#save(new.G.pos,file = "new_G_b.RData")
#load("new_G_b.RData")
#plot(old.G.b[,2],new.G.b[,2])
#combine old.G and new.G as a dataframe
#G.b = c(old.G.b[,2],new.G.b[,2])
#G.b = cbind(G.b, rep(c("old","new"),each=p*spn))
#G.b = as.data.frame(G.b)
#colnames(G.b) = c("covariance","model")
# 1st triat

#variance
par(mfrow=c(1,2))
plot(diag(new.G.pos),diag(old.G.pos),main = "Variance")
abline(0,1)
#covariance
diag(new.G.pos)=NA
diag(old.G.pos)=NA
plot(old.G.pos,new.G.pos, main="Covariance")
abline(0,1)
par(mfrow=c(1,1))
#each_trait_G_pair = function(i){
#G.b = cbind(old.G.b[,i],new.G.b[,i])
#G.b = as.data.frame(G.b)
#colnames(G.b) = c("Old","New")
#library(ggplot2)
#ggplot(G.b,aes(y = as.numeric(Old),x=as.numeric(New)))+
  #geom_dotplot(binwidth = 15)
#  geom_point()+
#  geom_abline(intercept = 0,slope = 1)+
#  ggtitle(paste0(BSFG_state$traitnames[i]))
#}
#pdf("each_trait_G_pair.pdf")
#par(mfrow = c(3,3))
#for (i in 1:p){
#  each_trait_G_pair(i)
#}
#dev.off()
#each_trait_G_pair(1)
#each_trait_G_pair(2)
#each_trait_G_pair(3)
#each_trait_G_pair(4)
#each_trait_G_pair(5)

# 3.delta z
source('~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG/BSFG_functions.R')
load(BSFGOld)
#old.traces.G = G_Traces_Comp(BSFG_state)
old.G = G_Matrix_Comp(BSFG_state)
p = BSFG_state$run_variables$p
# posterior mean
# do this in case the value of k differ for this two models and there in no lambda in BSFG_state_fixed case
Lambda = BSFG_state$Posterior$Lambda
load(BSFGNew)
load("BSFG_state6.RData")
load(BSFGStar)
BSFG_state$Posterior$Lambda = Lambda
#new.traces.G = G_Traces_Comp(BSFG_state)
new.G = G_Matrix_Comp(BSFG_state)

# The distance of two matrix(in predicted phenotype value/level)
#selection gradient
beta = 0.3609
beta.v = c(0,beta,rep(0,16))
#traces.traits = matrix(,p,spn)

#old.G in list form
old.traces.traits = sapply(old.G,function(x)x%*%beta.v)
old.traits.posmean = apply(old.traces.traits,1,mean)
new.traces.traits = sapply(new.G,function(x)x%*%beta.v)
new.traits.posmean = apply(new.traces.traits,1,mean)
star.traces.traits = sapply(star.G,function(x)x%*%beta.v)
star.traits.posmean = apply(star.traces.traits,1,mean)
#boxplot(old.traces_traits[2,])
boxplot(t(old.traces.traits))
#boxplot(new.traces_traits[2,])
boxplot(t(new.traces.traits))

#plot(x=rep(1:18,2),y=c(old.traits.posmean,new.traits.posmean),col=rep(c(1,2),each=18),type = "o")
save(traces_traits,file = 'traces_traits.RData')
# For each G matrix
setwd("~/Runcie Lab/SparseFactorMixedModel_v2/MCMC/19/19/Lambda1.5_delta2shape3")
load("G_BSFG_sc.RData")

ggplot(pos.fh2.b, aes(x = factor, y = Fh2, color=model,fill=model))+
  geom_boxplot(alpha = 0.5)+
  scale_color_manual(breaks = c("old", "new", "star"),
                     values=c("green", "red", "blue"))+
  scale_fill_manual(breaks = c("old", "new", "star"),
                    values=c("green", "red", "blue"))+
  ggtitle("The distribution of variance (Latent traits)")+
  scale_color_discrete(breaks = c("old", "new", "star"),labels=c("LconG3", "LhighG3", "LlowG3"))+
  scale_fill_discrete(breaks = c("old", "new", "star"),labels=c("LconG3", "LhighG3", "LlowG3"))+
  ylab("Variance")+
  geom_errorbar(aes(ymin=len-se, ymax=len+se), width=.1, position=pd)

#4. E(beta(G1-G2)beta)
source('~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG/BSFG_functions.R')
load(BSFGOld)
#old.traces.G = G_Traces_Comp(BSFG_state)
old.G = G_Matrix_Comp(BSFG_state)
p = BSFG_state$run_variables$p
# posterior mean
# do this in case the value of k differ for this two models and there in no lambda in BSFG_state_fixed case
Lambda = BSFG_state$Posterior$Lambda
load(BSFGNew)
load(BSFGStar)
BSFG_state$Posterior$Lambda = Lambda
#new.traces.G = G_Traces_Comp(BSFG_state)
load("BSFG_state6.RData")
new.G = G_Matrix_Comp(BSFG_state)
load("BSFG_state7.RData")
star.G = G_Matrix_Comp(BSFG_state)
beta = 0.3609
beta.v = as.vector(c(0,beta,rep(0,16)))
# calcualte mean diff response value
diff.traces = mapply(function(x,y) t(beta.v)%*%t(x-y)%*%(x-y)%*%beta.v, old.G,new.G)
boxplot(diff.traces)
mean(sqrt(diff.traces))

#5. PCA decomposition

source('~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG/BSFG_functions.R')
load(BSFGOld)
#old.traces.G = G_Traces_Comp(BSFG_state)
old.G = G_Matrix_Comp(BSFG_state)
p = BSFG_state$run_variables$p
traitnames = BSFG_state$traitnames
#load(BSFGNew)
load("BSFG_state2.RData")
new.G = G_Matrix_Comp(BSFG_state)
load("BSFG_state3.RData")
star.G = G_Matrix_Comp(BSFG_state)

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



save(traces.eigen.df.sum , file = "eigen567.RData")
save(traces.eigen.df.sum, file = "eigen51417.RData")

library(ggplot2)
pd <- position_dodge(width = 0.6)
load("eigen567.RData")

# size of each model
library(dplyr)
traces.eigen.df.sum %>%
  group_by(model) %>%
  summarize(size = sum(response))

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



load("eigen51417.RData")
# size of each model
library(dplyr)
traces.eigen.df.sum %>%
  group_by(model) %>%
  summarize(size = sum(response))

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

# 6. angle and length of deltaZ

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
old.theta = dist_within_theta(old.traces.traits)

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
old.new.theta = dist_between_theta(old.traces.traits,new.traces.traits)
# K-S test
ks.test(old.theta,old.new.theta) #different

old.star.theta = dist_between_theta(old.traces.traits,star.traces.traits)
new.star.theta = dist_between_theta(new.traces.traits,star.traces.traits)

theta.df = data.frame(old.theta,new.theta,star.theta, old.new.theta, old.star.theta, new.star.theta)

library(reshape2)
theta.df = melt(theta.df, value.name = "Degree")
names(theta.df) = c("labels", "theta")
ggplot(theta.df, aes(x = theta, col=labels))+
  geom_density(aes(fill=labels),alpha=0.1)
#+scale_color_manual(breaks = c("old.theta","new.theta","star.theta","old.new.theta","old.star.theta","new.star.theta"),
#                     values=rainbow(6))
  
# length
length_deltaZ = function(traces.traits){
    posmean.deltaZ = apply(traces.traits,1,mean)
    deltaZ.length = sqrt(sum(posmean.deltaZ * posmean.deltaZ)) 
    return(deltaZ.length)
}
old.deltaZ.length = length_deltaZ(old.traces.traits)
old.deltaZ.length
new.deltaZ.length = length_deltaZ(new.trace s.traits)
new.deltaZ.length
star.deltaZ.length = length_deltaZ(star.traces.traits)
star.deltaZ.length
