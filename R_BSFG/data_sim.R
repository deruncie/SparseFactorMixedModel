library(pedantics)
#library(R.matlab)


prepare_simulation2 = function(A_chol,A,Lambda,factor_h2s,E_a_prec,resid_Y_prec,name){
  n = dim(A_chol)[1]  #A_chol is  Cholesky decomposition of A
  p = nrow(Lambda) 
  
  r = n
  # n = 2*n
  
  k = ncol(Lambda)
  F_a = t(A_chol) %*% sweep(matrix(rnorm(n*k),n),2,sqrt(factor_h2s),'*') # make sure MN(0,A,sigma)
  F_r = sweep(matrix(rnorm(n*k),n),2,sqrt(1-factor_h2s),'*') # factor_h2s is the variance of factor
  F = F_a + F_r
  
  E_a_act = t(A_chol) %*% sweep(matrix(rnorm(n*p),n),2,sqrt(1/E_a_prec),'*') #
  E_r_act = sweep(matrix(rnorm(n*p),n),2,sqrt(1/resid_Y_prec),'*')
  
  Y = F %*% t(Lambda) + E_a_act + E_r_act
  
  X = matrix(1,nrow=1,ncol=n)
  B = rep(0,p)
  Z_1 = diag(r)[,rep(1:r,n/r)] 
  Va = colSums(Lambda^2) * factor_h2s + 1/E_a_prec  #diagonal elements in G
  Ve = colSums(Lambda^2) * (1-factor_h2s) + 1/resid_Y_prec #diagonal elements in R
  h2 = Va/(Va + Ve)   # heribility
  G = Lambda %*% diag(factor_h2s) %*% t(Lambda) + diag(1/E_a_prec) 
  R = Lambda %*% diag(1-factor_h2s) %*% t(Lambda) + diag(1/resid_Y_prec)
 # h2 = diag(G)/(diag(G)+diag(R))  
  
  # Z_2 = rbind(rep(c(1,0,0,0),n/4),rep(c(0,1,0,0),n/4),rep(c(0,0,1,0),n/4),rep(c(0,0,0,1),n/4))
  
  # plot(h2,type="l")
  
  print(name)
  setup = list(Y=Y,F_act=F,E_a_act=E_a_act,F_a_act = F_a,Z_1=Z_1,A=A,X=X,n=n,r=r,p=p,G=G,R=R,
               gen_factor_Lambda = Lambda,error_factor_Lambda = Lambda,h2=h2,B=B, factor_h2s=factor_h2s,name=name,pops)
  # ????????????????????gen_factor_Lambda ; error_factor_Lambda
  return(setup)
  #do.call(writeMat,c("setup.mat",setup))
  #save(pedigree,file="pedigree.Robj")
}

 
#==================================================================
counter = 0

#sample_size simulations
num_reps = 5
p = 100
k=10


#the same lambda
Lambda = matrix(0,p,k)
numeff = sample((p/30):(p/4),k,replace=T)
for(h in 1:k){   # even though numeff[h] is not a whole number, by applying sample, it only use the int part. It doesn't round off.
  Lambda[sample(1:p,numeff[h]),h] = rnorm(numeff[h])
}
Lambda = Lambda[,order(-diag(t(Lambda) %*% Lambda))]
cols=1:k
E_a_prec = rep(1/0.2,p)
resid_Y_prec = rep(1/0.2,p)
# pedigree from pengjuan's dataset to test if diff populations are indep
setwd("~/Runcie Lab/pedigree")
data=read.csv('IDDammSirelnUpdate(1).csv',header = TRUE)
source("~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG/setup_pedigree.R")
# find way to change this according to the size(different pops)

ori.setup = list()
ori.setup[[pops[1]]] = setup_pedigree(data,GenerationCode=3,LineCode=c(0,1),fixedfactor=NA) # no fixed factor because we seperate the dataset based on generation
ori.setup[[pops[2]]] = setup_pedigree(data,GenerationCode=2,LineCode=1,fixedfactor=NA)
ori.setup[[pops[3]]] = setup_pedigree(data,GenerationCode=1,LineCode=1,fixedfactor=NA)

setwd("~/Runcie Lab/SparseFactorMixedModel_v2/R_BSFG/Simulation")
 #sample_size = list(Small = c(50,5), Medium = c(100,10), Large = c(500,10))
# set seed for repeatability
set.seed(1)
pops = c("LConG1","LConG2","LConG3")
for(rep in 1:num_reps){
  # number of individuals that have the same mother:
  # c(50,5) in this case means that for 50 fathers, each of them have 5 children. There are total 250 individuals.
  # and for each individual, it has distinct mother
  setup = list()
  for(pop in pops){
    # nM = sample_size[[size]][1]	
    # nI = sample_size[[size]][2]	
    # pedigree = data.frame(ind=nM*nI + nM + 1:(nM*nI),dam=1:(nM*nI) + nM, sire = gl(nM,nI))
    # pedigree<-fixPedigree(pedigree)  #A
    # children = !is.na(pedigree[,3])
    # 
    # #generate A matrix as 2* kinship matrix from whole pedigree
    # A = 2*kinship(pedigree[,1], pedigree[,2], pedigree[,3])[children,children]
    
    # use pedigree from pengjuan's dataset
    A = ori.setup[[pop]]$A
    A_chol = chol(A)
    
    # this should be different for each pops
    factor_h2s = runif(k) # k numbers in [0-1]
    # order by descending order
    factor_h2s = factor_h2s[order(factor_h2s,decreasing = TRUE)]
    factor_h2s[seq(1,k,by=2)]=0  # make 5 of them to be zero
    
    g_cols = factor_h2s>0  # will this step makes lambda matrix to be different in pops?????
    #Lambda =Lambda[do.call("order", unname(split(-abs(Lambda[,cols[factor_h2s>0]]), col(Lambda[,cols[factor_h2s>0]])))),]
    # reordering lambda
    
    name = paste(pop,rep,sep='_')
    setup[[pop]]= prepare_simulation2(A_chol,A,Lambda,factor_h2s,E_a_prec,resid_Y_prec,name)
   
  }
  counter = counter+1
  dir = paste('Sim',counter,sep="_")
  try(dir.create(dir))
  setwd(dir)
  save(setup,file = 'setup.RData')
  setwd('..')
}



