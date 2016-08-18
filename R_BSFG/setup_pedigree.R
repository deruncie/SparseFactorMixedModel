setup_pedigree = function(data=data,LineCode){
  ## create pedigree table:
  # we use individuals from lines (F1 to F3 generation) and the parental ind.
  linecode_pos = which(data$LineCode==0|data$LineCode==LineCode)
  ped = data[linecode_pos, 1:3]
  Y = data[linecode_pos,-(1:6)]
  n = nrow(Y)
  Z_1 = diag(1,n,n)
  
  for(i in 1:3) ped[,i] = as.factor(ped[,i])
  #head(ped)
  #dim(ped)
  library(MCMCglmm)
  invA <- inverseA(ped)
  A <- solve(invA$Ainv)
  #str(A)
  #max(A)
  #min(A)
  A <- matrix(data=A@x, nrow=n, ncol=n, byrow = FALSE)
  #Amatrix
  #write.csv(Amatrix, file = "Amatrix.csv")
  

  #X = as.factor(data$Generation[data$LineCode<=1])
  X = as.factor(data$Generation[linecode_pos])
  #dummyX = matrix(c(rep(1,n),rep(0,3*n)),n,4)
  #dummyX[which(X==0),2] = 1 
  #dummyX[which(X==1),3] = 1
  #dummyX[which(X==2),4] = 1
  #X = dummyX
  #rm(dummyX)
  
  #reg=lm(as.matrix(Y)~0+X,na.action = na.pass)
  
  X = model.matrix(~X)
  #reg1=lm(as.matrix(Y)~X)
  #B_act_1 = reg1$coefficients
  reg=lm(as.matrix(Y)~X)
  B_act = reg$coefficients
  
  traitnames = names(data)[7:ncol(data)]
  setup = list(Y,A,Z_1,X,B_act,traitnames)
  names(setup) = c("Y","A","Z_1","X","B_act","traitnames")
  save(setup,file=sprintf("setup_LineCode%d.RData",LineCode))
  return(setup)
}
