
setup_pedigree = function(data=data,GenerationCode,LineCode,fixedfactor,randomfactor=TRUE){  
  
  generation_pos <- data$Generation %in% GenerationCode
  LineCode_pos <- data$LineCode%in%LineCode|data$LineCode==0
  data_pos <- (generation_pos&LineCode_pos)|data$LineCode==0
  nodes_names <- as.factor(as.character(data$animal[data_pos]))
  
  Y <- data[data_pos,-(1:6)]
  ped <- data[,1:3]
  n <- nrow(Y)
  
  Z_1 <- diag(1,n,n)
  for(i in 1:3) ped[,i] <- as.factor(ped[,i])
  library(MCMCglmm)
  if(all(GenerationCode==c(0,1,2,3))&all(LineCode==c(0,1,2,3))){
    invA <- inverseA(ped)
  }else{invA <- inverseA(ped,nodes<- nodes_names)}
  A <- solve(invA$Ainv)
  A <- matrix(A@x, nrow<-n, ncol<-n, byrow <- FALSE)
  
  if (randomfactor == TRUE){
    Z_2 <- model.matrix(~0+as.factor(data[data_pos,2]))
  }else{Z_2=matrix(0,nr = n,nc = 0)}
  
  if(!is.na(fixedfactor)){
  if(fixedfactor=="Generation"){
    X = as.factor(data$Generation[data_pos])
    X = model.matrix(~X,data = data[data_pos,])
  }else if(fixedfactor=="Generation,Line"){
    X1 = as.factor(data$Generation[data_pos])
    X2 = as.factor(as.character(data$Line[data_pos]))
    X = model.matrix(~X1+X2,data = data[data_pos,])
  }
    reg=lm(as.matrix(Y)~0+X, na.action= na.omit)
    B_act = reg$coefficients
  }else{X=matrix(0,nr=n,nc=0);B_act=matrix(0,nr=0,nc=ncol(Y))}
  
  
  traitnames <- names(data)[7:ncol(data)]
  setup <- list(Y,X,B_act,A,Z_1,Z_2,traitnames)
  names(setup) <- c("Y","X","B_act","A","Z_1","Z_2","traitnames")
  
  save(setup,file="setup.RData")
  return(setup)
}
