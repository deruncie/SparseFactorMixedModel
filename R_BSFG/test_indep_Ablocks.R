setwd("~/Runcie Lab/pedigree")
data=read.csv('IDDammSirelnUpdate2.csv',header = TRUE)
# 5
GenerationCode = c(0,1,2,3)
LineCode = c(1)


# 6
GenerationCode=3
LineCode = c(0,2)


# 7
GenerationCode=3
LineCode = c(0,3)

ped <- data[,1:3]
for(i in 1:3) ped[,i] <- as.factor(ped[,i])
library(MCMCglmm)
invA <- inverseA(ped)
A <- solve(invA$Ainv)

data$Gen_Code = paste(data$Generation,data$LineCode,sep='_')
o = order(data$Gen_Code)

A_order = A[o,o]
image(A_order>0)


generation_pos <- data$Generation %in% GenerationCode
LineCode_pos <- data$LineCode%in%LineCode|data$LineCode==0
data_pos <- which((generation_pos&LineCode_pos)|data$LineCode==0)
#nodes_names <- as.factor(as.character(data$animal[data_pos]))
n = length(data_pos)

# reorder the data according to specific rows in each model
data.trans = data
data.trans[1:n,] = data[data_pos,]
data.trans[185:(184+n),] = data[data_pos,]
data.trans[366:(365+n),] = data[data_pos,]


ped <- data.trans[,1:3]
ped <- data[,1:3]
for(i in 1:3) ped[,i] <- as.factor(ped[,i])
library(MCMCglmm)
invA2 <- inverseA(ped)#, node = data.trans$animal[1:545])
invA <- inverseA(ped, node = data.trans$animal[1:545])
A <- solve(invA$Ainv)
A <- matrix(A@x, nrow<-n, ncol<-n, byrow <- FALSE)
