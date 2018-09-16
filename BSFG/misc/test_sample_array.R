
sample_array = function(data = NA, dim, dimnames = NULL){
  X = array(data,dim[c(2,3,1)],dimnames)
  class(X) = c('sample_array',class(X))
  X
}

`[<-.sample_array` = function(X,i,j,k,value){
  if(missing(j) && missing(k) && length(i) == 1){
    if(!storage.mode(X) == 'double') storage.mode(X) = 'double'
    if(!storage.mode(value) == 'double') storage.mode(value) = 'double'
    addSample(X,value,dim(X),i)
    return(X)
  } else{
    class(X) = 'array'
    X[j,k,i] <- value
    class(X) = c('sample_array','array')
    return(X)
  }
}

`[.sample_array` = function(X,i,j,k) {
  if(missing(i)) i = 1:dim(X)[1]
  if(missing(j)) j = 1:dim(X)[2]
  if(missing(k)) k = 1:dim(X)[3]
  # `[.default`(X,j,k,i)
  # .Primitive('[')(X,i=j,j=k,i)
  # NextMethod('[',X,i=i,j,j,k)
  # recover()
  class(X) = class(X)[-1]
  if(!missing(i) && (length(i) == 1 || sum(i) == 1)) {
    res = X[,,i]
    return(res[j,k,drop=FALSE])
  } else{
    res = X[j,k,i,drop=FALSE]
    res = aperm(res,c(3,1,2))
    if(dim(res)[2] == 1) {
      res = res[,1,]
      return(res)
    }
    if(dim(res)[3] == 1) {
      res = res[,,1]
      return(res)
    }
    return(res)
  }
}
`[.array` = function(X,i,j,k) {
  print(c(i,j,k))
  NextMethod('[',X,i=i,j,k)
}

`dimnames.sample_array` = function(X) {
  x = NextMethod()
  x[c(3,1,2)]
}
`dimnames<-.sample_array` = function(X,x){
  # print('start')
  if(length(x) != 3) stop("require three dimensions")
  attr(X,'dimnames') = x[c(2,3,1)]
}

dim.sample_array = function(X){
  dim = NextMethod()
  dim[c(3,1,2)]
}

print.sample_array = function(X){
  print(aperm(X,c(3,1,2)))
}


Xa = array(0,c(50,10000,100))
Xs = sample_array(0,c(50,10000,100))
Xs2 = Xs
class(Xs2) = 'array'
attr(Xs2,'class') = 'array'

Lambda = matrix(rnorm(10000*100),10000)
rbenchmark::benchmark(
Xa[50,,] <- matrix(rnorm(10000*100),10000),
# Xs[50,,] <- Lambda,
Xs2[,,50] <- matrix(rnorm(10000*100),10000),
addSample(Xs,matrix(rnorm(10000*100),10000),dim(Xs),50),
replications=100)
