## See Thompson and Meyer 1986 for Mathematical treatment.

p = 3
n = 1

G = diag(1,p) + 1
R = diag(1,p) +1
K = diag(.5,n) + .5


g = sort(rexp(p),decreasing = T)
r = runif(p)
r[p] = r[p] * 2
U = svd(matrix(rnorm(p*p),p))$u
G = U %*% diag(g) %*% t(U)
R = U %*% diag(r) %*% t(U)

# G = diag(g)
# R = diag(r)

vars = function(G,R,K){
  p = ncol(G)
  n = ncol(K)
  S1 = solve(diag(1,n)/R[1,1]+solve(K)/G[1,1])
  S2 = solve(kronecker(solve(R),diag(1,n)) + kronecker(solve(G),solve(K)))[1:n,1:n]
  S3 = (kronecker(U,diag(1,n)) %*% solve(kronecker(diag(1/r),diag(1,n)) + kronecker(diag(1/g),solve(K))) %*% kronecker(t(U),diag(1,n)))[1:n,1:n]
  (kronecker(U[1,,drop=F],diag(1,n)) %*% solve(kronecker(diag(1/r),diag(1,n)) + kronecker(diag(1/g),solve(K))) %*% kronecker(t(U[1,,drop=F]),diag(1,n)))

  matrix(rowSums(sapply(1:p,function(j) U[,j,drop=F] %*% (1/(1/r[j] + 1/g[j])) %*% t(U[,j,drop=F]))),p)
  sum(sapply(1:p,function(j) U[1,j]^2 * (1/(1/r[j] + 1/g[j])) ))
  sum(sapply(1:p,function(j) U[1,j]^2 * r[j] ))

  U[1,,drop=F] %*% diag(1/(1/r + 1/g)) %*% t(U[1,,drop=F])
  1/(1/(U[1,,drop=F] %*% diag(r) %*% t(U[1,,drop=F])) + 1/(U[1,,drop=F] %*% diag(g) %*% t(U[1,,drop=F])))

  return(c(S1[1],S2[1]))
}

# ## Rvar
# x = 1:10
# res = sapply(x,function(x) {
#   R[1,1] = x
#   vars(G,R,K)
# })
# plot(x,log2(res[2,]/res[1,]))


## p
ps = seq(2,22,length=20)
res = sapply(ps,function(p) {

  G = diag(1,p) + .7
  R = diag(1,p) + .7
  K = diag(.5,n) + .5
  R[1,1] = 2
  # G[1,2]=G[2,1] = .5*G[1,1]
  vars(G,R,K)
})

plot(ps,log2(res[2,]/res[1,]))
# plot(ps,res[2,])
