library(GridLMM)
library(BSFG)
library(microbenchmark)

n = 200
b = 20

X = rstdnorm_mat(n,b)
R = tcrossprod(rstdnorm_mat(n,n))+diag(1,n)
cholL_R = t(chol(R))
cholL_R_inv = solve(cholL_R)

microbenchmark(
# Si1 <- solve(tcrossprod(X) + R),
Si2 <- chol(crossprod(t(t(cholL_R_inv)) %*% X) + diag(1,b)),
# Si2 <- t(cholL_R_inv) %*% solve(tcrossprod(cholL_R_inv %*% X) + diag(1,n)) %*% cholL_R_inv,
# Si3 <- crossprod(forwardsolve(chol_update_L(cholL_R,X,rep(1,b)),diag(1,n))),
Si3 = chol_update_L(cholL_R,X,rep(1,b)),
times=10)
