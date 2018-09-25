library(BSFG)
library(Matrix)
n = 3
p = 5
k = 2
zs = lapply(1:(p+k),function(x) tcrossprod(rstdnorm_mat(n,n)))

f = matrix(1:(n*k),ncol=k)
y = matrix(1:(n*p),ncol=p)

cov_vf = bdiag(zs[[p+1]],zs[[p+2]])
cov_vy = bdiag(zs[1:p])

cov_vft = cov_vf[c(t(f)),c(t(f))]
cov_vyt = cov_vy[c(t(y)),c(t(y))]

I = diag(1,n)
L = rstdnorm_mat(p,k)

X = kronecker(I,L)
t(X) %*% solve(cov_vyt) %*% X + solve(cov_vft)

cov_vyt + X %*% cov_vft %*% t(X)

qrL = qr(L)
qr.qty(qrL,L)
R = drop0(qr.qty(qrL,L)[1:k,],tol=1e-10)
y2 = qr.qty(qrL,t(y))[1:k,]

Q = qr.Q(qrL)
sigma2 = bdiag(lapply(1:n,function(x) t(Q))) %*% cov_vyt %*% bdiag(lapply(1:n,function(x) (Q)))
