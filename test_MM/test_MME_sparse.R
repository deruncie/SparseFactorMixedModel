library(Matrix)
library(SparseM)
library(microbenchmark)
library(Rcpp)
library(RcppArmadillo)
sourceCpp('MME_solve.cpp')
n = 200
A1 = sapply(1:5,function(x) sample(c(0,1),size=n,prob = c(.9,.1),replace=T)); A1 = A1 %*% t(A1) + diag(1,n)
A2 = sapply(1:5,function(x) sample(c(0,1),size=n,prob = c(.9,.1),replace=T)); A2 = A2 %*% t(A2) + diag(1,n)
A = A1+A2

x = rnorm(n)

solve1 = function(A,x){
	solve(A) %*% x
}
solve2 = function(A,x) {
	solve(A,x)
}
solve3 = function(cA,x) {
	x_star = forwardsolve(cA,x)
	forwardsolve(cA,backsolve(cA,x))
}
solve4 = function(cA,x) {
	x_star = c(as.matrix(solve(t(cA),x)))
	solve((cA),x_star)
}
solve5 = function(cAi,x) {
	Ai = cAi %*% t(cAi)
	Ai %*% x
}
solve6 = function(cA,x) {
	solve(cA,x)
}
cA = chol(A)
cA1 = SparseM::chol(A)
cA2 = chol(Matrix(A))
cA2i = solve(cA2)
cA3 = Cholesky(Matrix(A))
As = Matrix(A,sparse=T)
A1s = Matrix(A1,sparse=T)

max(abs(solve1(A,x) - solve2(A,x)))
max(abs(solve1(A,x) - solve3(cA1,x)))
max(abs(solve1(A,x) - solve4(cA2,x)))
max(abs(solve1(A,x) - solve5(cA2i,x)))
max(abs(solve1(A,x) - solve6(cA3,x)))
max(abs(solve1(A,x) - solve_sp_c(cA2,x)))
max(abs(solve1(A,x) - solve_c(cA,x)))

microbenchmark(solve1(A,x),solve2(A,x),solve4(cA2,x),solve5(cA2,x),solve6(cA3,x),solve_c(cA,x))
microbenchmark(chol(10*As1),update(cA3,10*As))


max(abs(solve4(chol(10*A1s),x) - solve6(update(cA3,10*A1s),x)))
microbenchmark(solve4(chol(10*A1s),x),solve6(update(cA3,10*A1s),x))
microbenchmark(det(As),det(update(cA3,As),x)^2)

max(abs(solve4(chol(A1),x) - solve_c(cA,x)))