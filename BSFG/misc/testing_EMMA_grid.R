library(emma)
rs <- emma.REML.t(emmadat$ys,emmadat$xs,emmadat$K)

h2s = rs$vgs/(rs$vgs + rs$ves)
boxplot(h2s)

j=22
ps = sapply(1:dim(emmadat$xs)[1],function(i) {
	vg = rs$vgs[i,j]
	ve = rs$ves[i,j]
	Sigma = vg*emmadat$K + ve*diag(1,nrow(emmadat$K))
	R = chol(Sigma)
	res = solve(t(R),cbind(c(emmadat$ys[j,]),1,c(emmadat$xs[i,])))
	Rity = res[,1]
	Ritx = res[,-1]
	lm1 = lm(Rity~0+Ritx)
	return(summary(lm1)$coef[2,4])
	})
plot(ps,rs$ps[,j],log='xy');abline(0,1)

ngrid=7
K = emmadat$K
center_K = function(K){
  n = nrow(K)
  J = matrix(1,n,1)
  M = diag(1,n) - J %*% solve(t(J) %*% J) %*% t(J)
  centered_K = M %*% K %*% t(M)
  rownames(centered_K) = rownames(K)
  colnames(centered_K) = colnames(K)
  centered_K
}
K = center_K(K)
# K[K<0]=0
# K = as.matrix(drop0(Matrix(abs(K)),tol = 1e-2))
# K = diag(1,nrow(K))
ps2 = sapply(1:dim(emmadat$xs)[1],function(i) {
	vg = rs$vgs[i,j]
	ve = rs$ves[i,j]
	vp = vg+ve
	h2 = vg/(vg+ve)
	# h2 = h2 / 1.5
	# h2 = 0.95
	# h2 = round(h2,digits=1)
	h2 = round(h2*ngrid)/ngrid
	vg = h2*vp
	ve = (1-h2)*vp
	Sigma = vg*K + ve*diag(1,nrow(emmadat$K))
	R = chol(Sigma)
	res = solve(t(R),cbind(c(emmadat$ys[j,]),1,c(emmadat$xs[i,])))
	Rity = res[,1]
	Ritx = res[,-1]
	lm1 = lm(Rity~0+Ritx)
	return(summary(lm1)$coef[2,4])
	})
plot(ps,ps2,log='xy');abline(0,1)