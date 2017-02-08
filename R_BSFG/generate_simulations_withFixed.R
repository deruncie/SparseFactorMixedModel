# simulator for various fixed effect models


simulate_data = function(
							Lambda,			# p x K matrix
							X,				# n x b matrix
							Z1,				# n x r matrix
							A1_chol,		# r x r upper triangle chol matrix
							F_vars,			# K x 3 matrix of fixed, A and E vars for factors
							E_h2,			# p vector of residual heritabilities
							mu,				# p vector of global means
							B,				# b x p residual fixed effects
							B_F	,			# b x K factor fixed effects
							sigma2_Resid	# p x 1 non-factor variances
						) {
  n = nrow(X)
  r = ncol(Z1)
  p = nrow(Lambda)
  K = ncol(Lambda)
  
  E_A1 = t(A1_chol) %*% matrix(rnorm(r*p),r,p) %*% diag(E_h2)
  F_A1 = t(A1_chol) %*% matrix(rnorm(r*k),r,k) %*% diag(F_vars[,2])
  
  E_R = matrix(rnorm(n*p),n,p) %*% diag(1-E_h2)
  F_R = matrix(rnorm(n*k),n,k) %*% diag(F_vars[,3])
  
  E = X %*% B + Z1 %*% E_A1 + E_R
  F = X %*% B_F + Z1 %*% F_A1 + F_R
  Y = matrix(1,n,1) %*% matrix(mu,1,p) + F %*% t(Lambda) + E %*% diag(sqrt(sigma2_Resid))
  
  return(Y)
}

save_sim = function(
						sim_ID = sim_ID,
						Lambda = Lambda,
						X = X,
						Z1 = Z1,
						A1_chol = A1_chol,
						F_vars = F_vars,
						E_h2 = E_h2,
						mu = mu,
						B = B,
						B_F = B_F,
						sigma2_Resid = sigma2_Resid
					){

	F_h2 = F_vars[,2]/(F_vars[,2]+F_vars[,3])
	sigma2_Resid = rgamma(p,3,6)

	Y = simulate_data(Lambda,X,Z1,A1_chol,F_vars,E_h2,mu,B,B_F,sigma2_Resid)
	sim_data = list(
					n      = n,
					k      = k,
					p      = p,
					b      = ncol(X),
					r      = n_Lines,
					Y      = Y,
					X      = X,
					Z1     = Z1,
					A1     = A1,
					B      = B,
					B_F    = B_F,
					mu     = mu,
					F_vars = F_vars,
					F_h2   = F_h2,
					E_h2   = E_h2,
					Lambda = Lambda,
					G      = Lambda %*% diag(F_h2) %*% t(Lambda) + diag(sigma2_Resid*E_h2),
					R      = Lambda %*% diag(1-F_h2) %*% t(Lambda) + diag(sigma2_Resid*(1-E_h2))
			)

	try(dir.create(sim_ID))
	save(sim_data,file=paste(sim_ID,'Sim_data.RData',sep='/'))
}

sim_reps = 1

n_Lines = 100
n_Reps = 3
n_TRT = 2
n = n_Lines*n_Reps*n_TRT
k = 4
p = 100

# design matrices
sample_info = data.frame(Line = gl(n_Lines,n_Reps*n_TRT),TRT = rep(c(0,1),each = n_Reps), Rep = c(1:n_Reps))
Z1 = model.matrix(~Line+0,sample_info)
X = model.matrix(~TRT+0,sample_info)
A1 = diag(1,ncol(Z1))
A1_chol = chol(A1)

Lambda = cbind(
	c(rep(1,floor(p/3)),rep(0,p-floor(p/3))),
	c(rep(0,floor(p/6)),rep(1,floor(p/4)),rep(0,p-floor(p/6)-floor(p/4))),
	c(rep(0,floor(p/3)),rep(1,floor(p/4)),rep(0,p-floor(p/4)-floor(p/3))),
	c(rep(0,floor(p/3)+floor(p/4)),rep(1,floor(p/8)),rep(0,p-floor(p/3)-floor(p/4)-floor(p/8))),
	c(rep(0,floor(p/3)+floor(p/4)),rep(0,floor(p/8)),rep(1,floor(p/8)),rep(0,p-floor(p/3)-floor(p/4)-2*floor(p/8))),
	c(rep(0,floor(p/3)+floor(p/4)),rep(0,2*floor(p/8)),rep(1,floor(p/8)),rep(0,p-floor(p/3)-floor(p/4)-3*floor(p/8))),
	c(rep(0,floor(p/3)+floor(p/4)),rep(0,3*floor(p/8)),rep(1,p-floor(p/3)-floor(p/4)-3*floor(p/8)))
	)
Lambda = Lambda[,1:k]
image(Lambda)
colSums(Lambda > 0)

for(rep in 1:sim_reps){

	# Set 1. Dense fixed effects, moderate heritability
	sim_ID = paste0('Sim_FE1_',rep)
	mu = rnorm(p)
	B = matrix(rep(c(1,-1),p/2),nr = 1)
	B_F = matrix(rep(0,k),nr = 1)
	F_vars = cbind(rep(0,k),sample(seq(.1,.7,length=k)))
	F_vars = cbind(F_vars,1-rowSums(F_vars))
	E_h2 = runif(p,0,.7)

	save_sim(sim_ID, Lambda, X, Z1, A1_chol, F_vars, E_h2, mu, B, B_F, sigma2_Resid)

	# Set 2. Factor fixed effects
		# first factor is 100% fixed. Remainder have 0% fixed

	sim_ID = paste0('Sim_FE2_',rep)
	mu = rnorm(p)
	B = matrix(0,nr = 1,nc = p)
	B_F = matrix(c(1,rep(0,k-1)),nr = 1)
	F_vars = cbind(c(1,rep(0,k-1)),c(0,sample(seq(.1,.7,length=k-1))))
	F_vars = cbind(F_vars,1-rowSums(F_vars))
	F_vars[1,3] = 0
	E_h2 = runif(p,0,.7)

	save_sim(sim_ID, Lambda, X, Z1, A1_chol, F_vars, E_h2, mu, B, B_F, sigma2_Resid)



	# Set 3. 1 joint fixed and random factor at 50% each. Rest only random

	sim_ID = paste0('Sim_FE3_',rep)
	mu = rnorm(p)
	B = matrix(0,nr = 1,nc = p)
	which_fixed = sample(1:k,1)
	B_F = matrix(0,nr = 1,nc = k)
	B_F[1,which_fixed] = 1
	F_vars = cbind(rep(0,k),sample(seq(.1,.7,length=k)))
	F_vars = cbind(F_vars,1-rowSums(F_vars))
	F_vars[which_fixed,1] = 0.5
	F_vars[which_fixed,-1] = 0.5 * F_vars[which_fixed,-1]

	save_sim(sim_ID, Lambda, X, Z1, A1_chol, F_vars, E_h2, mu, B, B_F, sigma2_Resid)

	# Set 4. 2 joint fixed and random factors at 50% each. Rest only random

	sim_ID = paste0('Sim_FE4_',rep)
	mu = rnorm(p)
	B = matrix(0,nr = 1,nc = p)
	which_fixed = sample(1:k,2)
	B_F = matrix(0,nr = 1,nc = k)
	B_F[1,which_fixed] = 1
	F_vars = cbind(rep(0,k),sample(seq(.1,.7,length=k)))
	F_vars = cbind(F_vars,1-rowSums(F_vars))
	F_vars[which_fixed,1] = 0.5
	F_vars[which_fixed,-1] = 0.5 * F_vars[which_fixed,-1]

	save_sim(sim_ID, Lambda, X, Z1, A1_chol, F_vars, E_h2, mu, B, B_F, sigma2_Resid)

}