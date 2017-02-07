p = ncol(BSFG_state$data$Y)
n = nrow(BSFG_state$data$Y)
Posterior = BSFG_state$Posterior
k = nrow(Posterior$Lambda) / p

j = 1
i_p = (j-1)*p + 1:p
i_n = (j-1)*n + 1:n

Lj = Posterior$Lambda[i_p,]
Fj = Posterior$F[i_n,]

par(mfrow=c(1,1))


par(mfrow=c(4,4))
for(j in 1:15){
i_p = (j-1)*p + 1:p
i_n = (j-1)*n + 1:n

Lj = Posterior$Lambda[i_p,]
Fj = Posterior$F[i_n,]
plot(apply(Fj,2,var),apply(abs(Lj),2,max))
abline(v=1)
}
par(mfrow=c(1,1))
plot(apply(abs(Lj),2,max),type='l')
plot(apply(Fj,2,var),type='l')
acf(apply(Fj,2,var))

varFs = sapply(1:k,function(j) {
	i_n = (j-1)*n + 1:n
	apply(Posterior$F[i_n,],2,var)
	}
	)
boxplot(varFs)

varPs = sapply(1:p,function(pp){
	i = pp+((1:k)-1)*p
	colSums(Posterior$Lambda[i,]^2) + 1/Posterior$E_a_prec[pp,] + 1/Posterior$resid_Y_prec[pp,]
})
boxplot(varPs)
hist(colMeans(varPs))
qqnorm(varPs[5,])

BSFG_state_1 = BSFG_state

BSFG_state = clear_Posterior(BSFG_state)
BSFG_state$current_state$F = .5*BSFG_state$current_state$F
BSFG_state$current_state$F_a = .5*BSFG_state$current_state$F_a
BSFG_state$current_state$F_px = .5*BSFG_state$current_state$F_px

BSFG_state$current_state$delta[1] = .5*BSFG_state$current_state$delta[1]
BSFG_state$current_state$Plam = .5*BSFG_state$current_state$Plam