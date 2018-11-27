library(BSFG)
library(GridLMM)
library(lme4qtl)
library(lme4)


n = 200
r = 2
p = 1000
X = rstdnorm_mat(n,p)
rownames(X) = 1:nrow(X)
K = tcrossprod(X)
K = K/mean(diag(K))
# diag(K) = diag(K)+1
# K = K/2

sK = svd(K)
U = sK$u
d = sK$d


h2 = 0.1
s2e = 30
beta = c(1,rep(0,p-1))

data = data.frame(ID = rownames(X))
data_full = data[rep(1:n,each=r),,drop=FALSE]
data_full$geno = data_full$ID
Z = model.matrix(~0+ID,data_full)
Z = Z[,paste0('ID',data$ID)]
data$Ut1 = t(U) %*% matrix(1,n,1)
D = diag(d)
rownames(D) = colnames(D) = data$ID

j = 1
res = c()
for(i in 1:100){
  print(i)
  r1 = U %*% (sqrt(d) * rnorm(n))
  y = X[as.character(data_full$ID),,drop=FALSE] %*% beta + sqrt(s2e)* Z %*% (sqrt(h2) * r1)+ sqrt(s2e)*sqrt(1-h2)*rnorm(n*r)
  ymean = aggregate(y~ID,data_full,FUN=mean)
  ymean = ymean[match(data$ID,ymean$ID),2]
  data$Uty_mean = t(U) %*% ymean
  m0 = relmat_lmer(y ~ (1|ID)+(1|geno),data_full,relmat = list(geno = K))
  data$Uty_fit = t(U) %*% (ranef(m0)[[1]][data$ID,1] + ranef(m0)[[2]][data$ID,1])
  m1f = relmatLmer(Uty_mean~t(U) %*% X[,j] + Ut1 + (1|ID),data,relmat = list(ID = D),REML=F)
  m1r = relmatLmer(Uty_mean~0 + Ut1 + (1|ID),data,relmat = list(ID = D),REML=F)
  m2f = relmatLmer(Uty_fit ~ t(U) %*% X[,j] + Ut1 + (1|ID),data,relmat = list(ID = D),REML=F)
  m2r = relmatLmer(Uty_fit ~ + Ut1 + (1|ID),data,relmat = list(ID = D),REML=F)
  r1 = anova(m1f,m1r)$Pr[2]
  r2 = anova(m2f,m2r)$Pr[2]
  res = rbind(res,c(r1,r2,m2f@optinfo$conv$opt+m2r@optinfo$conv$opt))
}
i = res[,3]==0
qqplot(-log10(res[i,1]),-log10(res[i,2]));abline(0,1)
