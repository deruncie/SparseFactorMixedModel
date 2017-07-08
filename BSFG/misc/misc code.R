files = list.files()
new_files = sub(' (deruncie@ucdavis.edu)','',files,fixed=T)
for(i in 1:length(files)){
  system(sprintf('mv "%s" "%s"',files[i],new_files[i]))

}



cl <- makeCluster(1, type = "PSOCK")

BSFG_state$current_state$Lambda_ncores = 8
BSFG_state$current_state$B_ncores = 8
BSFG_state$current_state$B_F_ncores = 8

f = function() return(matrix(rnorm(1:1e8),1e6))


for(i in 1:15){
as = 2^seq(-5,5,length=100)
res=sapply(as,function(a) {
  with(BSFG_state$current_state,{
    sl = sum(dnorm(a*Lambda[,i],0,1/sqrt(Plam[,i]),log=T))
    fr = F - BSFG_state$data_matrices$X_F %*% B_F - U_F
    sr = sum(dnorm(1/a*F[,i],0,sqrt((1-F_h2[i])/tot_F_prec[i]),log=T))
    sl+sr
  })
})
plot(log2(as),log2(-res+max(res)+5),main=i)
}
