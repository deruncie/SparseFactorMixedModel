CovToCor = function(X){
	# Normalize a covariance matrix into a correlation matrix
    corX=diag(1/sqrt(diag(X))) %*% X %*% diag(1/sqrt(diag(X)));
    return(corX)
}

trace_plot = function(data,main = NULL,ylim = NULL){
  if(is.null(ylim)) ylim = range(data)
  plot(NA,NA,xlim = c(0,ncol(data)),ylim = ylim,main= main,xlab = "iteration")
	for(i in 1:nrow(data)) lines(data[i,],col=i)
}

draw_simulation_diagnostics = function(sp_num,run_parameters,run_variables,Posterior,Lambda,F_h2,E_a_prec,resid_Y_prec){

    devices = dev.list()
    while(length(devices) < 4){
        if(.Platform$OS.type != "windows") {
    	   quartz()
        } else {
            windows()
        }
    	devices = dev.list()
    }


    p = run_variables$p
    E_h2 = 1-F_h2;
    G_Lambda = sweep(Lambda,2,sqrt(F_h2),'*')
    E_Lambda = sweep(Lambda,2,sqrt(E_h2),'*')
    actual_G_Lambda = run_parameters$setup$gen_factor_Lambda
    actual_E_Lambda = run_parameters$setup$error_factor_Lambda
    G_est = G_Lambda %*% t(G_Lambda) + diag(1/E_a_prec)
    E_est = E_Lambda %*% t(E_Lambda) + diag(1/resid_Y_prec)
    G_act = run_parameters$setup$G
    E_act = run_parameters$setup$R
    h2s_act = run_parameters$setup$h2

    dev.set(devices[1])
    par(mfrow=c(3,3))

    cors = abs(cor(Lambda,actual_E_Lambda))
    cors=rbind(cors,apply(cors,2,max))
    cors = cbind(cors,apply(cors,1,max))
    image(cors[,ncol(cors):1],zlim=c(0,1),main = "cor of Lambda vs act_E_lambda")

    image(CovToCor(G_est),zlim=c(-1,1),main = "G_est")
    image(CovToCor(E_est),zlim=c(-1,1), main = "E_est")

    clim=max(0.1,min(.75,max(max(abs(G_Lambda)))))
    clims=c(-clim,clim)
    image(t(G_Lambda),zlim=clims,main = "G_Lambda")

    image(t(actual_G_Lambda),zlim=clims,main = "act_G_Lambda")

    clim=max(0.1,min(.75,max(max(abs(E_Lambda)))))
    clims=c(-clim,clim)
    image(t(E_Lambda),zlim=clims,main = "E_lambda")

    plot(c(G_est),c(G_act),main = "G_est vs G_act");abline(0,1)

    plot(c(E_est),c(E_act),main = "E_est vs E_act");abline(0,1)

    plot(F_h2,type='l',ylim=c(0,1),main = "F_h2 vs E_h2")
    lines(E_h2,col=2)

    if(sp_num>1){

    	dev.set(devices[2])
        # Figure of trace plots of the largest elements of each column of the factor loading matrix
        f2_row=4;
        f2_col=4;
        par(mfrow=c(f2_row,f2_col))
        for(k in 0:(min(2*f2_row,nrow(Posterior$Lambda)/p)-1)) {
            o = order(-abs(rowMeans(Posterior$Lambda[(1:p)+k*p,max(1,sp_num-100):sp_num])))
            traces = Posterior$Lambda[o[1:5]+k*p,1:sp_num]
            trace_plot(traces)
            
        }

    	dev.set(devices[3])
    	# Figure of trace plots and histograms of the factor heritablities
        f4_row=4;
        f4_col=4;
        par(mfrow=c(f4_row,f4_col))
        for(k in 1:min(2*f4_row,nrow(Posterior$F_h2))) {
            h2s = Posterior$F_h2[k,1:sp_num]
            if(sum(h2s[!is.na(h2s)])==0) {
                next
            }
            plot(h2s,type='l',main = "plot of h2s")
            hist(h2s,breaks=100,xlim=c(-0.1,1),main = "histogram of h2s")
        }

        k = nrow(Posterior$Lambda)/p;
        h2s = Posterior$F_h2[,1:sp_num]
        G_Lambdas = array(0,dim = dim(Posterior$Lambda))
        Lambda_est = matrix(0,p,k)
        G_est = E_est = matrix(0,p,p)
        for(j in 1:sp_num) {
            Lj = matrix(Posterior$Lambda[,j],p,k);
            h2j = Posterior$F_h2[,j];
            G_Lj = Lj %*%  diag(sqrt(h2j));
            G_Lambdas[,j] = c(G_Lj)
            Gj = G_Lj %*%  t(G_Lj) + diag(1/Posterior$E_a_prec[,j])
            G_est = G_est + Gj/sp_num
            
            E_Lj = Lj  %*% diag(1-h2j) %*%  t(Lj) + diag(1/Posterior$resid_Y_prec[,j])
            E_est = E_est + E_Lj/sp_num;
            Lambda_est = Lambda_est + matrix(Posterior$Lambda[,j],p,k)/sp_num;
        }
        G_Lambda = matrix(rowMeans(G_Lambdas),p,k)

        dev.set(devices[4])
        par(mfrow=c(3,3))
        cors = abs(cor(G_Lambda,actual_E_Lambda))
    	cors=rbind(cors,apply(cors,2,max))
    	cors = cbind(cors,apply(cors,1,max))
    	image(cors[,ncol(cors):1],zlim=c(0,1),main = "cor of G_lambda vs act_E_lambda")

    	image(CovToCor(G_est),zlim=c(-1,1),main = "G_est")
    	image(CovToCor(E_est),zlim=c(-1,1),main = "E_est")

    	clim=max(0.1,min(.75,max(max(abs(G_Lambda)))))
    	clims=c(-clim,clim)
    	image(t(G_Lambda),zlim=clims,main = "G_lambda")

    	image(t(actual_G_Lambda),zlim=clims,main = "act_G_lambda")

    	plot(c(G_est),c(G_act),main = "G_est vs G_act");abline(0,1)
        plot(c(E_est),c(E_act),main = "E_est vs E_act");abline(0,1)

        plot(diag(G_est)/(diag(G_est)+diag(E_est)),h2s_act,xlim=c(0,1),ylim=c(0,1),main = "h2s_est vs h2s_act")
        abline(0,1)

    	F_h2 = rowMeans(Posterior$F_h2[,1:sp_num])
    	E_h2 = 1-F_h2
    	plot(F_h2,type='l',ylim=c(0,1),main = "F_h2 vs E_h2")
    	lines(E_h2,col=2)
    }
}

draw_results_diagnostics = function(sp_num,params,run_variables,Lambda, F_h2, Posterior,E_a_prec,resid_Y_prec,traitnames){
    # devices = dev.list()
    # while(length(devices) < 7){
    #     if(.Platform$OS.type != "windows") {
    #        quartz()
    #     } else {
    #         windows()
    #     }
    #     devices = dev.list()
    # }

    p = run_variables$p
    E_h2 = 1-F_h2;
    G_Lambda = sweep(Lambda,2,sqrt(F_h2),'*')
    E_Lambda = sweep(Lambda,2,sqrt(E_h2),'*')
    G_est = G_Lambda %*% t(G_Lambda) + diag(1/E_a_prec)
    E_est = E_Lambda %*% t(E_Lambda) + diag(1/resid_Y_prec)
 
    # Figure 1: accuracy of current parameter estimates
    #dev.set(devices[1])
    pdf('plotting_diagnostics.pdf')
    par(mfrow=c(2,2))
    # plot 1: estimated number of important factors   
    factor_variances = diag(t(Lambda) %*% Lambda)
    factor_variances = sort(factor_variances,decreasing=T)
    plot(cumsum(factor_variances)/sum(factor_variances),type='l',main = "cumsum importance of factor")

    # plots 2-3: visualize current genetic and residual correlation matrices
    image(CovToCor(G_est),zlim=c(-1,1), main = "G_est")
    image(CovToCor(E_est),zlim=c(-1,1),main = "E_est")

    # plot 4: Factor heritablities
    plot(F_h2,type='l',ylim=c(0,1), main = "F_h2 vs E_h2",col="blue")
    lines(E_h2,col="red")
    legend("topright",legend = c("F_h2","E_h2"),col = c("blue","red"),text.col = c("blue","red"),bty = "n",pch = 1)
    #dev.off()
    
    if(sp_num>1){

        #dev.set(devices[2])
        # Figure of trace plots of the largest elements of each column of the factor loading matrix
        #png('trace plots column of the factor loading matrix.png',width = )
        f2_row=4;
        f2_col=4;
        par(mfrow=c(f2_row,f2_col))
        
        # for(k in 0:(min(2*f2_row,nrow(Posterior$Lambda)/p)-1)) {
        for(k in 0:(f2_row*f2_col-1)){
            if(k >= nrow(Posterior$Lambda)/p) next
            o = order(-abs(rowMeans(Posterior$Lambda[(1:p)+k*p,max(1,sp_num-100):sp_num])))
            traces = Posterior$Lambda[o[1:5]+k*p,1:sp_num]
            trace_plot(traces,main = sprintf('l_%d,.',k),ylim = c(-1,1)*max(abs(traces)))
            abline(h=0)            
        }
        #dev.off()
        
        #dev.set(devices[3])
        # Figure of trace plots and histograms of the factor heritablities
        f4_row=4;
        f4_col=4;
        #pdf('trace plots and histograms of the factor heritablities.pdf')
        par(mfrow=c(f4_row,f4_col))
                
        for(k in 1:min(2*f4_row,nrow(Posterior$F_h2))) {
            h2s = Posterior$F_h2[k,1:sp_num]
            if(sum(h2s[!is.na(h2s)])==0) {
                next
            }
            plot(h2s,type='l', main = "plot of h2s")
            hist(h2s,breaks=100,xlim=c(-0.1,1),main = "histogram of h2s")
        }

        k = nrow(Posterior$Lambda)/p;
        h2s = Posterior$F_h2[,1:sp_num]
        G_Lambdas = array(0,dim = dim(Posterior$Lambda))
        Lambda_est = matrix(0,p,k)
        G_est = E_est = matrix(0,p,p)
        traces_G = matrix(,p*(p+1)/2,sp_num)
        traces_G_cor = matrix(,p*(p+1)/2,sp_num)
        traces_E = matrix(,p*(p+1)/2,sp_num)
   
        for(j in 1:sp_num) {
            Lj = matrix(Posterior$Lambda[,j],p,k)
            h2j = Posterior$F_h2[,j]
            G_Lj = Lj %*%  diag(sqrt(h2j))
            G_Lambdas[,j] = c(G_Lj)
            Gj = G_Lj %*%  t(G_Lj) + diag(1/Posterior$E_a_prec[,j])
            G_est = G_est + Gj/sp_num
            library(gdata)
            traces_G[,j] = lowerTriangle(Gj,diag = TRUE)
            traces_G_cor[,j] = lowerTriangle(CovToCor(Gj),diag = TRUE)
            
            E_Lj = Lj  %*% diag(1-h2j) %*%  t(Lj) + diag(1/Posterior$resid_Y_prec[,j])
            E_est = E_est + E_Lj/sp_num;
            Lambda_est = Lambda_est + matrix(Posterior$Lambda[,j],p,k)/sp_num;
            traces_E[,j] = lowerTriangle(E_Lj,diag = TRUE)
        }
        G_Lambda = matrix(rowMeans(G_Lambdas),p,k)
        #dev.off()
        save(G_est,file = "G_est.RData")
        # Figure of posterior mean estimates
        #dev.set(devices[4])
        #pdf('Figure of posterior mean estimates.pdf')
        par(mfrow=c(2,2))
        #plot 1-2: visualize posterior mean genetic and residual correlation matrices
        image(CovToCor(G_est),zlim=c(-1,1),main = "G_est")
        image(CovToCor(E_est),zlim=c(-1,1),main = "E_est")

        # plot 3: Lambda matrix
        clim=max(0.1,min(.75,max(max(abs(Lambda)))))
        clims=c(-clim,clim)
        image(t(Lambda),zlim=clims, main="transpose of lambda")

        # plot 4: posterior mean factor heritabilities
        F_h2 = rowMeans(Posterior$F_h2[,1:sp_num])
        E_h2 = 1-F_h2
        plot(F_h2,type='l',ylim=c(0,1),main = "F_h2 vs E_h2",col="blue")
        lines(E_h2,col="red")
        legend("topright",legend = c("F_h2","E_h2"),col = c("blue","red"),text.col = c("blue","red"),bty = "n",pch = 1)
        #dev.off()
        
        #plot 5: trace plot of G_est
        #dev.set(devices[5])
        #pdf('trace plot of G_est.pdf')
        par(mfrow = c(2,1))
        trace_plot(traces_G,main = "trace plot of Gj")
        trace_plot(traces_G_cor,main = "trace plot of Gj_cor")
        #dev.off()
        
        #plot 6: distance matrix plot 
        #dev.set(devices[6])
        #pdf('distance matrix plot.pdf')
        par(mfrow = c(1,2))
        Gcov2 = G_est^2
        rownames(Gcov2) = traitnames
        plot(hclust(dist(1-Gcov2)), main = "distance matrix plot of Gcov2")
        
        Gcor = CovToCor(G_est) 
        Gcor2 = Gcor^2
        rownames(Gcor2) = traitnames
        plot(hclust(dist(1-Gcor2)),main = "distance matrix plot of Gcor2")
        #dev.off()
      
        #plot 7: distribution of lambda for each trait 
        #pdf('Lambda_columns_posterior.pdf')
        lam_col = 2
        lam_row = 2
        par(mfrow = c(lam_row,lam_col))
        Lambda = Posterior$Lambda
        
        for(i in 1:k) {
          par(mar=c(10,3,3,1))
          Li_samples = t(Lambda[(i-1)*p + 1:p,])
          colnames(Li_samples) = BSFG_state$traitnames
          boxplot(Li_samples,ylim = c(-1,1)*max(abs(Li_samples)),
                    main=sprintf("boxplot of Lambda_%d by traits",i),
                  las = 2,cex.axis=.6)
          abline(h=0)
        }
        #dev.off()
        
        #plot 8: distribution of upper triangle Gcor 
        #dev.set(devices[7])
        #pdf('boxplot of G_cor.pdf')
        par(mfrow = c(1,2))
        med = apply(t(traces_G_cor), 2, function(x)quantile(x,.5))
        L = apply(t(traces_G_cor), 2, function(x)quantile(x,.05))
        U = apply(t(traces_G_cor), 2, function(x)quantile(x,.95))
        require(plotrix)
        
        plotCI(1:(p*(p+1)/2),med,ui=U,li=L,main = "90% C.I of G_cor")
        boxplot(t(traces_G_cor),main="boxplot of G_cor")
        dev.off()
  }
}



plot_diagnostics = function(BSFG_state){
  burn           = BSFG_state$run_parameters$burn
  thin           = BSFG_state$run_parameters$thin
  start_i        = BSFG_state$current_state$nrun
  params         = BSFG_state$params
  run_variables  = BSFG_state$run_variables
  Lambda         = BSFG_state$current_state$Lambda
  F_h2           = BSFG_state$current_state$F_h2
  Posterior      = BSFG_state$Posterior
  E_a_prec       = BSFG_state$current_state$E_a_prec
  resid_Y_prec   = BSFG_state$current_state$resid_Y_prec
  traitnames     = BSFG_state$traitnames
  draw_iter      = BSFG_state$run_parameters$draw_iter
#   for(i in start_i+(1:n_samples)){
#   if( (i-burn) %% thin == 0 && i > burn) {
#       sp_num = (i-burn)/thin    
#   }
#     if(i %% draw_iter  == 0) {
  sp_num = ncol(BSFG_state$Posterior$Lambda)    
    draw_results_diagnostics(sp_num,params,run_variables,Lambda, F_h2, Posterior,E_a_prec,resid_Y_prec,traitnames)
    
  # }
}
# load("BSFG_state.RData")
# plot_diagnostics(BSFG_state)

ComparingGMatrix_plot = function(target){
  #load data from the original population
  load("BSFG_state.RData")
  spn = dim(BSFG_state$Posterior[[target]])[2]
  n   = dim(BSFG_state$data_matrices$Y)[1]
  pos = BSFG_state$Posterior[[target]][,spn]
  if (target!="F_h2"){
    pos = matrix(pos,nr=n)
    k   = nrow(BSFG_state$Posterior[[target]])/n
  }
  
  #load data from new population
  load("BSFG_fixedlambda.RData")
  n   = dim(BSFG_state$data_matrices$Y)[1]
  star = BSFG_state$Posterior[[target]][,spn]
  if (target!="F_h2"){
    star = matrix(star,nr=n)
    pdf(sprintf("comparing_%s_densityplot.pdf",target))
    for(i in 1:k){
      plot(density(pos[,i]),main = sprintf("%s %d",target,i),col = "blue",type = "l",xlab = "#obs")
      # plot(density(F_a_pos[,i]),main = sprintf("%d",i),col = "blue",type = "l",ylim = c(min(F_a_pos)-5,max(F_a_pos)+5),xlab = "#obs")
      lines(density(star[,i]),col="red",type = "l")
      abline
      legend("topright",legend = c("original","new"),col = c("blue","red"),text.col = c("blue","red"),bty = "n",pch = 1)
    }
    dev.off()
  }else{
    pdf(sprintf("comparing_%s_densityplot.pdf",target))
    plot(density(pos),main = sprintf("%s",target),col = "blue",type = "l",xlab = "#obs")
    # plot(density(F_a_pos[,i]),main = sprintf("%d",i),col = "blue",type = "l",ylim = c(min(F_a_pos)-5,max(F_a_pos)+5),xlab = "#obs")
    lines(density(star),col="red",type = "l")
    abline
    legend("topright",legend = c("original","new"),col = c("blue","red"),text.col = c("blue","red"),bty = "n",pch = 1)
    dev.off()
  }
}


