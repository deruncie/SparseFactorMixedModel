CovToCor = function(X){
	# Normalize a covariance matrix into a correlation matrix
    corX=diag(1/sqrt(diag(X))) %*% X %*% diag(1/sqrt(diag(X)));
    return(corX)
}

trace_plot = function(data,main = NULL){
	plot(NA,NA,xlim = c(0,ncol(data)),ylim = range(data),main= main)
	for(i in 1:nrow(data)) lines(data[i,],col=i)
}

draw_simulation_diagnostics = function(sp_num,run_parameters,run_variables,Posterior,Lambda,F_h2,E_a_prec,resid_Y_prec){

    devices = dev.list()
    while(length(devices) < 4){
    	quartz()
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
    image(cors[,ncol(cors):1],zlim=c(0,1))

    image(CovToCor(G_est),zlim=c(-1,1))
    image(CovToCor(E_est),zlim=c(-1,1))

    clim=max(0.1,min(.75,max(max(abs(G_Lambda)))))
    clims=c(-clim,clim)
    image(t(G_Lambda),zlim=clims)

    image(t(actual_G_Lambda),zlim=clims)

    clim=max(0.1,min(.75,max(max(abs(E_Lambda)))))
    clims=c(-clim,clim)
    image(t(E_Lambda),zlim=clims)

    plot(c(G_est),c(G_act));abline(0,1)

    plot(c(E_est),c(E_act));abline(0,1)

    plot(F_h2,type='l',ylim=c(0,1))
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
            plot(h2s,type='l')
            hist(h2s,breaks=100,xlim=c(-0.1,1))
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
    	image(cors[,ncol(cors):1],zlim=c(0,1))

    	image(CovToCor(G_est),zlim=c(-1,1))
    	image(CovToCor(E_est),zlim=c(-1,1))

    	clim=max(0.1,min(.75,max(max(abs(G_Lambda)))))
    	clims=c(-clim,clim)
    	image(t(G_Lambda),zlim=clims)

    	image(t(actual_G_Lambda),zlim=clims)

    	plot(c(G_est),c(G_act));abline(0,1)
        plot(c(E_est),c(E_act));abline(0,1)

        plot(diag(G_est)/(diag(G_est)+diag(E_est)),h2s_act,xlim=c(0,1),ylim=c(0,1))
        abline(0,1)

    	F_h2 = rowMeans(Posterior$F_h2[,1:sp_num])
    	E_h2 = 1-F_h2
    	plot(F_h2,type='l',ylim=c(0,1))
    	lines(E_h2,col=2)
    }
}

draw_results_diagnostics = function(sp_num,params,Lambda, F_h2, Posterior){
    devices = dev.list()
    while(length(devices) < 4){
        quartz()
        devices = dev.list()
    }

    p = run_variables$p
    E_h2 = 1-F_h2;
    G_Lambda = sweep(Lambda,2,sqrt(F_h2),'*')
    E_Lambda = sweep(Lambda,2,sqrt(E_h2),'*')
    G_est = G_Lambda %*% t(G_Lambda) + diag(1/E_a_prec)
    E_est = E_Lambda %*% t(E_Lambda) + diag(1/resid_Y_prec)
 
    # Figure 1: accuracy of current parameter estimates
    dev.set(devices[1])
    par(mfrow=c(2,2))
    # plot 1: estimated number of important factors   
    factor_variances = diag(t(Lambda) %*% Lambda)
    factor_variances = sort(factor_variances,decreasing=T)
    plot(cumsum(factor_variances)/sum(factor_variances),type='l')

    # plots 2-3: visualize current genetic and residual correlation matrices
    image(CovToCor(G_est),zlim=c(-1,1))
    image(CovToCor(E_est),zlim=c(-1,1))

    # plot 4: Factor heritablities
    plot(F_h2,type='l',ylim=c(0,1))
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
            plot(h2s,type='l')
            hist(h2s,breaks=100,xlim=c(-0.1,1))
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


        # Figure of posterior mean estimates
        dev.set(devices[4])
        par(mfrow=c(2,2))
        #plot 1-2: visualize posterior mean genetic and residual correlation matrices

        image(CovToCor(G_est),zlim=c(-1,1))
        image(CovToCor(E_est),zlim=c(-1,1))

        # plot 3: Lambda matrix
        clim=max(0.1,min(.75,max(max(abs(Lambda)))))
        clims=c(-clim,clim)
        image(t(Lambda),zlim=clims)

        # plot 4: posterior mean factor heritabilities
        F_h2 = rowMeans(Posterior$F_h2[,1:sp_num])
        E_h2 = 1-F_h2
        plot(F_h2,type='l',ylim=c(0,1))
        lines(E_h2,col=2)
    }
}