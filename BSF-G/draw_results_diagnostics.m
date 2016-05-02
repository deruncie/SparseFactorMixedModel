function draw_results_diagnostics(sp_num,params,Lambda, F_h2, Posterior)

p = params.p;
fig1=figure(1);
f1_row=3;
f1_col=2;
clf
p_view = 1:p;

E_h2 = 1-F_h2;
G_Lambda = bsxfun(@times,Lambda,sqrt(F_h2'));
E_Lambda = bsxfun(@times,Lambda,sqrt(E_h2'));

%Figure 1: accuracy of current parameter estimates
%plot 1: estimated number of important factors
subplot(f1_row,f1_col,1)
factor_variances = diag(Lambda'*Lambda);
factor_variances = sort(factor_variances,'descend');
plot(cumsum(factor_variances)/sum(factor_variances));

%plot 3: visualize estimated genetic factor loading matrix
subplot(f1_row,f1_col,3)
clim=max(.1,min(.75,max(max(abs(G_Lambda)))));
clims=[-clim clim];
imagesc(G_Lambda(p_view,:),clims)
colorbar;

%plot 5: visualize estimated residual factor loading matrix
subplot(f1_row,f1_col,5)
clim=max(.1,min(.75,max(max(abs(E_Lambda)))));
clims=[-clim clim];
imagesc(E_Lambda(p_view,:),clims)
colorbar;

%plot 6: Plot of factor heritabilities
subplot(f1_row,f1_col,6)
plot(F_h2,'Color','r')
hold on
plot(E_h2,'Color','b')
hold off
ylim([0 1])


if sp_num>10,

    %Figure of trace plots of the largest elements of each column of the factor loading matrix
    fig2 = figure(2);
    set(fig2,'Position',[0 800,600,800]);
    clf;
    f2_row=8;
    f2_col=3;
    for k=0:(min(f2_col*f2_row,size(Posterior.Lambda,1)/p)-1)
        subplot(f2_row,f2_col,k+1)
        [~,o] = sort(-abs(mean(Posterior.Lambda((1:p)+k*p,max(1,sp_num-100):sp_num),2)));
        plot(Posterior.Lambda(o(1:5)+k*p,1:sp_num)')
    end

    %Figure of trace plots and histograms of the factor heritablities
    fig4 = figure(4);
    set(fig4,'Position',[800 800,600,800]);
    clf
    f4_row=8;
    f4_col=3;
    for k=1:min(f4_col*f4_row,size(Posterior.F_h2,1))
        h2s = Posterior.F_h2(k,1:sp_num);
        if sum(h2s(~isnan(h2s)))==0
            continue
        end
        subplot(f4_row,2*f4_col,2*k-1)
        plot(h2s)
%         ylim([0 max(sds(~isnan(sds)))])
        ylim([0 1])
        subplot(f4_row,2*f4_col,2*k)
        hist(h2s,100)
        xlim([0 1])
    end               

    figure(fig1);  
    
    k = size(Posterior.Lambda,1)/p;
    h2s = Posterior.F_h2(:,1:sp_num);     
    G_Lambdas = zeros(size(Posterior.Lambda));
    G_est = zeros(p,p);
    E_est = zeros(p,p);
    for j=1:sp_num
        Lj = reshape(Posterior.Lambda(:,j),p,k);
        h2j = Posterior.F_h2(:,j);
        G_Lj = Lj * diag(sqrt(h2j));
        G_Lambdas(:,j) = G_Lj(:);
        Gj = G_Lj * G_Lj' + diag(1./Posterior.E_a_prec(:,j));
        G_est = G_est + Gj./sp_num;
        
        E_Lj = Lj * diag(1-h2j)* Lj' + diag(1./Posterior.resid_Y_prec(:,j));
        E_est = E_est + E_Lj./sp_num;
    end
    G_Lambda = reshape(mean(G_Lambdas,2),p,k);
    
    %plot 2: visualize posterior mean genetic correlation matrix
    subplot(f1_row,f1_col,2)
    clims=[-1 1];
    imagesc(CovToCor(G_est),clims)
 
    %plot 4: visualize posterior mean genetic factor loading matrix
    subplot(f1_row,f1_col,4)
    clim=min(.75,max(max(abs(G_Lambda))));
    clims=[-clim clim];
    imagesc(G_Lambda(p_view,:),clims)
    lx = xlim;
    colorbar;  
    
    %plot 6: posterior mean factor heritabilities
    subplot(f1_row,f1_col,6)
    plot(mean(h2s,2),'Color','r');
    ylim([0 1])
    hold on
    plot(1-mean(h2s,2),'Color','b');
    hold off
    colorbar
    xlim(lx)

     drawnow()
    saveas(fig2,'Diagnostics_autoCorr_Lambda.jpg','jpg')
    saveas(fig4,'Diagnostics_factor_selection.jpg','jpg')
end
drawnow()
saveas(fig1,'Diagnostics_Lambda.jpg','jpg')
        
end

