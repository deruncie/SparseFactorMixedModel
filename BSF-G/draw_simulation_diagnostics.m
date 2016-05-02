function draw_simulation_diagnostics(sp_num,params,Posterior,Lambda,F_h2,E_a_prec,resid_Y_prec,gen_factor_Lambda,error_factor_Lambda,G_act,E_act,h2)
% Factors,genetic_effects,resid,...
%         Posterior,gen_factor_Lambda,error_factor_Lambda,G,R,h2)
%draw some diagnostic plots    
    
p = params.p;
E_h2 = 1-F_h2;
G_Lambda = bsxfun(@times,Lambda,sqrt(F_h2'));
E_Lambda = bsxfun(@times,Lambda,sqrt(E_h2'));
actual_G_Lambda = gen_factor_Lambda;
actual_E_Lambda = error_factor_Lambda;

G_est = G_Lambda*G_Lambda' + diag(1./E_a_prec);

E_est = E_Lambda*E_Lambda' + diag(1./resid_Y_prec);
 

fig1=figure(1);
f1_row=3;
f1_col=3;
clf
p_view = 1:p;

%Figure 1: accuracy of current parameter estimates
%plot 1: correlation between each estimated factor and the true factors.
%Outer row and column shows maximum correlation for each fitted or actual
%factor, respectively.
subplot(f1_row,f1_col,1)
clims=[0 1];        
cors = abs(Cor(Lambda',[actual_E_Lambda]'));
cors=[cors;max(cors)];
cors=[cors max(cors')'];
imagesc(cors,clims)
colorbar

if sp_num <= 1,
    k=min(size(actual_E_Lambda,2),size(Lambda,2));
    max_cors = zeros(k,1);
    cors = cors(1:end-1,1:end-1);
    for j=1:k,
        max_cors(j) = find(cors(j,:) == max(cors(j,:)));
    end
    [~,o]=sort(max_cors);
    L=Lambda(:,o);

    figure(6)
    for j=1:min(k,5*5),
        subplot(5,5,j)
        plot(actual_E_Lambda(:,j),L(:,j),'.')
        axis equal
        line(xlim,xlim)
        line(-xlim,xlim)
    end
end



figure(fig1);

%plot 2: visualize estimated genetic covariance matrix (as correlation matrix)
subplot(f1_row,f1_col,2)
clims=[-1 1];
imagesc(CovToCor(G_est),clims)

%plot 3: visualize estimated residual covariance matrix (as correlation matrix)
subplot(f1_row,f1_col,3)
clims=[-1 1];
imagesc(CovToCor(E_est),clims)

%plot 4: visualize estimated genetic factor loading matrix
subplot(f1_row,f1_col,4)
clim=max(0.1,min(.75,max(max(abs(G_Lambda)))));
clims=[-clim clim];
imagesc(G_Lambda(p_view,:),clims)
colorbar;

%plot 5: visualize actual genetic factor loading matrix
subplot(f1_row,f1_col,5)
imagesc([actual_G_Lambda(p_view,:) zeros(length(p_view),1)],clims)
colorbar;

%plot 6: visualize estimated residual factor loading matrix
subplot(f1_row,f1_col,6)
clim=max(0.1,min(.75,max(max(abs(E_Lambda)))));
clims=[-clim clim];
imagesc(E_Lambda(p_view,:),clims)
colorbar;

%plot 7: Element-wise comparison of estimated and actual G matrices
subplot(f1_row,f1_col,7)
plot(G_est,G_act,'.')
line(xlim,xlim)

%plot 8: Element-wise comparison of estiamted and actual R matrices
subplot(f1_row,f1_col,8)
plot(triu(E_est+G_est,1),triu(E_act+G_act,1),'.')
line(xlim,xlim)

%plot 9: Plot of factor heritabilities
subplot(f1_row,f1_col,9)
plot(F_h2,'Color','r');
ylim([0 1])
hold on
plot(E_h2,'Color','b');
hold off

if sp_num>1,

    %Figure of trace plots of the largest elements of each column of the factor loading matrix
    fig2 = figure(2);
    set(fig2,'Position',[0 800,600,800]);
    clf;
    f2_row=8;
    f2_col=2;
    for k=0:(min(2*f2_row,size(Posterior.Lambda,1)/p)-1)
        subplot(f2_row,f2_col,k+1)
        [~,o] = sort(-abs(mean(Posterior.Lambda((1:p)+k*p,max(1,sp_num-100):sp_num),2)));
        plot(Posterior.Lambda(o(1:5)+k*p,1:sp_num)')
    end

    %Figure of trace plots and histograms of the factor heritablities
    fig4 = figure(4);
    set(fig4,'Position',[0 800,600,800]);
    clf
    f4_row=8;
    f4_col=4;
    for k=1:min(2*f4_row,size(Posterior.F_h2,1))
        h2s = Posterior.F_h2(k,1:sp_num);
        if sum(h2s(~isnan(h2s)))==0
            continue
        end
        subplot(f4_row,f4_col,2*k-1)
        plot(h2s)
        ylim([0 1])
        subplot(f4_row,f4_col,2*k)
        hist(h2s,100)
        xlim([-0.1 1])
    end        

    %Figure similar to figure 1, but posterior mean parameter estimates 
    fig3 = figure(3);
    f1_row=3;
    f1_col=2;
    clf;      

    k = size(Posterior.Lambda,1)/p;
    h2s = Posterior.F_h2(:,1:sp_num);     
    G_Lambdas = zeros(size(Posterior.Lambda));
    Lambda_est = zeros(p,k);
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
        Lambda_est = Lambda_est + reshape(Posterior.Lambda(:,j),p,k)./sp_num;
    end
    G_Lambda = reshape(mean(G_Lambdas,2),p,k);
    
    %plot 1: visualize estimated genetic covariance matrix (as correlation matrix).
    subplot(f1_row,f1_col,1)
    clims=[-1 1];
    imagesc(CovToCor(G_est),clims)
 
    %plot 4: visualize estimated genetic factor loading matrix
    subplot(f1_row,f1_col,4)
    clim=min(.75,max(max(abs(G_Lambda))));
    clims=[-clim clim];
    imagesc(G_Lambda(p_view,:),clims)
    lx = xlim;
    colorbar;  
    
    %plot 2: visualize actual genetic factor loading matrix
    subplot(f1_row,f1_col,2)
    imagesc([actual_G_Lambda(p_view,:) zeros(length(p_view),1)],clims)
    colorbar;
 
    %plot 3: Element-wise comparison of estiamted and actual G matrices
    subplot(f1_row,f1_col,3)
    plot(triu(G_est,1),triu(G_act,1),'.')
    line(xlim,xlim)
 
    %plot 5: Plot of trait heritabilities
    subplot(f1_row,f1_col,5)    
    plot(diag(G_est)./(diag(G_est)+diag(E_est)),h2,'.','Color','r')
    line(xlim,xlim)
    
    %plot 6: Plot of factor heritabilities
    subplot(f1_row,f1_col,6)
    plot(mean(h2s,2),'Color','r');
    ylim([0 1])
    hold on
    plot(1-mean(h2s,2),'Color','b');
    hold off
    colorbar
    xlim(lx)

    clims=[0 1];        
    cors = abs(Cor(Lambda_est',[error_factor_Lambda]'));
    cors=[cors;max(cors)];
    cors=[cors max(cors')'];
    

    %plot estimated factor loadings against actual factor loadings for matched factors
    figure(6)
    k=min(size(actual_E_Lambda,2),size(Lambda_est,2));
    max_cors = zeros(k,1);
    cors = cors(1:end-1,1:end-1);
    for j=1:k,
        max_cors(j) = find(cors(:,j) == max(cors(:,j)));
    end
    [~,o]=sort(max_cors);
    L=error_factor_Lambda(:,o);

    for j=1:min(k,5*5),
        subplot(5,5,j)
        plot(L(:,j),Lambda_est(:,j),'.')
        axis equal
        line(xlim,xlim)
        line(-xlim,xlim)
    end


    drawnow()
    drawnow()
    saveas(fig2,'Diagnostics_autoCorr_Lambda.jpg','jpg')
    saveas(fig3,'Diagnostics_posterior_mean.jpg','jpg')
    saveas(fig4,'Diagnostics_factor_selection.jpg','jpg')
end
drawnow()
saveas(fig1,'Diagnostics_Lambda.jpg','jpg')

end

