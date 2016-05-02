%% prepare output directory
clear all

a=clock;s=1e6*mod(a(6),1)   %reset random seeds
rand('seed',s)
randn('seed',s)

clear Posterior %clear previous runs

rep = 1 %can be set by shell at runtime for repeated runs.
folder=strcat('rep',num2str(rep))

mkdir(folder);
cd(folder);

%if needed, add path to model scripts and functions
%addpath(genpath('../../scripts'))

%% Initialize run parameters and prior hyper-parameters
% run parameters 
draw_iter=100;
burn=1000;sp=1000;thin=10;
b0=1;b1=0.0005;
epsilon=1e-2;
h2_divisions = 100;

% prior hyperparamters

k_init = 20; % initial number of factors to initialize
as=2; bs=1/10; % inverse gamma hyperparameters for model residuals, as well as non-factor random effects
df = 3; % degrees of freedom for t-distribution of ARD prior on factor loadings
ad1 = 2.1; bd1 = 1/20; % inverse gamma hyperparamters for first factor shrinkage multiplier (/delta_1)
ad2 = 3; bd2 = 1; % % inverse gamma hyperparamters for remaining factor shrinkage multiplier (/delta_i, i \in 2...k)

save('priors','as','bs', ...
'k_init','df','ad1','bd1','ad2','bd2')
priors = load('priors'); % this is done so that priors can be set in the shell at runtime, with the above commented out.


%% Run Gibbs sampler
c1=fix(clock)
[Posterior, params] = fast_BSF_G_sampler(burn,sp,thin,b0,b1,h2_divisions,epsilon,priors,draw_iter);
c2=fix(clock)
c2-c1

%% Calculate posterior means of G, R, P, Lambda, factor_h2s, h2s
p = size(Posterior.ps,1); sp_num = size(Posterior.Lambda,2);
VY = ones(p,1)';
k = size(Posterior.Lambda,1)/p;
h2s = Posterior.G_h2(:,1:sp_num);     
G_Lambdas = zeros(size(Posterior.Lambda));

Lambda_est = zeros(p,k);
G_est = zeros(p,p);
P_est = zeros(p,p);
factor_h2s_est = zeros(k,1);
for j=1:sp_num
    Lj = bsxfun(@times,reshape(Posterior.Lambda(:,j),p,k),1./sqrt(VY'));
    h2j = Posterior.G_h2(:,j);
    
    Pj = Lj*Lj' + diag(1./(VY'.*Posterior.ps(:,j)))+ diag(1./(VY'.*Posterior.resid_ps(:,j)));
    Gj = Lj*diag(h2j)*Lj'  + diag(1./(VY'.*Posterior.ps(:,j)));
    
    Lambda_est = Lambda_est + Lj./sp_num;
    P_est = P_est  + Pj./sp_num;
    G_est = G_est + Gj./sp_num;
    
    factor_h2s_est = factor_h2s_est + h2j./sp_num;
end

cors = abs(Cor(Lambda_est',params.Lambda'));

k=size(params.Lambda,2);
o = zeros(k,1);
for j=1:k,
    o(j) = find(cors(:,j) == max(cors(:,j)));
end
% [~,o]=sort(max_cors);
Lambda_est_ordered = Lambda_est(:,o);
factor_h2s_est = factor_h2s_est(o);

factor_angles = zeros(k,1);
for j=1:k,
    a = Lambda_est_ordered(:,j);
    b = params.Lambda(:,j);
    costheta = dot(a,b)/(norm(a)*norm(b));
    theta = acos(costheta);
    factor_angles(j) = pi/2-abs(pi/2-theta);
end

%    Figure of trace plots and histograms of the factor heritablities
fig4 = figure(4);
set(fig4,'Position',[0 800,600,800]);
clf
f4_row=8;
f4_col=4;
for ji=1:min(4*8/2,length(o))
    j = o(ji);
    h2s = Posterior.G_h2(j,:);
    if sum(h2s(~isnan(h2s)))==0
        continue
    end
    subplot(f4_row,f4_col,2*ji-1)
    plot(h2s)
    ylim([0 1])
    subplot(f4_row,f4_col,2*ji)
    hist(h2s,100)
    xlim([-0.1 1])
end        

 
posterior_mean.G = G_est;
posterior_mean.P = P_est;
posterior_mean.Lambda = Lambda_est;
posterior_mean.factor_order = o;
posterior_mean.factor_h2s_est = factor_h2s_est;
posterior_mean.factor_angles = factor_angles;

save('Posterior_mean','posterior_mean')

%% Produce plot of results
figure(5)
p = size(Posterior.ps,1); sp_num = size(Posterior.Lambda,2);
VY = ones(p,1)';
k = size(Posterior.Lambda,1)/p;
h2s = Posterior.G_h2(:,1:sp_num);     
G_Lambdas = zeros(size(Posterior.Lambda));
G_est = zeros(p,p);
R_est = zeros(p,p);
for j=1:sp_num
    Lj = bsxfun(@times,reshape(Posterior.Lambda(:,j),p,k),1./sqrt(VY'));
    h2j = Posterior.G_h2(:,j);
    G_Lj = Lj * diag(sqrt(h2j));
    G_Lambdas(:,j) = G_Lj(:);
    Gj = G_Lj * G_Lj' + diag(1./(VY'.*Posterior.ps(:,j)));%
    G_est = G_est + Gj./sp_num;

    E_Lj = Lj * diag(1-h2j)* Lj' + diag(1./(VY'.*Posterior.resid_ps(:,j)));
    R_est = R_est + E_Lj./sp_num;
end
G_Lambda = reshape(mean(G_Lambdas,2),p,k);

f1_row = 1; f1_col=2;
subplot(f1_row,f1_col,1)
clims=[-1 1];
imagesc(CovToCor(G_est),clims)
axis square

subplot(f1_row,f1_col,2)
clim=min(.75,max(max(abs(G_Lambda))));
clims=[-clim clim];
imagesc(G_Lambda,clims)
lx = xlim;
colorbar;  

% export_fig('transparent','append','nocrop','/Users/der7/Documents/Statistics/Sparse_factor_G_matrix/Genetics_paper/BGSFM/Instructions/Sim.pdf')

cd('..')

%%
exit
