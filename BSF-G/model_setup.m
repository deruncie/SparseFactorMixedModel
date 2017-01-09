% Daniel Runcie
%
% Gibbs sampler for genetic covariance estimation based on mixed effects
% model, with missing data
%
% Modified for Mark Blows to add another hierarchical level with
% transcripts nested within probes, and all covariances modeled at the
% transcript level.
%
% Based on:
% 
% Runcie and Mukherjee (2013) Dissecting high-dimensional traits
% with Bayesian sparse factor analysis of genetic covariance matrices.
% GENETICS.
%
% (c) July 30, 2013
%
% code based on original provided by Anirban Bhattacharya
%
        
%% set Random Seed: Needed if your version of MATLAB always starts with the
% same seed for the Random number generator, and you want to re-run the
% chain with different starting values. Depending on your version of
% MATLAB, the code to set seeds may vary. This is absed on MATLAB2013a

% a=clock;s=1e6*mod(a(6),1)   %reset random seeds
% rand('seed',s)
% randn('seed',s)
% s = RandStream('mt19937ar','Seed','shuffle');
s = RandStream('mt19937ar','Seed',1);  % for repeatability
RandStream.setGlobalStream(s);

%% Make a new folder to hold the analyses.

rep = 1 %can be set by shell at runtime for repeated runs.
folder=strcat('rep',num2str(rep))

mkdir(folder);
cd(folder);

%if needed, add path to model scripts and functions
%addpath(genpath('../../scripts'))
%% Set priors and model parameters. See MA_sampler_init() for details.

params.b0=1;params.b1=0.0005;
params.epsilon=1e-2;
params.prop = 1.00;
params.h2_divisions = 100;
params.save_freq = 100;

priors.k_init              =   20;
priors.resid_Y_prec_shape  =   2;
priors.resid_Y_prec_rate   =   1/10;
priors.E_a_prec_shape      =   2;
priors.E_a_prec_rate       =   1/10;
priors.W_prec_shape        =   2;
priors.W_prec_rate         =   1/10;
priors.Lambda_df           =   3;
priors.delta_1_shape       =   2;
priors.delta_1_rate        =   1/20;
priors.delta_2_shape       =   3;
priors.delta_2_rate        =   1;
priors.h2_priors_factors   =   (1+[params.h2_divisions-2;zeros(params.h2_divisions-1,1)])/(-2+2*params.h2_divisions);

save('Priors','priors')
%% Initialize Chain, prep runs
[data_matrices,params,priors,current_state,Posterior,simulation] = fast_BSFG_sampler_init(priors,params);  
[Posterior,params] = clear_Posterior(Posterior,params);
start_i = 0;
burn = 100;
thin=2;
draw_iter = 200;

current_state.randStream    =   RandStream.getGlobalStream;
current_state.rand_state    =   RandStream.getGlobalStream.State;

%% optional: To load from end of previous run, run above code, then run these lines:
load('current_state')
load('Posterior')
load('Priors')
start_i = size(Posterior.Lambda,2);

%% Run Gibbs sampler. Run in smallish chunks. Output can be used to re-start chain where it left off.
c1=fix(clock);
sp = 200;
for i =1:10
    [Posterior,current_state,params,start_i] = fast_BSFG_sampler(data_matrices,...
                                    start_i,draw_iter,burn,sp,thin,simulation, ...
                                    params,priors,Posterior,current_state);
	c2=fix(clock);
    c2-c1 
    %The first run (i=1) includes the specified burnin period. After this,
    %set burn = 0. However, if convergence has not yet been achieved, call
    %[Posterior,params] = clear_Posterior(Posterior,params); 
    %burn length = start_i;
    burn = 0;       
end
%%
n = size(data_matrices.Y_full,1);
figure(11);clf;
i=1;
labels = randperm(7);
males = sum(bsxfun(@times,data_matrices.Z_1,labels),2)-1;
males = data_matrices.X_f(:,3)
boxplot(Posterior.F((i-1)*n+(1:n),:)','colorgroup',males,'boxstyle','filled');
%%
P_Fa = Posterior.F_a;
for i = 1:size(P_Fa,2),
    for(j = 1:size(Posterior.F_h2,1)),
        if Posterior.F_h2(j,i) == 0,
            P_Fa((j-1)*7+(1:7),i)=NaN;
        end
    end
end


%%
r = size(data_matrices.Z_1,2);
figure(12);clf;
i=2;
males = labels;
% males = sum(bsxfun(@times,data_matrices.Z_1,randperm(7)),2)-1
% males = data_matrices.X_f(:,1)
boxplot(P_Fa((i-1)*r+(1:r),:)','colorgroup',males,'boxstyle','filled');line(xlim,[0,0])


%%
i=11;
figure(22)
p = size(data_matrices.Y_full,2);
boxplot(Posterior.Lambda((i-1)*p+(1:50),:)');line(xlim,[0,0])

%%
boxplot(Posterior.F_h2')

%% Export Posterior to CSV
fields = fieldnames(Posterior);
maxI = min(10000,size(Posterior.Lambda,2));
for i = 1:size(fields,1),
    field = fields{i}
    mat = getfield(Posterior,field);
    if size(mat,2) > 1000,
        mat = mat(:,1:maxI);
    end
    csvwrite(strcat('Posterior_',field,'.csv'),mat)
end
save('Params','params')

%% Calculate posterior means of G, R, P, Lambda, factor_h2s, h2s
p = size(Posterior.resid_Y_prec,1); sp_num = size(Posterior.Lambda,2);
k = size(Posterior.Lambda,1)/p;

Lambda_est = zeros(p,k);
G_est = zeros(p,p);
P_est = zeros(p,p);
factor_h2s_est = zeros(k,1);
for j=1:sp_num
    Lj = reshape(Posterior.Lambda(:,j),p,k);
    h2j = Posterior.F_h2(:,j);
    
    Pj = Lj*Lj' + diag(1./Posterior.E_a_prec(:,j))+ diag(1./Posterior.resid_Y_prec(:,j));
    Gj = Lj*diag(h2j)*Lj'  + diag(1./Posterior.E_a_prec(:,j));
    
    Lambda_est = Lambda_est + Lj./sp_num;
    P_est = P_est  + Pj./sp_num;
    G_est = G_est + Gj./sp_num;
    
    factor_h2s_est = factor_h2s_est + h2j./sp_num;
end

o = 1:p;
if simulation
    k=size(params.Lambda,2);
    cors = abs(Cor(Lambda_est',params.Lambda'));
    o = zeros(k,1);
    for j=1:k,
        o(j) = find(cors(:,j) == max(cors(:,j)));
    end
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
    posterior_mean.factor_order = o;
    posterior_mean.factor_angles = factor_angles;
end

%    Figure of trace plots and histograms of the factor heritablities
fig4 = figure(4);
set(fig4,'Position',[0 800,600,800]);
clf
f4_row=8;
f4_col=4;
for ji=1:min(4*8/2,length(o))
    j = o(ji);
    F_h2s = Posterior.F_h2(j,:);
    if sum(F_h2s(~isnan(F_h2s)))==0
        continue
    end
    subplot(f4_row,f4_col,2*ji-1)
    plot(F_h2s)
    ylim([0 1])
    subplot(f4_row,f4_col,2*ji)
    hist(F_h2s,100)
    xlim([-0.1 1])
end        

 
posterior_mean.G = G_est;
posterior_mean.P = P_est;
posterior_mean.Lambda = Lambda_est;
posterior_mean.factor_h2s_est = factor_h2s_est;

save('Posterior_mean','posterior_mean')

%% Produce plots of results
figure(5)
p = size(Posterior.resid_Y_prec,1); sp_num = sum(Posterior.Lambda(1,:) ~= 0);p_view = 1:p;
    
k = size(Posterior.Lambda,1)/p;
F_h2s = Posterior.F_h2(:,1:sp_num);     
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


f1_row = 1;f1_col=2;
%plot 1: visualize posterior mean genetic correlation matrix
subplot(f1_row,f1_col,1)
clims=[-1 1];
imagesc(CovToCor(G_est),clims)

%plot 2: visualize posterior mean genetic factor loading matrix
subplot(f1_row,f1_col,2)
clim=min(.75,max(max(abs(G_Lambda))));
clims=[-clim clim];
imagesc(G_Lambda(p_view,:),clims)
lx = xlim;
colorbar;  

% export_fig function (http://www.mathworks.com/matlabcentral/fileexchange/23629-exportfig)
% outputs pdfs that work better than regular MATLAB figure export,
% especially for images.
% export_fig('transparent','append','nocrop','results_plot.pdf')
