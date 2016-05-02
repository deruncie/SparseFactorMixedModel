function [Posterior, params] = fast_BSF_G_sampler(burn,sp,thin,b0,b1,h2_divisions,epsilon,priors,draw_iter)
% -- Daniel Runcie -- %
%
% Gibbs sampler for genetic covariance estimation based on mixed effects
% model, with missing data
%
% Based on:
% 
% Runcie and Mukherjee (2013) Dissecting high-dimensional traits
% with Bayesian sparse factor analysis of genetic covariance matrices.
% GENETICS.
%
% (c) April 22, 2013
%
% code modified from original provided by Anirban Bhattacharya
%
%         
% This function implements the BSF-G partially collapsed Gibbs sampler. It loads all input data
% and matrices from setup.mat in the current directory. Priors and control
% parameters are passed to the function.
%
% setup.mat is a struct with at least:
%     Y: data matrix
%     X: fixed effect design matrix
%     Z_1: random effect design matrix for factor model
%     Z_2: additional random effect design matrix
%     A: Additive genetic relationship matrix
%     
%    For analysis of Ayroles et al 2009 data, can include:
%     Ayroles_results: struct holding gene names and correlations estimated in that paper
%     
%    For analysis of simulations:
%     U_act: r x p matrix of true breeding values
%     E_act: n x p matrix of true model residuals
%     gen_factor_Lambda: p x k_G matrix of true genetic factor loadings
%     error_factor_Lambda: p x k matrix of true residual factor loadings
%     h2: p x 1 vector of true heritabilities
%     factor_h2s: k x 1 vector of true latent factor heritabilities
%     G, R: p x p matrix of true genetic and residual covariances
% 
% The function takes the following inputs:
%     burn: number of burnin samples
%     sp: total number of samples to collect
%     thin: thinning rate of chain
%     b0,b1: parameters controlling rate of adaptation of factor model size
%     h2_divisions: number of discrete steps for each factor heritability parameter
%     epsilon: truncation point for factor loadings during adaptation
%     draw_iter: frequency of updating diagnostic plots
%     priors: struct holding various prior hyperparameters: 
%         k_init: initial number of factors to initialize
%         as, bs: inverse gamma hyperparameters for model residuals, as well as non-factor random effects
%         df: degrees of freedom for t-distribution of ARD prior on factor loadings
%         ad1, bd1: inverse gamma hyperparamters for first factor shrinkage multiplier (/delta_1)
%         ad2, bd2: inverse gamma hyperparamters for remaining factor shrinkage multiplier (/delta_i, i \in 2...k)
% 
% The function output is the struct, Posterior, with the following fields:
%     Lambda, U, no_f, ps, resid_ps, delta, G_h2: matrices with each column a 
%         (vectorized if necessary) posterior sample of the:
%             Lambda: factor loading matrix
%             U: genetic random effect matrix
%             no_f: number of significant factors
%             ps: genetic residual precision after accounting for latent factors
%             resid_ps: phenotypic residual precision
%             delta: column shrinkage parameters for Lambda
%             G_h2: factor heritabilties
%     B, d, W: matrices of posterior means of fixed effect (B), 
%         residual genetic (d), or 2nd random effect (W) coefficients. Some may be empty
% 
% Several diagnostic plots are produced during the run. 
%     Their interpretation is described within the source codes:
%         draw_simulation_diagnostics.m: For simulated data with known true values
%         draw_results_diagnostics.m: Otherwise
%         
    
clear Y
clear Z_1
clear Z_2
clear X

global Y   %n x p matrix of phenotypic data
global Z_1   %n x r incidence matrix for additive genetic effects
global Z_2   %n x r2 incidence matrix for another set of random effects
global X   %n x b design matrix of fixed effects


nrun=burn+sp*thin; % number of posterior samples
k_min = 1e-1; % minimum factor loading size to report in running status
prop = 1.00;  % proportion of redundant elements within columns necessary to drop column

% ------read data--------%
load('../setup.mat');

%Determine if 'setup.mat' contains output of a simulation, based on if
%known factor loadings are included. Affects plotting functions
simulation = true;
if isempty(who('gen_factor_Lambda'))
    simulation=false;
end
simulation

%normalize Y to have zero mean and unit variances among observed values,
%allowing for NaNs.
[n,p]=size(Y);
Y_full=Y;
Mean_Y = zeros(1,size(Y,2));
VY = zeros(1,size(Y,2));
for j=1:p
    Mean_Y(j) = mean(Y(~isnan(Y(:,j)),j));
    VY(j) = var(Y(~isnan(Y(:,j)),j));
    if isnan(VY(j))
        VY(j)=1;
    end
end
Y = bsxfun(@minus,Y,Mean_Y);                     
Y = bsxfun(@times,Y,1./sqrt(VY));       

%determine if a design matrix (X) exists (is loaded from setup.mat). If
%not, make a dummy X-matrix with no columns.
if ~exist('X','var')
    X=zeros(0,n);
end
if size(X,2) ~= n
    X=zeros(0,n);
end

%Determine if a second random effects design matrix exists. If not, make a
%dummy matrix
if ~exist('Z_2','var')
    Z_2 = zeros(0,n);
end
if size(Z_2,2) ~= n
    Z_2 = zeros(0,n);
end
    
%calculate model dimensions
r=size(Z_1,1);
r2=size(Z_2,1);
b=size(X,1);


% --- Initialize variables --- %
%residual parameters. This structure holds the priors hyperparamters for
%the gamma prior on the model residual variances. It also holds the current
%estimate of the residual precision
resid.as=priors.as;
resid.bs=priors.bs;
resid.Y=Y;
resid.p=p;
resid.ps=gamrnd(resid.as,1/resid.bs,p,1);           %residual precision


%Factors. This struct holds all information about the latent factors,
%including the current number of factors, priors hyperparameters for the
%factor Loadings, as well as current values of the factor loadings, their
%precision, the factor scores, and the genetic heritability of each factor
k = priors.k_init;      % initial number of factors
df=priors.df;           % prior degrees of freedom for hierarchical t-prior on loadings 
ad1=priors.ad1;         % priors on delta
bd1=priors.bd1;
ad2=priors.ad2;
bd2=priors.bd2;
Factors.r=n;            
Factors.n=n;
Factors.p=p;
Factors.k=k;
Factors.df=df;
Factors.ad1=ad1;
Factors.bd1=bd1;
Factors.ad2=ad2;
Factors.bd2=bd2;
Factors.psijh = gamrnd(df/2,2/df,[p,k]);    %individual loadings precisions
Factors.delta = [gamrnd(ad1+10,1/bd1);gamrnd(ad2,1/bd2,[k-1,1])];  % components of tauh
Factors.tauh = cumprod(Factors.delta);                              %extra shrinkage of each loading column
Factors.Plam = bsxfun(@times,Factors.psijh,Factors.tauh');          %total precision of each loading
Factors.Lambda = zeros(p,k) + randn(p,k).*reshape(sqrt(1./Factors.Plam),p,k);   %factor loadings
Factors.h2 = rand(k,1);                                             %factor heritability
Factors.h2_divisions = h2_divisions;                                %discretizations of heritability

Factors.num=0;
Factors.no_f=zeros(sp,1);
Factors.nofout = k*ones(nrun,1);

%genetic_effects. This structure holds information about latent genetic
%effects. U is latent genetic effects on factor traits. d is genetic
%effects on residuals of the factor traits. Plus prior hyperparameters for
%genetic effect precisions
genetic_effects.n=size(Z_1,1);
as = priors.as;
bs = priors.bs;
genetic_effects.as=as;
genetic_effects.bs=bs;
genetic_effects.ps = gamrnd(as,1/bs,p,1);
genetic_effects.U =  bsxfun(@times,randn(k,r),sqrt(Factors.h2));
genetic_effects.d = bsxfun(@times,randn(p,r),1./sqrt(genetic_effects.ps));

%interaction_effects. Similar to genetic_effects structure except for
%additional random effects that do not contribute to variation in the
%factor traits
as = priors.as;
bs = priors.bs;
interaction_effects.as=as;
interaction_effects.bs=bs;
interaction_effects.ps = gamrnd(as,1/bs,p,1);
interaction_effects.mean = zeros(p,r2);
interaction_effects.n = r2;
interaction_effects.W = bsxfun(@times,randn(p,r2),1./sqrt(interaction_effects.ps));
interaction_effects.W_out = zeros(p,r2);

%fixed_effects hold B
fixed_effects.b=b;
fixed_effects.cov = zeros(b,b); %inverse covariance of fixed effects
fixed_effects.mean = zeros(p,b); %mean of fixed effects
fixed_effects.B = randn(p,b); %current estimate of fixed effects

Factors.scores = genetic_effects.U*Z_1 + bsxfun(@times,randn(k,n),sqrt(1-Factors.h2)); %initialize factor scores

%Posterior holds Posterior samples and Posterior means
Posterior.Lambda = zeros(0,sp);
Posterior.no_f = zeros(sp,1);
Posterior.ps = zeros(p,sp);
Posterior.resid_ps = zeros(size(resid.ps,1),sp);
Posterior.B = zeros(p,size(X,1));
Posterior.U = zeros(0,sp);
Posterior.d = zeros(p,size(Z_1,1));
Posterior.W = zeros(p,size(Z_2,1));
Posterior.delta = zeros(0,sp);
Posterior.G_h2 = zeros(0,sp);

%save run parameters and hyperparameters

params.p=p;
params.n=n;
params.r=r;
params.Mean_Y=Mean_Y;
params.VY=VY;
params.b0=b0;
params.b1=b1;
params.epsilon=epsilon;
params.prop=prop;
params.as=priors.as;
params.bs=priors.bs;
params.df=priors.df;
params.ad1=priors.ad1;
params.bd1=priors.bd1;
params.ad2=priors.ad2;
params.bd2=priors.bd2;
params.burn=burn;
params.thin=thin;
params.sp=sp;
params.nrun=nrun;
params.h2_divisions = h2_divisions;

if simulation,
    params.U_act = U_act;
    params.Lambda = error_factor_Lambda;
    params.h2 = h2;
    params.G = G;
    params.R = R;
    params.B = B;
    params.factor_h2s = factor_h2s;
    params.name = name;
end

%precalculate some matrices
%invert the random effect covariance matrices
Ainv = inv(A);
A_2_inv = eye(size(Z_2,1)); %Z_2 random effects are assumed to have covariance proportional to the identity. Can be modified.

%pre-calculate transformation parameters to diagonalize aI + bZAZ for fast
%inversion: inv(aI + bZAZ) = 1/b*u*diag(1./(s+a/b))*u'
%uses singular value decomposition of ZAZ for stability when ZAZ is low
%rank
ZAZ=Z_1'*A*Z_1;
[u,s,~]=svd(ZAZ);
eig_ZAZ.vectors = u;
eig_ZAZ.values = diag(s);


%fixed effects + random effects 1
%diagonalize mixed model equations for fast inversion: 
%inv(a*blkdiag(fixed_effects.cov,Ainv) + b*[X; Z_1][X; Z_1]') = Q*diag(1./(a.*s1+b.*s2))*Q'
Design=[X; Z_1];
Design2 = Design*Design';
[~,~,q,S1,S2] = gsvd(cholcov(blkdiag(fixed_effects.cov,Ainv)),cholcov(Design2));
svd_Design_Ainv.Q = inv(q)';
svd_Design_Ainv.s1 = diag(S1'*S1);
svd_Design_Ainv.s2 = diag(S2'*S2);
Qt_Design = svd_Design_Ainv.Q'*Design;      

%random effects 2
%as above, but for random effects 2. Here, fixed effects will be conditioned on, not sampled simultaneously. Otherwise identical.
Design=Z_2;
Design2 = Design*Design';
[~,~,q,S1,S2] = gsvd(cholcov(A_2_inv),cholcov(Design2));
svd_Z2_2_A2inv.Q = inv(q)';
svd_Z2_2_A2inv.s1 = diag(S1'*S1);
svd_Z2_2_A2inv.s2 = diag(S2'*S2);
Qt_Z2 = svd_Z2_2_A2inv.Q'*Design;


%genetic effect variances of factor traits
% diagonalizing a*Z_1*Z_1' + b*Ainv for fast inversion
%diagonalize mixed model equations for fast inversion: 
% inv(a*Z_1*Z_1' + b*Ainv) = Q*diag(1./(a.*s1+b.*s2))*Q'
%similar to fixed effects + random effects 1 above, but no fixed effects.
ZZt = Z_1*Z_1';
[~,~,q,S1,S2] = gsvd(cholcov(ZZt),cholcov(Ainv));
svd_ZZ_Ainv.Q = inv(q)';
svd_ZZ_Ainv.s1 = diag(S1'*S1);
svd_ZZ_Ainv.s2 = diag(S2'*S2);

%------start gibbs sampling-----%
sp_num=0;
tic
for i = 1:nrun
    
   %fill in missing phenotypes
    %conditioning on everything else
    phenMissing = isnan(Y_full);
    if sum(sum(phenMissing))>0
        meanTraits = fixed_effects.B*X +  genetic_effects.d*Z_1 ...
            + interaction_effects.W*Z_2 + Factors.Lambda*Factors.scores;
        meanTraits = meanTraits';        
        resids = bsxfun(@times,randn(size(Y_full)),1./sqrt(resid.ps'));
        Y(phenMissing) = meanTraits(phenMissing) + resids(phenMissing);
    end
    
  %sample Lambda
    %conditioning on W, X, F, marginalizing over D
    Ytil = Y'-fixed_effects.B*X - interaction_effects.W*Z_2;
    [Factors] = sample_lambda( Ytil,Factors, resid,genetic_effects,eig_ZAZ );
    
  %sample fixed effects + random effects 1 ([B;D])
    %conditioning on W, F, L
    Ytil = Y'-interaction_effects.W*Z_2 - Factors.Lambda*Factors.scores;
    N = genetic_effects.n + fixed_effects.b;
    [location_sample] = sample_means( Ytil,Qt_Design,N, ...
            resid, genetic_effects.ps, svd_Design_Ainv );
    fixed_effects.B = location_sample(:,1:fixed_effects.b);
    genetic_effects.d = location_sample(:,fixed_effects.b+1:fixed_effects.b+genetic_effects.n);
        
  %sample random effects 2
    %conditioning on B, D, F, L
    Ytil = Y'-fixed_effects.B*X-genetic_effects.d*Z_1 - Factors.Lambda*Factors.scores;
    N = interaction_effects.n;
    if N>0,
        [location_sample] = sample_means( Ytil,Qt_Z2,N, ...
                resid, interaction_effects.ps, svd_Z2_2_A2inv );
        interaction_effects.W = location_sample;
    end
    
  %sample factor h2
    %conditioning on F, marginalizing over U
    Factors = sample_h2s_discrete(Factors,eig_ZAZ);
    
  %sample genetic effects (U)
    %conditioning on F, Factor h2
    genetic_effects = sample_Us(Factors,genetic_effects,svd_ZZ_Ainv);
           
  %sample F
    %conditioning on U, Lambda, B, D, W, factor h2s
    Ytil = Y'-fixed_effects.B*X - genetic_effects.d*Z_1-interaction_effects.W*Z_2;
    Factors = sample_factors_scores( Ytil, Factors,resid,genetic_effects );

  % -- Update ps -- %
    Lambda2 = Factors.Lambda.^2;    
    Factors.psijh = gamrnd(Factors.df/2 + 0.5,2./(Factors.df + bsxfun(@times,Lambda2,Factors.tauh')));
    
    %continue from previous Y residual above
    Ytil = Ytil - Factors.Lambda*Factors.scores;
    n = size(Y,1);
    resid.ps=gamrnd(resid.as + 0.5*n,1./(resid.bs+0.5*sum(Ytil.^2,2)));  %model residual precision
    n = genetic_effects.n;
    genetic_effects.ps=gamrnd(genetic_effects.as + 0.5*n,1./(genetic_effects.bs+0.5*sum(genetic_effects.d.^2,2))); %random effect 1 (D) residual precision
    n = interaction_effects.n;
    interaction_effects.ps=gamrnd(interaction_effects.as + 0.5*n,1./(interaction_effects.bs+0.5*sum(interaction_effects.W.^2,2))); %random effect 2 (W) residual precision
 
  %------Update delta & tauh------%
    [delta,tauh] = sample_delta( Factors,Lambda2 );
    Factors.delta = delta;
    Factors.tauh = tauh;
    
  %---update precision parameters----%
    Factors.Plam = bsxfun(@times,Factors.psijh,Factors.tauh');

  % ----- adapt number of factors to samples ----%
    [ Factors,genetic_effects] = update_k( Factors,genetic_effects,b0,b1,i,epsilon,prop );

  % -- save sampled values (after thinning) -- %
    if mod(i,thin)==0 && i > burn
                
        sp_num = (i-burn)/thin;
        
        Posterior =  save_posterior_samples( sp_num,params, ...
            Posterior,resid,fixed_effects,genetic_effects,Factors,interaction_effects);
        
        if mod(sp_num,100)==0
            save('Posterior','Posterior','params')
        end

    end


    % -- provide run diagnostics and plots -- %
   if mod(i,draw_iter) == 0   
        directory=strread(pwd,'%s','delimiter','/');
        disp(directory(end))
        disp(i)
        Factors.nofout(i)-Factors.num
        elapsed=toc;
        %output some running statistics on the current factors and their
        %genetic variances
        [Factors.delta,[1:Factors.k]' Factors.h2,sum(Factors.scores'.^2)'./(size(Y,1)-1),sum(genetic_effects.U'.^2)'./(size(Z_1,1)-1)]
        disp(strcat('Time remaining:',num2str((nrun-i) * (elapsed/i) * 1/60)));

        %make some plots of some running statistics
        if simulation
           draw_simulation_diagnostics(i,sp_num,params,Factors,genetic_effects,resid,...
        Posterior,gen_factor_Lambda,error_factor_Lambda,G,R,h2)
          else
           draw_results_diagnostics(i,sp_num,params,Factors,Posterior)
        end
   end
end
toc
save('Posterior','Posterior','params')

end

