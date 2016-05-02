function [Posterior,current_state,params,i] = fast_BSFG_sampler(data_matrices,...
                                    start_i,draw_iter,burn,sp,thin,simulation, ...
                                    params,priors,Posterior,current_state)
% -- Daniel Runcie -- %
%
% Gibbs sampler for genetic covariance estimation based on mixed effects
% model, with missing data
%
% Code for:
% 
% Runcie and Mukherjee (2013) Dissecting high-dimensional traits
% with Bayesian sparse factor analysis of genetic covariance matrices.
% GENETICS.
%
% (c) July 30, 2013
%
% code based on original provided by Anirban Bhattacharya
%
%         
% This function implements the BSF-G partially collapsed Gibbs sampler.
% Variable notation follows the Runcie and Mukherjee manuscript as closely
% as possible.
%
% All input and initialization functions are carried out by the function
% MA_sampler_init. See the documentation of that function for details.
% 
% The sampler is designed to do a short-medium length run and then return
% the state of the chain. After assessing the progress of the chain,
% (optionally), the Posterior matrices can be reset, and the chain
% re-started where it left off. Right now, the random seed is not saved. Probably should
% be if I can figure out how.
% 
% This function takes the following inputs:
%     data_matrices: struct holding the (should be imutable) data, design and incidence matrices:
%           Y_full: Full probe data, may include NaNs. n x p
%           X      fixed-effect design matrix (optional)
%           Z_1    random effects 1 incidence matrix n x r1
%           Z_2    random effects 2 incidence matrix n x r2 (optional)
%     start_i: iteration number of the end of the last run.
%     draw_iter: frequency of updating diagnostic plots
%     burn: number of burnin samples
%     sp: total number of samples to collect
%     thin: thinning rate of chain
%     simulation: boolean. Is this a simulation?
%     params: struct with chain parameters.
%     priors: struct with all relevant prior hyperparameters
%     Posterior: struct with posterior matrices, or running posterior means. 
%            Note: not sure how well posterior means work after re-starting chain. Probably not well.
%     current_state: current (initial) conditions of all model parameters
% 
% Several diagnostic plots are produced during the run. 
%     Their interpretation is described within the source codes:
%         draw_simulation_diagnostics.m: For simulated data with known true values
%         draw_results_diagnostics.m: Otherwise
%         

% ----------------------------------------------- %
% ----------------Load data matrices------------- %
% ----------------------------------------------- %

Y_full  = data_matrices.Y_full;
Z_1     = data_matrices.Z_1;
Z_2     = data_matrices.Z_2;
X       = data_matrices.X;

% ----------------------------------------------- %
% ----------------Load priors-------------------- %
% ----------------------------------------------- %

resid_Y_prec_shape  =   priors.resid_Y_prec_shape;
resid_Y_prec_rate   =   priors.resid_Y_prec_rate;
E_a_prec_shape      =   priors.E_a_prec_shape;
E_a_prec_rate       =   priors.E_a_prec_rate;
W_prec_shape        =   priors.W_prec_shape;
W_prec_rate         =   priors.W_prec_rate;
Lambda_df           =   priors.Lambda_df;
delta_1_shape       =   priors.delta_1_shape;
delta_1_rate        =   priors.delta_1_rate;
delta_2_shape       =   priors.delta_2_shape;
delta_2_rate        =   priors.delta_2_rate;
h2_priors_factors   =   priors.h2_priors_factors;

% ----------------------------------------------- %
% ----------------Load current state------------- %
% ----------------------------------------------- %
resid_Y_prec        =   current_state.resid_Y_prec;
F                   =   current_state.F;
Lambda              =   current_state.Lambda;
E_a                 =   current_state.E_a;
W                   =   current_state.W;
Lambda_prec         =   current_state.Lambda_prec;
delta               =   current_state.delta;
tauh                =   current_state.tauh;
E_a_prec            =   current_state.E_a_prec;
Plam                =   current_state.Plam;
F_h2                =   current_state.F_h2;
F_a                 =   current_state.F_a;
B                   =   current_state.B;

% ----------------------------------------------- %
% -----------Reset Global Random Number Stream--- %
% ----------------------------------------------- %
RandStream.setGlobalStream(current_state.randStream);
stream = RandStream.getGlobalStream();
stream.State = current_state.rand_state;

% ----------------------------------------------- %
% -----------Load pre-calcualted matrices-------- %
% ----------------------------------------------- %
h2_divisions                        =   params.h2_divisions;
invert_aI_bZAZ                      =   params.invert_aI_bZAZ;   
invert_aPXA_bDesignDesignT          =   params.invert_aPXA_bDesignDesignT;  
invert_aZZt_Ainv                    =   params.invert_aZZt_Ainv;
invert_aPXA_bDesignDesignT_rand2    =   params.invert_aPXA_bDesignDesignT_rand2;   
Ainv                                =   params.Ainv;
A_2_inv                             =   params.A_2_inv;


% ----------------------------------------------- %
% ----------------Set up run--------------------- %
% ----------------------------------------------- %
%     b0,b1: parameters controlling rate of adaptation of factor model size
%     h2_divisions: number of discrete steps for each factor heritability parameter
%     epsilon: truncation point for factor loadings during adaptation
nrun        = burn+sp*thin; % number of posterior samples
b0          = params.b0;b1=params.b1;
epsilon     = params.epsilon;
prop        = params.prop;
save_freq   = params.save_freq;


% ----------------------------------------------- %
% -----------Update run length parameters-------- %
% ----------------------------------------------- %

params.burn = params.burn+burn;
params.thin = thin;
params.sp   = params.sp+sp;
params.nrun = params.nrun+nrun;


% ----------------------------------------------- %
% ---Extend posterior matrices for new samples--- %
% ----------------------------------------------- %

Posterior.Lambda        = [Posterior.Lambda,zeros(size(Posterior.Lambda,1),sp)];
Posterior.F             = [Posterior.F,zeros(size(Posterior.F,1),sp)];
Posterior.F_a           = [Posterior.F_a,zeros(size(Posterior.F_a,1),sp)];
Posterior.delta         = [Posterior.delta,zeros(size(Posterior.delta,1),sp)];
Posterior.F_h2          = [Posterior.F_h2,zeros(size(Posterior.F_h2,1),sp)];
Posterior.resid_Y_prec  = [Posterior.resid_Y_prec,zeros(size(Posterior.resid_Y_prec,1),sp)];
Posterior.E_a_prec      = [Posterior.E_a_prec,zeros(size(Posterior.E_a_prec,1),sp)];
Posterior.W_prec        = [Posterior.W_prec,zeros(size(Posterior.W_prec,1),sp)];

% ----------------------------------------------- %
% --------------start gibbs sampling------------- %
% ----------------------------------------------- %

sp_num = 0;
tic
for i = start_i+(1:nrun)
   
 % -----fill in missing phenotypes----- %
    %conditioning on everything else
    Y = Y_full;
    phenMissing = isnan(Y_full);
    if sum(sum(phenMissing))>0
        meanTraits = X*B + F*Lambda' + Z_1*E_a + Z_2*W;
        resids = bsxfun(@times,randn(size(Y_full)),1./sqrt(resid_Y_prec'));
        Y(phenMissing) = meanTraits(phenMissing) + resids(phenMissing);
    end
    
    
% -----Sample Lambda------------------ %
    %conditioning on W, B, F, marginalizing over E_a
    Y_tilde = Y - X*B - Z_2*W;
    Lambda = sample_Lambda( Y_tilde,F,resid_Y_prec, E_a_prec,Plam,invert_aI_bZAZ );
   
    
 % -----Sample B and E_a--------------- %
    %conditioning on W, F, Lambda
    Y_tilde = Y - F*Lambda' - Z_2*W;
    b = size(B,1);
    r = size(E_a,1);
    [location_sample] = sample_means( Y_tilde, ...
            resid_Y_prec, E_a_prec, invert_aPXA_bDesignDesignT );
    B   = location_sample(:,1:b)';
    E_a = location_sample(:,b+(1:r))';

    
 % -----Sample W ---------------------- %
    %conditioning on B, E_a, F, Lambda
    Y_tilde = Y - X*B - Z_1*E_a - F*Lambda';
    if size(Z_2,2) > 0,
        [location_sample] = sample_means( Y_tilde, ...
            resid_Y_prec, W_prec, invert_aPXA_bDesignDesignT_rand2 );
        W = location_sample';
    end    
  
    
 % -----Sample F_h2-------------------- %
    %conditioning on F, marginalizing over F_a
    F_h2 = sample_h2s_discrete(F,h2_divisions,h2_priors_factors,invert_aI_bZAZ);
 
    
 % -----Sample F_a--------------------- %
    %conditioning on F, F_h2
    F_a = sample_F_a(F,Z_1,F_h2,invert_aZZt_Ainv);
   
    
 % -----Sample F----------------------- %
    %conditioning on B,F_a,E_a,W,Lambda, F_h2
    Y_tilde = Y - X*B - Z_1*E_a - Z_2*W;
    F = sample_factors_scores( Y_tilde, Z_1,Lambda,resid_Y_prec,F_a,F_h2 );

    
 % -----Sample Lambda_prec------------- %
    Lambda2 = Lambda.^2;    
    Lambda_prec = gamrnd(Lambda_df/2 + 0.5,2./(Lambda_df + bsxfun(@times,Lambda2,tauh')));
    
    
 % -----Sample resid_Y_prec------------ %
    Y_tilde = Y - X*B - F*Lambda' - Z_1*E_a - Z_2*W;
    resid_Y_prec = gamrnd(resid_Y_prec_shape + 0.5*size(Y_tilde,1),1./(resid_Y_prec_rate+0.5*sum(Y_tilde.^2,1)'));  %model residual precision
    
    
 % -----Sample E_a_prec---------------- %
    E_a_prec = gamrnd(E_a_prec_shape + 0.5*size(E_a,1),1./(E_a_prec_rate+0.5*diag(E_a' * Ainv * E_a))); %random effect 1 (D) residual precision
    
    
 % -----Sample W_prec------------------ %
    W_prec = gamrnd(W_prec_shape + 0.5*size(E_a,1),1./(W_prec_rate+0.5*diag(W' * A_2_inv * W))); %random effect 1 (D) residual precision
 
    
 % -----Sample delta, update tauh------ %
    [delta,tauh] = sample_delta( delta,tauh,Lambda_prec,delta_1_shape,delta_1_rate,delta_2_shape,delta_2_rate,Lambda2 );
    
    
 % -----Update Plam-------------------- %
    Plam = bsxfun(@times,Lambda_prec,tauh');    

    
 % -- adapt number of factors to samples ---%
    [F,Lambda,F_a,F_h2,Lambda_prec,Plam,delta,tauh] = update_k( F,Lambda,F_a,F_h2,Lambda_prec,Plam,delta,tauh,Z_1,Lambda_df,delta_2_shape,delta_2_rate,b0,b1,i,epsilon,prop );

    
 % -- save sampled values (after thinning) -- %
    if mod(i-params.burn,thin)==0 && i > params.burn
            
        sp_num = (i-params.burn)/thin;
        
        Posterior = save_posterior_samples( sp_num,params, ...
            Posterior,Lambda,F,F_a,B,W,E_a,delta,F_h2,resid_Y_prec,E_a_prec,W_prec);
        
        if mod(sp_num,save_freq)==0
            save('Posterior','Posterior','params')
        end

    end


 % -- provide run diagnostics and plots -- %
   if mod(i,draw_iter) == 0   
        directory=strread(pwd,'%s','delimiter','/');
        disp(directory(end))
        disp(i)
        elapsed=toc;
        %output some running statistics on the current factors and their genetic variances
        % each row represents a factor:
        %    delta[i], [i], F_h2[i], V_P(F_i), V_A(F_i)
        [delta,[1:size(F,2)]',F_h2,sum(F.^2)'./(size(F,1)-1),sum(F_a.^2)'./(size(F_a,1)-1)]
        disp(strcat('Time remaining:',num2str((nrun-i) * (elapsed/i) * 1/60)));

        %make some plots of some running statistics
        if simulation
           gen_factor_Lambda    = params.gen_factor_Lambda;
           error_factor_Lambda  = params.Lambda;
           G_act                = params.G;
           E_act                = params.R;
           h2                   = params.h2;
           draw_simulation_diagnostics(sp_num,params,Posterior,Lambda,F_h2,E_a_prec,resid_Y_prec,gen_factor_Lambda,error_factor_Lambda,G_act,E_act,h2)
           % optional to also draw diagnostics for fixed effects. 10: fixed
           % effects on factors; 11: fixed effects on genes.
%            figure(10);plot(X_f*current_state.F_b*current_state.Lambda',X_f*params.B_f*params.Lambda','.');line(xlim,xlim)
%            figure(11);plot(X*current_state.B,X*params.B,'.');line(xlim,xlim)
          else
           draw_results_diagnostics(sp_num,params,Lambda, F_h2, Posterior)
        end
   end
end
toc
save('Posterior','Posterior','params')


% ----------------------------------------------- %
% ------------Save state for restart------------- %
% ----------------------------------------------- %

current_state.resid_Y_prec  =   resid_Y_prec;
current_state.Lambda_prec   =   Lambda_prec;
current_state.delta         =   delta;
current_state.tauh          =   tauh;
current_state.Plam          =   Plam;
current_state.Lambda        =   Lambda;
current_state.F_h2          =   F_h2;
current_state.E_a_prec      =   E_a_prec;
current_state.W_prec        =   W_prec;
current_state.F_a           =   F_a;
current_state.F             =   F;
current_state.E_a           =   E_a;
current_state.B             =   B;
current_state.W             =   W;
current_state.rand_state    =   RandStream.getGlobalStream.State;
save('current_state','current_state','params')
end

