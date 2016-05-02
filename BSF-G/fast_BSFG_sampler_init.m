function [data_matrices,params,priors,current_state,Posterior,simulation] = fast_BSFG_sampler_init(priors,params)       
%   Preps everything for an analysis with the function fast_BSFG_sampler. 
%       Loads data from setup.mat
%       Decides if the data is from a simulation
%       Uses prior to find initial values for all parameters
%       Initializes Posterior struct
%       Saves run parameters
%
%   Input file setup.mat should contain:
%     Y      gene expression data data: n x p
%     X      fixed-effect design matrix (optional)
%     Z_1    random effects 1 incidence matrix n x r1
%     Z_2    random effects 2 incidence matrix n x r2 (optional)
%     A      kinship matrix of lines: r1 x r1
%     optional:
%         U_act
%         gen_factor_Lambda
%         error_factor_Lambda
%         h2_act
%         G
%         R
%         B
%         factor_h2s
%         name

% ----------------------- %
% ------read data-------- %
% ----------------------- %
    load('../setup.mat');

    %Determine if 'setup.mat' contains output of a simulation, based on if
    %known factor loadings are included. Affects plotting functions
    simulation = true;
    if isempty(who('gen_factor_Lambda'))
        simulation=false;
    end
    simulation
    if simulation,
        if(isempty(who('B_act')))
            B_act = zeros(size(X,2),size(Y,2));
        end
        params.U_act                = U_act;
        params.gen_factor_Lambda    = gen_factor_Lambda;
        params.Lambda               = error_factor_Lambda;
        params.h2                   = h2;
        params.G                    = G;
        params.R                    = R;
        params.B                    = B_act;
        params.factor_h2s           = factor_h2s;
        params.name                 = name;
    end

    %normalize Y to have zero mean and unit variances among observed values,
    %allowing for NaNs.
    [n,p] = size(Y);
    
    r     = size(Z_1,2);
    Mean_Y = zeros(1,size(Y,2));
    VY = zeros(1,size(Y,2));
    for j=1:p
        Mean_Y(j) = mean(Y(~isnan(Y(:,j)),j));
        VY(j) = var(Y(~isnan(Y(:,j)),j));
        if isnan(VY(j))
            VY(j)=1;
        end
    end
    % Don't remove the mean and standardize the variance if it's a
    % simulation because this makes comparing to the simulated values more
    % difficult. 
    % Do we want to do this for real data, or let the user do it?
    if simulation
        VY = ones(size(VY));
    else
        Y = bsxfun(@minus,Y,Mean_Y);                     
        Y = bsxfun(@times,Y,1./sqrt(VY));    
    end
    
    Y_full = Y;
    
        
    %determine if a design matrix (X) exists (is loaded from setup.mat). If
    %not, make a dummy X-matrix with no columns.
    if ~exist('X','var')
        X=zeros(n,0);
    end
    if size(X,1) ~= n
        X=zeros(n,0);
    end
    b = size(X,2);
    
    %Determine if a second random effects design matrix exists. If not, make a
    %dummy matrix
    if ~exist('Z_2','var')
        Z_2 = zeros(n,0);
    end
    if size(Z_2,1) ~= n
        Z_2 = zeros(n,0);
    end
    r2 = size(Z_2,2);
    
    data_matrices.Y_full    = Y_full;
    data_matrices.Z_1       = Z_1;
    data_matrices.Z_2       = Z_2;
    data_matrices.X         = X;
            
% ----------------------------- %
% -----Initialize variables---- %
% ----------------------------- % 

   % --- transcript-level model
    % p-vector of probe residual precisions. 
	%  Prior: Gamma distribution for each element
    %       shape = resid_Y_prec_shape
    %       rate = resid_Y_prec_rate
    resid_Y_prec_shape  = priors.resid_Y_prec_shape;
    resid_Y_prec_rate   = priors.resid_Y_prec_rate;
    resid_Y_prec        = gamrnd(resid_Y_prec_shape,1/resid_Y_prec_rate,p,1);   
   
    % Factors:
    %  initial number of factors
    k = priors.k_init;

    % Factor loading precisions (except column penalty tauh).
	 %  Prior: Gamma distribution for each element. 
     %       shape = Lambda_df/2
     %       rate = Lambda_df/2
     %    Marginilizes to t-distribution with Lambda_df degrees of freedom
     %    on each factor loading, conditional on tauh
    Lambda_df   = priors.Lambda_df;
    Lambda_prec = gamrnd(Lambda_df/2,2/Lambda_df,[p,k]);
    
    % Factor penalty. tauh(h) = \prod_{i=1}^h \delta_i
	 %  Prior: Gamma distribution for each element of delta
     %     delta_1:
     %       shape = delta_1_shape
     %       rate = delta_1_rate
     %     delta_2 ... delta_m:
     %       shape = delta_2_shape
     %       rate = delta_2_rate
    delta_1_shape  = priors.delta_1_shape;
    delta_1_rate   = priors.delta_1_rate;
    delta_2_shape  = priors.delta_2_shape;
    delta_2_rate   = priors.delta_2_rate;
    delta          = [gamrnd(delta_1_shape,1/delta_1_rate);gamrnd(delta_2_shape,1/delta_2_rate,[k-1,1])];
    tauh           = cumprod(delta);
    
    % Total Factor loading precisions Lambda_prec * tauh
    Plam = bsxfun(@times,Lambda_prec,tauh');
    
    % Lambda - factor loadings
     %   Prior: Normal distribution for each element.
     %       mu = 0
     %       sd = sqrt(1/Plam)
    Lambda = zeros(p,k) + randn(p,k).*reshape(sqrt(1./Plam),p,k);
    
    % g-vector of specific precisions of genetic effects. 
	%  Prior: Gamma distribution for each element
    %       shape = E_a_prec_shape
    %       rate = E_a_prec_rate
    E_a_prec_shape = priors.E_a_prec_shape;
    E_a_prec_rate  = priors.E_a_prec_rate;
    E_a_prec       = gamrnd(E_a_prec_shape,1/E_a_prec_rate,p,1);
    
    % Genetic effects not accounted for by factors.
    %   Prior: Normal distribution on each element.
    %       mean = 0
    %       sd = 1./sqrt(E_a_prec)' on each row
    E_a = bsxfun(@times,randn(r,p),1./sqrt(E_a_prec)');
        
    % Latent factor heritabilties. h2 can take h2_divisions values
    %   between 0 and 1.
     %   Prior: 0.5: h2=0, .05: h2 > 0. 
    F_h2 = rand(k,1);
    
    % Genetic effects on the factor scores.
	%  Prior: Normal distribution for each element
    %       mean = 0
    %       sd = sqrt(F_h2') for each row.
    F_a = bsxfun(@times,randn(r,k),sqrt(F_h2'));
    
    % Full Factor scores. Combination of genetic and residual variation on
    % each factor.
	%  Prior: Normal distribution for each element
    %       mean = Z_1 * F_a
    %       sd = sqrt(1-F_h2') for each row.
    F = Z_1 * F_a + bsxfun(@times,randn(n,k),sqrt(1-F_h2'));
    
    
    % g-vector of specific precisions of genetic effects. 
	%  Prior: Gamma distribution for each element
    %       shape = E_a_prec_shape
    %       rate = E_a_prec_rate
    W_prec_shape = priors.W_prec_shape;
    W_prec_rate  = priors.W_prec_rate;
    W_prec       = gamrnd(W_prec_shape,1/W_prec_rate,p,1);
    
    % Genetic effects not accounted for by factors.
    %   Prior: Normal distribution on each element.
    %       mean = 0
    %       sd = 1./sqrt(E_a_prec)' on each row
    W = bsxfun(@times,randn(r2,p),1./sqrt(W_prec)');

    % Fixed effect coefficients.
	%  Prior: Normal distribution for each element
    %       mean = 0
    %       sd = sqrt(1/fixed_effects_prec)
    B                 = randn(b,p);

    
% ----------------------- %
% -Initialize Posterior-- %
% ----------------------- %

    Posterior.Lambda        = zeros(0,0);
    Posterior.F_a           = zeros(0,0);
    Posterior.F             = zeros(0,0);
    Posterior.delta         = zeros(0,0);
    Posterior.F_h2          = zeros(0,0);
    Posterior.resid_Y_prec  = zeros(p,0);
    Posterior.E_a_prec      = zeros(p,0);
    Posterior.W_prec        = zeros(p,0);
    Posterior.B             = zeros(b,p);
    Posterior.W             = zeros(r2,p);
    Posterior.E_a           = zeros(r,p);

% ----------------------- %
% ---Save initial values- %
% ----------------------- %

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


% ------------------------------------ %
% ----Precalculate some matrices------ %
% ------------------------------------ %

    %invert the random effect covariance matrices
    Ainv = inv(A);
    A_2_inv = eye(size(Z_2,2)); %Z_2 random effects are assumed to have covariance proportional to the identity. Can be modified.

    %pre-calculate transformation parameters to diagonalize aI + bZAZ for fast
    %inversion: inv(aI + bZAZ) = 1/b*U*diag(1./(s+a/b))*U'
    %uses singular value decomposition of ZAZ for stability when ZAZ is low
    %rank
%     XZ = [X_f Z_1];
%     [U,S,~]          = svd(XZ*blkdiag(1e6*eye(b_f),A)*XZ');
    [U,S,~]          = svd(Z_1*A*Z_1');
    invert_aI_bZAZ.U = U;
    invert_aI_bZAZ.s = diag(S);   
    
   
    %fixed effects + random effects 1
    %diagonalize mixed model equations for fast inversion: 
    %inv(a*blkdiag(fixed_effect_prec*eye(b),Ainv) + b*[X Z_1]'[X Z_1]) = U*diag(1./(a.*s1+b.*s2))*U'
    Design=[X Z_1];
    Design2 = Design'*Design;
%     [~,~,U,S1,S2] = gsvd(cholcov(blkdiag(fixed_effect_prec*eye(b),Ainv)),cholcov(Design2));
    [~,~,U,S1,S2] = gsvd(cholcov(blkdiag(zeros(b,b),Ainv)),cholcov(Design2));
    invert_aPXA_bDesignDesignT.U = inv(U)';
    invert_aPXA_bDesignDesignT.s1 = diag(S1'*S1);
    invert_aPXA_bDesignDesignT.s2 = diag(S2'*S2);
    invert_aPXA_bDesignDesignT.Design_U = Design*invert_aPXA_bDesignDesignT.U; 
        
    %random effects 2
    %diagonalize mixed model equations for fast inversion: 
    %inv(a*A_2_inv + b*Z_2'Z_2]) = U*diag(1./(a.*s1+b.*s2))*U'
    Design=Z_2;
    Design2 = Design'*Design;
    [~,~,U,S1,S2] = gsvd(cholcov(A_2_inv),cholcov(Design2));
    invert_aPXA_bDesignDesignT_rand2.U = inv(U)';
    invert_aPXA_bDesignDesignT_rand2.s1 = diag(S1'*S1);
    invert_aPXA_bDesignDesignT_rand2.s2 = diag(S2'*S2);
    invert_aPXA_bDesignDesignT_rand2.Design_U = Design*invert_aPXA_bDesignDesignT_rand2.U;  
    
    %genetic effect variances of factor traits
    % diagonalizing a*Z_1'*Z_1 + b*Ainv for fast inversion
    %diagonalize mixed model equations for fast inversion: 
    % inv(a*Z_1'*Z_1 + b*Ainv) = U*diag(1./(a.*s1+b.*s2))*U'
    %similar to fixed effects + random effects 1 above, but no fixed effects.
    ZZt = Z_1'*Z_1;
    [~,~,U,S1,S2] = gsvd(cholcov(ZZt),cholcov(Ainv));
    invert_aZZt_Ainv.U = inv(U)';
    invert_aZZt_Ainv.s1 = diag(S1'*S1);
    invert_aZZt_Ainv.s2 = diag(S2'*S2);


% ----------------------------- %
% ----Save run parameters------ %
% ----------------------------- %

    params.p                = p;
    params.n                = n;
    params.r                = r;
    params.r2               = r2;
    params.Mean_Y           = Mean_Y;
    params.VY               = VY;
    params.burn             = 0;
    params.thin             = 0;
    params.sp               = 0;
    params.nrun             = 0;
    params.Ainv             = Ainv;
    params.A_2_inv          = A_2_inv;
    params.invert_aI_bZAZ               = invert_aI_bZAZ;  
    params.invert_aPXA_bDesignDesignT   = invert_aPXA_bDesignDesignT;       
    params.invert_aZZt_Ainv             = invert_aZZt_Ainv;
    params.invert_aPXA_bDesignDesignT_rand2   = invert_aPXA_bDesignDesignT_rand2; 

end
