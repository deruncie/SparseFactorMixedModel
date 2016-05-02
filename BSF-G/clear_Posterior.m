function [Posterior,params] = clear_Posterior(Posterior,params)
%resets Posterior samples if burnin was not sufficient
    params.burn = params.burn + params.thin*size(Posterior.Lambda,2);
    p = size(Posterior.resid_Y_prec,1);
    b = size(Posterior.B,1);
    n = size(Posterior.W,1);
    r = size(Posterior.E_a,1);
    r2 = size(Posterior.W,1);
    
    Posterior.Lambda        = zeros(size(Posterior.Lambda,1),0);
    Posterior.F             = zeros(size(Posterior.F,1),0);
    Posterior.F_a           = zeros(size(Posterior.F_a,1),0);
    Posterior.delta         = zeros(size(Posterior.delta,1),0);
    Posterior.F_h2          = zeros(size(Posterior.F_h2,1),0);
    Posterior.resid_Y_prec  = zeros(p,0);
    Posterior.E_a_prec      = zeros(p,0);
    Posterior.W_prec        = zeros(p,0);
    Posterior.B             = zeros(b,p);
    Posterior.W             = zeros(r2,p);
    Posterior.E_a           = zeros(r,p);
end
