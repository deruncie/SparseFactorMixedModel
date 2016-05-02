function [ Posterior ] =save_posterior_samples( sp_num,params, ...
      Posterior,resid,fixed_effects,genetic_effects,Factors,interaction_effects)
%save posteriors. Re-scale samples back to original variances.
  
%save factors
sp = params.sp;
VY = params.VY;       
Lambda = bsxfun(@times,Factors.Lambda,sqrt(VY'));     %re-scale by Y variances
G_h2 = Factors.h2;
U = genetic_effects.U;     
delta = Factors.delta;
genetic_ps = genetic_effects.ps./VY';
resid_ps = resid.ps./VY';

%save factor samples
Lambda = Lambda(:,1:Factors.k);
if numel(Lambda) > size(Posterior.Lambda,1),
    %expand factor sample matrix if necessary
    Posterior.Lambda = [Posterior.Lambda; zeros(numel(Lambda)-size(Posterior.Lambda,1),sp)];
    Posterior.U = [Posterior.U; zeros(numel(U)-size(Posterior.U,1),sp)];
    Posterior.delta = [Posterior.delta; zeros(numel(delta)-size(Posterior.delta,1),sp)];
    Posterior.G_h2 = [Posterior.G_h2; zeros(numel(G_h2)-size(Posterior.G_h2,1),sp)];
end                
Posterior.Lambda(1:numel(Lambda),sp_num) = Lambda(:);
Posterior.delta(1:numel(delta),sp_num) = delta;
Posterior.G_h2(1:numel(G_h2),sp_num) = G_h2;
Posterior.U(1:numel(U),sp_num) = U(:);

Posterior.ps(:,sp_num) = genetic_ps;
Posterior.resid_ps(:,sp_num) = resid_ps;

%save B,U,W
Posterior.B = (Posterior.B*(sp_num-1) + bsxfun(@times,fixed_effects.B,sqrt(VY')))./sp_num;
Posterior.d = (Posterior.d*(sp_num-1) + bsxfun(@times,genetic_effects.d,sqrt(VY')))./sp_num;
Posterior.W = (Posterior.W*(sp_num-1) + bsxfun(@times,interaction_effects.W,sqrt(VY')))./sp_num;

if mod(sp_num,100)==0
    save('Posterior','Posterior','params')
end

end

