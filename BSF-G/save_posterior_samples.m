function [ Posterior ] = save_posterior_samples( sp_num,~, ...
      Posterior,Lambda,F,F_a,B,W,E_a,delta,F_h2,resid_Y_prec,E_a_prec,W_prec)
%save posteriors. Full traces are kept of the more interesting parameters.
%Only the posterior means are kept of less interesting parameters. These
%should be correctly calculated over several re-starts of the sampler.

sp = size(Posterior.Lambda,2);

%save factor samples
if numel(Lambda) > size(Posterior.Lambda,1),
    %expand factor sample matrix if necessary
    Posterior.Lambda = [Posterior.Lambda; zeros(numel(Lambda)-size(Posterior.Lambda,1),sp)];
    Posterior.F   = [Posterior.F; zeros(numel(F)-size(Posterior.F,1),sp)];
    Posterior.F_a = [Posterior.F_a; zeros(numel(F_a)-size(Posterior.F_a,1),sp)];
    Posterior.delta = [Posterior.delta; zeros(numel(delta)-size(Posterior.delta,1),sp)];
    Posterior.F_h2 = [Posterior.F_h2; zeros(numel(F_h2)-size(Posterior.F_h2,1),sp)];
end                
Posterior.Lambda(1:numel(Lambda),sp_num) = Lambda(:);
Posterior.F(1:numel(F),sp_num)     = F(:);
Posterior.F_a(1:numel(F_a),sp_num) = F_a(:);
Posterior.delta(1:numel(delta),sp_num) = delta;
Posterior.F_h2(1:numel(F_h2),sp_num) = F_h2;

Posterior.resid_Y_prec(:,sp_num) = resid_Y_prec;
Posterior.E_a_prec(:,sp_num) = E_a_prec;
Posterior.W_prec(:,sp_num) = W_prec;

%save B,U,W
Posterior.B = (Posterior.B*(sp_num-1) + B)./sp_num;
Posterior.E_a = (Posterior.E_a*(sp_num-1) + E_a)./sp_num;
Posterior.W = (Posterior.W*(sp_num-1) + W)./sp_num;

end

